#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "CandidateDBReader.h"
#include "CandidateDBWriter.h"
#include "DBReader.h"

#include <algorithm>
#include <atomic>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#ifdef OPENMP
#include <omp.h>
#endif

namespace {

constexpr size_t COVERAGE_BIN_COUNT = 65536;
constexpr size_t MAX_REPORT_ROWS = 50;

// Load speciesId -> genome size from <dbDir>/species2genomeSize.tsv
// (col0 = speciesId, col2 = genome size), matching Classifier::parseSp2GenomeSize.
std::unordered_map<TaxID, uint64_t> loadGenomeSizes(const std::string &dbDir) {
    std::unordered_map<TaxID, uint64_t> sp2genomeSize;
    const std::string fileName = dbDir + "/species2genomeSize.tsv";
    std::ifstream infile(fileName);
    if (!infile.is_open()) {
        std::cout << "Warning: could not open " << fileName
                  << "; coverage-based filtering will be skipped." << std::endl;
        return sp2genomeSize;
    }
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        const size_t tab1 = line.find('\t');
        if (tab1 == std::string::npos) continue;
        const size_t tab2 = line.find('\t', tab1 + 1);
        if (tab2 == std::string::npos) continue;

        // Exception-free parsing (the project is built with -fno-exceptions).
        const std::string spIdStr = line.substr(0, tab1);
        const std::string sizeStr = line.substr(tab2 + 1);
        char *endPtr = nullptr;
        const unsigned long long spId = std::strtoull(spIdStr.c_str(), &endPtr, 10);
        if (endPtr == spIdStr.c_str()) continue; // no parseable species id
        const unsigned long long genomeSize = std::strtoull(sizeStr.c_str(), &endPtr, 10);
        if (endPtr == sizeStr.c_str()) continue; // no parseable genome size
        sp2genomeSize[static_cast<TaxID>(spId)] = static_cast<uint64_t>(genomeSize);
    }
    return sp2genomeSize;
}

// -------- Filter method 0: average score + read count + genome coverage (adjusted evenness) --------
//
// Score and count use only each read's TOP hit (the species the read is assigned to),
// so they match report.tsv's per-species avg_score rather than averaging over all
// candidate occurrences (including reads the species only ranks second on). Removes a
// candidate species if its mean top-hit idScore is below --min-avg-score, if it is the
// top hit for fewer than --min-count reads, or (when the DB stores k-mer positions and a
// genome size is known) if its per-species genome coverage has an adjusted evenness
// below --min-adj-evenness.
std::unordered_set<TaxID> filterByScoreAndCoverage(
    CandidateDBReader &reader,
    size_t entryCount,
    const std::unordered_map<TaxID, uint64_t> &sp2genomeSize,
    const LocalParameters &par)
{
    const float minAvgScore = par.minAvgScore;
    const float minAdjEvenness = par.minAdjEvenness;
    const uint64_t minCount = par.minCount < 0 ? 0 : static_cast<uint64_t>(par.minCount);
    const bool useAllHits = par.covUseAllHits != 0;

    std::cout << "Filter method          : top-hit score + read count + genome coverage" << std::endl;
    std::cout << "Min. average score     : " << minAvgScore
              << (minAvgScore <= 0.0f ? " (score filter disabled)" : "") << std::endl;
    std::cout << "Min. top-hit reads     : " << minCount
              << (minCount == 0 ? " (count filter disabled)" : "") << std::endl;
    std::cout << "Min. adjusted evenness : " << minAdjEvenness
              << (minAdjEvenness <= 0.0f ? " (coverage filter disabled)" : "") << std::endl;
    std::cout << "Coverage aggregation   : " << (useAllHits ? "all candidate hits" : "top hit per read only") << std::endl;

    // --- Pass 1: accumulate per-species scores and (when available) coverage ---
    // Single-threaded: coverage keeps a 64 Ki-bin histogram per species, so a
    // shared accumulator avoids large per-thread duplication.
    std::unordered_map<TaxID, double> speciesScoreSum;
    std::unordered_map<TaxID, uint64_t> speciesScoreCount;
    std::unordered_map<TaxID, std::vector<uint8_t>> sp2bins;
    std::unordered_map<TaxID, uint64_t> sp2readCnt;
    std::unordered_map<TaxID, uint64_t> sp2readLen;

    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (!reader.getByIndex(i, entry, 0)) {
            continue;
        }
        for (size_t idx = 0; idx < entry.candidates.size(); ++idx) {
            const SpeciesCandidate &candidate = entry.candidates[idx];
            // Average score and read count use only the read's top hit (index 0),
            // i.e. the species this read is assigned to. This mirrors report.tsv's
            // per-species avg_score and excludes reads the species only ranks
            // behind a better match on.
            if (idx == 0) {
                speciesScoreSum[candidate.speciesId] += candidate.idScore;
                speciesScoreCount[candidate.speciesId] += 1;
            }

            // Coverage uses all hits, or only the top-scoring hit (index 0),
            // and only when this candidate carries k-mer positions.
            if ((useAllHits || idx == 0) && !candidate.kmerPositions.empty()) {
                std::vector<uint8_t> &bins = sp2bins[candidate.speciesId];
                if (bins.empty()) {
                    bins.resize(COVERAGE_BIN_COUNT, 0);
                }
                for (const uint16_t pos : candidate.kmerPositions) {
                    if (bins[pos] < 255) {
                        bins[pos]++;
                    }
                }
                sp2readCnt[candidate.speciesId] += 1;
                sp2readLen[candidate.speciesId] += entry.queryLength;
            }
        }
    }

    // --- Decide which species survive ---
    // Only species that are a read's top hit at least once appear in
    // speciesScoreCount; a species that is never a top hit wins no reads and is
    // dropped implicitly (it is never added to keptSpecies).
    struct SpeciesDecision {
        TaxID spId;
        uint64_t topHitCount;   // reads where this species is the top hit
        double avgScore;
        double coverage;
        double adjEvenness;
        bool coverageEvaluated;
        bool removeScore;
        bool removeCount;
        bool removeEvenness;
    };

    std::unordered_set<TaxID> keptSpecies;
    keptSpecies.reserve(speciesScoreCount.size());
    std::vector<SpeciesDecision> decisions;
    decisions.reserve(speciesScoreCount.size());
    size_t removedByScore = 0;
    size_t removedByCount = 0;
    size_t removedByEvenness = 0;
    size_t coverageEvaluated = 0;
    size_t coverageMissingGenomeSize = 0;

    for (const auto &kv : speciesScoreCount) {
        const TaxID spId = kv.first;
        const uint64_t topHitCount = kv.second;
        const double avgScore = speciesScoreSum[spId] / static_cast<double>(topHitCount);

        bool removeScore = (minAvgScore > 0.0f) && (avgScore < static_cast<double>(minAvgScore));
        bool removeCount = (minCount > 0) && (topHitCount < minCount);

        bool removeEvenness = false;
        bool covEval = false;
        double coverage = 0.0;
        double adjEvenness = 0.0;
        const auto binIt = sp2bins.find(spId);
        if (binIt != sp2bins.end()) {
            const auto gsIt = sp2genomeSize.find(spId);
            if (gsIt != sp2genomeSize.end() && gsIt->second > 0) {
                const CovMetric metric = computeCoverageMetric(
                    binIt->second, sp2readCnt[spId], sp2readLen[spId], gsIt->second);
                covEval = true;
                coverage = metric.coverage;
                adjEvenness = metric.adjustedEvenness;
                ++coverageEvaluated;
                if (minAdjEvenness > 0.0f && metric.adjustedEvenness < static_cast<double>(minAdjEvenness)) {
                    removeEvenness = true;
                }
            } else {
                ++coverageMissingGenomeSize;
            }
        }

        // Report each species against the first criterion it fails (score, then
        // count, then evenness); a species is kept only if it passes all three.
        if (removeScore) ++removedByScore;
        else if (removeCount) ++removedByCount;
        else if (removeEvenness) ++removedByEvenness;

        if (!removeScore && !removeCount && !removeEvenness) {
            keptSpecies.insert(spId);
        }
        decisions.push_back({spId, topHitCount, avgScore, coverage, adjEvenness,
                             covEval, removeScore, removeCount, removeEvenness});
    }

    std::cout << "Distinct top-hit species   : " << speciesScoreCount.size() << std::endl;
    std::cout << "Coverage evaluated for : " << coverageEvaluated << " species" << std::endl;
    if (coverageMissingGenomeSize > 0) {
        std::cout << "Coverage skipped (no genome size) : " << coverageMissingGenomeSize << " species" << std::endl;
    }
    std::cout << "Species removed (score)    : " << removedByScore << std::endl;
    std::cout << "Species removed (count)    : " << removedByCount << std::endl;
    std::cout << "Species removed (evenness) : " << removedByEvenness << std::endl;
    std::cout << "Species kept               : " << keptSpecies.size() << std::endl;

    // Per-species report (capped to keep the log readable).
    std::sort(decisions.begin(), decisions.end(),
              [](const SpeciesDecision &a, const SpeciesDecision &b) {
                  return a.topHitCount > b.topHitCount;
              });
    std::cout << "species\ttopHitReads\tavgScore\tcoverage\tadjEvenness\tdecision" << std::endl;
    for (size_t r = 0; r < decisions.size() && r < MAX_REPORT_ROWS; ++r) {
        const SpeciesDecision &d = decisions[r];
        std::cout << d.spId << '\t' << d.topHitCount << '\t' << d.avgScore << '\t';
        if (d.coverageEvaluated) {
            std::cout << d.coverage << '\t' << d.adjEvenness << '\t';
        } else {
            std::cout << "-\t-\t";
        }
        if (d.removeScore) std::cout << "removed(score)";
        else if (d.removeCount) std::cout << "removed(count)";
        else if (d.removeEvenness) std::cout << "removed(evenness)";
        else std::cout << "kept";
        std::cout << std::endl;
    }
    if (decisions.size() > MAX_REPORT_ROWS) {
        std::cout << "... (" << (decisions.size() - MAX_REPORT_ROWS) << " more species)" << std::endl;
    }

    return keptSpecies;
}

// -------- Filter method 1: best-evidence (A) + uniqueness (B), combined conservatively (F) --------
//
// A (best-evidence): a real species has some reads that match it well even if the true
//   species is absent, whereas a spurious species never does. Count reads whose idScore
//   is >= --min-strong-score; require >= --min-strong-reads of them.
// B (uniqueness): a homology "hitchhiker" is never the sole best explanation of a read.
//   Count reads where the species is the unique top candidate (its score clearly beats
//   the runner-up, i.e. runner-up < top * --tie-ratio); require >= --min-unique-reads.
// F (conservative): keep the species if it passes EITHER criterion; remove only when it
//   fails BOTH, so the filter errs toward keeping borderline-real species.
std::unordered_set<TaxID> filterByEvidenceAndUniqueness(
    CandidateDBReader &reader,
    size_t entryCount,
    const LocalParameters &par)
{
    const float minStrongScore = par.minStrongScore;
    const uint64_t minStrongReads = par.minStrongReads < 0 ? 0 : static_cast<uint64_t>(par.minStrongReads);
    const uint64_t minUniqueReads = par.minUniqueReads < 0 ? 0 : static_cast<uint64_t>(par.minUniqueReads);
    const float tieRatio = par.tieRatio;

    std::cout << "Filter method          : best-evidence + uniqueness" << std::endl;
    std::cout << "Min. strong score      : " << minStrongScore << std::endl;
    std::cout << "Min. strong reads      : " << minStrongReads << std::endl;
    std::cout << "Min. unique-top reads  : " << minUniqueReads << std::endl;
    std::cout << "Tie ratio (uniqueness) : " << tieRatio << std::endl;

    struct EvidenceStats {
        uint64_t occurrences = 0;   // number of reads this species is a candidate for
        uint64_t strongReads = 0;   // reads with idScore >= minStrongScore
        uint64_t uniqueTop = 0;     // reads where it is the unique top candidate
        float maxScore = 0.0f;      // best idScore observed
    };
    std::unordered_map<TaxID, EvidenceStats> stats;

    // --- Pass 1: per-species evidence and uniqueness counts ---
    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (!reader.getByIndex(i, entry, 0)) {
            continue;
        }
        const std::vector<SpeciesCandidate> &cands = entry.candidates;
        if (cands.empty()) {
            continue;
        }

        // Candidates are stored best-first: the read's top candidate is unique when
        // the runner-up is clearly below it (not within the tie margin).
        const float topScore = cands[0].idScore;
        const float secondScore = (cands.size() > 1) ? cands[1].idScore : 0.0f;
        const bool topIsUnique = (cands.size() == 1) || (secondScore < topScore * tieRatio);

        for (size_t idx = 0; idx < cands.size(); ++idx) {
            const SpeciesCandidate &candidate = cands[idx];
            EvidenceStats &st = stats[candidate.speciesId];
            st.occurrences += 1;
            if (candidate.idScore > st.maxScore) {
                st.maxScore = candidate.idScore;
            }
            if (candidate.idScore >= minStrongScore) {
                st.strongReads += 1;
            }
            if (idx == 0 && topIsUnique) {
                st.uniqueTop += 1;
            }
        }
    }

    // --- Decide which species survive (keep if it passes EITHER criterion) ---
    struct SpeciesDecision {
        TaxID spId;
        uint64_t occurrences;
        uint64_t strongReads;
        uint64_t uniqueTop;
        float maxScore;
        bool passEvidence;
        bool passUnique;
    };

    std::unordered_set<TaxID> keptSpecies;
    keptSpecies.reserve(stats.size());
    std::vector<SpeciesDecision> decisions;
    decisions.reserve(stats.size());
    size_t removedCnt = 0;

    for (const auto &kv : stats) {
        const EvidenceStats &st = kv.second;
        const bool passEvidence = st.strongReads >= minStrongReads;
        const bool passUnique = st.uniqueTop >= minUniqueReads;
        const bool keep = passEvidence || passUnique;
        if (keep) {
            keptSpecies.insert(kv.first);
        } else {
            ++removedCnt;
        }
        decisions.push_back({kv.first, st.occurrences, st.strongReads, st.uniqueTop,
                             st.maxScore, passEvidence, passUnique});
    }

    std::cout << "Distinct species       : " << stats.size() << std::endl;
    std::cout << "Species removed        : " << removedCnt << std::endl;
    std::cout << "Species kept           : " << keptSpecies.size() << std::endl;

    std::sort(decisions.begin(), decisions.end(),
              [](const SpeciesDecision &a, const SpeciesDecision &b) {
                  return a.occurrences > b.occurrences;
              });
    std::cout << "species\treads\tstrongReads\tuniqueTop\tmaxScore\tdecision" << std::endl;
    for (size_t r = 0; r < decisions.size() && r < MAX_REPORT_ROWS; ++r) {
        const SpeciesDecision &d = decisions[r];
        std::cout << d.spId << '\t' << d.occurrences << '\t' << d.strongReads << '\t'
                  << d.uniqueTop << '\t' << d.maxScore << '\t';
        if (d.passEvidence || d.passUnique) {
            std::cout << "kept(" << (d.passEvidence ? "evidence" : "")
                      << (d.passEvidence && d.passUnique ? "+" : "")
                      << (d.passUnique ? "unique" : "") << ")";
        } else {
            std::cout << "removed";
        }
        std::cout << std::endl;
    }
    if (decisions.size() > MAX_REPORT_ROWS) {
        std::cout << "... (" << (decisions.size() - MAX_REPORT_ROWS) << " more species)" << std::endl;
    }

    return keptSpecies;
}

} // namespace

// filter-candidates
//
// Reads a species-candidate DB (produced by "create-candidates") and writes a NEW
// candidate DB with low-quality species removed. --filter-method selects the strategy:
//   0 (default): average score (--min-avg-score) + genome coverage adjusted evenness
//                (--min-adj-evenness, needs a candidate DB with k-mer positions).
//   1          : best-evidence (--min-strong-score/--min-strong-reads) + uniqueness
//                (--min-unique-reads, --tie-ratio), kept if either passes.
// The input DB is left untouched.
int filterCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.minAvgScore = 0.0f;
    par.minAdjEvenness = 0.5f;
    par.covUseAllHits = 1;
    par.minCount = 0;
    par.filterMethod = 0;
    par.minStrongScore = 0.7f;
    par.minStrongReads = 3;
    par.minUniqueReads = 2;
    par.tieRatio = 0.99;
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string inputDb = par.filenames[0];
    const std::string dbDir = par.filenames[1];
    const std::string outputDb = par.filenames[2];

    // Validate input candidate DB (data + index)
    if (!FileUtil::fileExists(inputDb.c_str())) {
        std::cerr << "Error: candidate DB " << inputDb << " is not found." << std::endl;
        return 1;
    }
    const std::string inputIndex = inputDb + ".index";
    if (!FileUtil::fileExists(inputIndex.c_str())) {
        std::cerr << "Error: candidate DB index " << inputIndex << " is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        std::cerr << "Error: database directory " << dbDir << " is not found." << std::endl;
        return 1;
    }
    if (inputDb == outputDb) {
        std::cerr << "Error: output DB must differ from the input DB." << std::endl;
        return 1;
    }

    // Make sure the output directory exists
    const std::string outParent = FileUtil::dirName(outputDb);
    if (!outParent.empty() && !FileUtil::directoryExists(outParent.c_str())) {
        FileUtil::makeDir(outParent.c_str());
    }

    const int threadCount = par.threads <= 0 ? 1 : par.threads;
#ifdef OPENMP
    omp_set_num_threads(threadCount);
#endif

    CandidateDBReader reader(inputDb, threadCount);
    if (!reader.open(DBReader<unsigned int>::SORT_BY_ID)) {
        std::cerr << "Error: failed to open candidate DB " << inputDb << std::endl;
        return 1;
    }
    const size_t entryCount = reader.size();
    std::cout << "Candidate DB entries   : " << entryCount << std::endl;

    // --- Select filter method and decide which species to keep ---
    std::unordered_set<TaxID> keptSpecies;
    if (par.filterMethod == 1) {
        keptSpecies = filterByEvidenceAndUniqueness(reader, entryCount, par);
    } else {
        const std::unordered_map<TaxID, uint64_t> sp2genomeSize = loadGenomeSizes(dbDir);
        keptSpecies = filterByScoreAndCoverage(reader, entryCount, sp2genomeSize, par);
    }

    // --- Pass 2: write the filtered candidate DB (method-agnostic) ---
    CandidateDBWriter writer(outputDb, static_cast<unsigned int>(threadCount), 0);
    writer.open();

    std::atomic<size_t> keptCandidateCnt{0};
    std::atomic<size_t> removedCandidateCnt{0};
    std::atomic<size_t> emptiedQueryCnt{0};

    if (entryCount > 0) {
#ifdef OPENMP
#pragma omp parallel default(none) shared(reader, writer, entryCount, keptSpecies, \
        keptCandidateCnt, removedCandidateCnt, emptiedQueryCnt)
#endif
        {
            int threadIdx = 0;
#ifdef OPENMP
            threadIdx = omp_get_thread_num();
#endif
            CandidateDBEntry entry;
            size_t keptLocal = 0;
            size_t removedLocal = 0;
            size_t emptiedLocal = 0;

#ifdef OPENMP
#pragma omp for schedule(dynamic, 64)
#endif
            for (size_t i = 0; i < entryCount; ++i) {
                if (!reader.getByIndex(i, entry, threadIdx)) {
                    continue;
                }

                // Rebuild a Query so we can reuse CandidateDBWriter's serializer.
                Query query;
                query.name = entry.queryName;
                query.queryLength = static_cast<int>(entry.queryLength);
                query.queryLength2 = 0;
                query.speciesCandidates.reserve(entry.candidates.size());
                for (SpeciesCandidate &candidate : entry.candidates) {
                    if (keptSpecies.count(candidate.speciesId) != 0) {
                        query.speciesCandidates.push_back(std::move(candidate));
                        ++keptLocal;
                    } else {
                        ++removedLocal;
                    }
                }
                if (query.speciesCandidates.empty()) {
                    ++emptiedLocal;
                }

                writer.writeQuery(entry.queryId, query, threadIdx);
            }

            keptCandidateCnt += keptLocal;
            removedCandidateCnt += removedLocal;
            emptiedQueryCnt += emptiedLocal;
        }
    }

    writer.close();
    reader.close();

    std::cout << "Candidates kept        : " << keptCandidateCnt.load() << std::endl;
    std::cout << "Candidates removed     : " << removedCandidateCnt.load() << std::endl;
    std::cout << "Queries left empty     : " << emptiedQueryCnt.load() << std::endl;
    std::cout << "Filtered candidate DB written to: " << outputDb << std::endl;
    return 0;
}
