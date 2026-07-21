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

// -------- Filter method 0: iterative top-hit average score + read count + genome coverage --------
//
// Iterative because removing a species reshuffles which species is each read's top hit.
// Each round assigns every read to its best surviving candidate and scores species over
// the reads they currently win (matching report.tsv's per-species avg_score, not an
// average over all candidate occurrences). A species is dropped when it wins no reads,
// when its mean top-hit idScore is below --min-avg-score, or when it is the top hit for
// fewer than --min-count reads; reads whose winner is dropped are rescued onto their
// next-best survivor. The loop repeats to a fixpoint. Coverage/evenness
// (--min-adj-evenness, needs k-mer positions + genome sizes) is a coarser outer gate:
// after score+count converge it is evaluated once on the survivors, and if it removes
// any species the score+count loop is re-converged. Kept species therefore satisfy all
// active thresholds simultaneously under one consistent assignment.
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

    // --- Load each read's candidate (speciesId, idScore) list once, best-first ---
    // Score and count are then iterated in RAM (no k-mer positions retained);
    // coverage is recomputed by re-reading the DB only when its gate runs.
    std::vector<std::vector<std::pair<TaxID, float>>> reads(entryCount);
    std::unordered_set<TaxID> alive;
    {
        CandidateDBEntry entry;
        for (size_t i = 0; i < entryCount; ++i) {
            if (!reader.getByIndex(i, entry, 0)) {
                continue;
            }
            std::vector<std::pair<TaxID, float>> &r = reads[i];
            r.reserve(entry.candidates.size());
            for (const SpeciesCandidate &candidate : entry.candidates) {
                r.emplace_back(candidate.speciesId, candidate.idScore);
                alive.insert(candidate.speciesId);
            }
        }
    }
    const size_t distinctSpecies = alive.size();

    // Per-species score/count from the most recent assignment (kept for reporting).
    std::unordered_map<TaxID, double> speciesScoreSum;
    std::unordered_map<TaxID, uint64_t> speciesScoreCount;
    size_t removedByScore = 0, removedByCount = 0, removedByNoReads = 0, removedByEvenness = 0;
    size_t scoreRounds = 0, coverageRounds = 0;

    // Iterate to a fixpoint: reassign each read to its best surviving candidate,
    // rescore, and drop species below threshold. Coverage/evenness is a coarser
    // outer gate that re-converges the score+count loop whenever it removes a
    // species (removal reshuffles which species is each read's top hit).
    while (true) {
        // ---- Inner loop: top-hit score + count to a fixpoint ----
        while (true) {
            speciesScoreSum.clear();
            speciesScoreCount.clear();
            for (const std::vector<std::pair<TaxID, float>> &r : reads) {
                for (const std::pair<TaxID, float> &hit : r) {
                    if (alive.count(hit.first) != 0) {
                        // First surviving candidate = this read's top hit.
                        speciesScoreSum[hit.first] += hit.second;
                        speciesScoreCount[hit.first] += 1;
                        break;
                    }
                }
            }

            std::vector<TaxID> toRemove;
            for (const TaxID sp : alive) {
                const auto cntIt = speciesScoreCount.find(sp);
                const uint64_t c = (cntIt == speciesScoreCount.end()) ? 0 : cntIt->second;
                if (c == 0) {                       // no longer any read's top hit
                    toRemove.push_back(sp);
                    ++removedByNoReads;
                    continue;
                }
                const double avg = speciesScoreSum[sp] / static_cast<double>(c);
                if (minAvgScore > 0.0f && avg < static_cast<double>(minAvgScore)) {
                    toRemove.push_back(sp);
                    ++removedByScore;
                    continue;
                }
                if (minCount > 0 && c < minCount) {
                    toRemove.push_back(sp);
                    ++removedByCount;
                    continue;
                }
            }
            if (toRemove.empty()) {
                break;
            }
            for (const TaxID sp : toRemove) {
                alive.erase(sp);
            }
            ++scoreRounds;
        }

        // ---- Outer gate: coverage / adjusted evenness on the current survivors ----
        if (minAdjEvenness <= 0.0f) {
            break; // coverage filter disabled
        }

        std::unordered_map<TaxID, std::vector<uint8_t>> sp2bins;
        std::unordered_map<TaxID, uint64_t> sp2readCnt;
        std::unordered_map<TaxID, uint64_t> sp2readLen;
        {
            CandidateDBEntry entry;
            for (size_t i = 0; i < entryCount; ++i) {
                if (!reader.getByIndex(i, entry, 0)) {
                    continue;
                }
                bool topTaken = false;
                for (const SpeciesCandidate &candidate : entry.candidates) {
                    if (alive.count(candidate.speciesId) == 0) {
                        continue; // removed species do not contribute
                    }
                    const bool isTop = !topTaken;
                    topTaken = true;
                    if ((useAllHits || isTop) && !candidate.kmerPositions.empty()) {
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
                    if (!useAllHits) {
                        break; // top surviving candidate only
                    }
                }
            }
        }

        std::vector<TaxID> toRemove;
        for (const TaxID sp : alive) {
            const auto binIt = sp2bins.find(sp);
            if (binIt == sp2bins.end()) {
                continue; // no k-mer positions -> coverage not evaluable
            }
            const auto gsIt = sp2genomeSize.find(sp);
            if (gsIt == sp2genomeSize.end() || gsIt->second == 0) {
                continue; // no genome size -> skip coverage for this species
            }
            const CovMetric metric = computeCoverageMetric(
                binIt->second, sp2readCnt[sp], sp2readLen[sp], gsIt->second);
            if (metric.adjustedEvenness < static_cast<double>(minAdjEvenness)) {
                toRemove.push_back(sp);
                ++removedByEvenness;
            }
        }
        if (toRemove.empty()) {
            break;
        }
        for (const TaxID sp : toRemove) {
            alive.erase(sp);
        }
        ++coverageRounds;
        // Removing species reshuffles top hits -> re-converge score+count.
    }

    // --- Summary ---
    std::cout << "Distinct species       : " << distinctSpecies << std::endl;
    std::cout << "Score/count rounds     : " << scoreRounds << std::endl;
    if (minAdjEvenness > 0.0f) {
        std::cout << "Coverage gate rounds   : " << coverageRounds << std::endl;
    }
    std::cout << "Species removed (score)     : " << removedByScore << std::endl;
    std::cout << "Species removed (count)     : " << removedByCount << std::endl;
    std::cout << "Species removed (no reads)  : " << removedByNoReads << std::endl;
    std::cout << "Species removed (evenness)  : " << removedByEvenness << std::endl;
    std::cout << "Species kept                : " << alive.size() << std::endl;

    // Per-species report of survivors under the final assignment (capped).
    std::vector<std::pair<TaxID, uint64_t>> survivors;
    survivors.reserve(alive.size());
    for (const TaxID sp : alive) {
        const auto cntIt = speciesScoreCount.find(sp);
        survivors.emplace_back(sp, cntIt == speciesScoreCount.end() ? 0 : cntIt->second);
    }
    std::sort(survivors.begin(), survivors.end(),
              [](const std::pair<TaxID, uint64_t> &a, const std::pair<TaxID, uint64_t> &b) {
                  return a.second > b.second;
              });
    std::cout << "species\ttopHitReads\tavgScore" << std::endl;
    for (size_t r = 0; r < survivors.size() && r < MAX_REPORT_ROWS; ++r) {
        const TaxID sp = survivors[r].first;
        const uint64_t c = survivors[r].second;
        const double avg = (c > 0) ? speciesScoreSum[sp] / static_cast<double>(c) : 0.0;
        std::cout << sp << '\t' << c << '\t' << avg << std::endl;
    }
    if (survivors.size() > MAX_REPORT_ROWS) {
        std::cout << "... (" << (survivors.size() - MAX_REPORT_ROWS) << " more kept species)" << std::endl;
    }

    return alive;
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
