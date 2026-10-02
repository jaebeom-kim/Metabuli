#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "CandidateDBReader.h"
#include "CandidateDBWriter.h"
#include "CandidatePresence.h"
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

constexpr size_t MAX_REPORT_ROWS = 50;


// Number of distinct species in a read's tie set: the top candidate plus any
// candidate within its tie margin. Mirrors chooseBestTaxonFromCandidates
// (score-dependent myTieRatio), so we can predict whether classify-candidates
// would resolve the read to a single species (size 1) or an LCA of >=2 species.
// Approximation: uses idScore only (not subScore/priority taxa), so callers
// should use the same --tie-ratio / --min-score they pass to classify-candidates.
size_t tieSetDistinctSpecies(const std::vector<SpeciesCandidate> &cands,
                             float tieRatio, float minScore) {
    float best = -1.0f;
    for (const SpeciesCandidate &c : cands) {          // candidates are best-first
        if (c.idScore >= minScore) { best = c.idScore; break; }
    }
    if (best < 0.0f) {
        return 0;
    }
    const float threshold = best * tieMarginFactor(best, tieRatio);
    std::unordered_set<TaxID> tied;
    for (const SpeciesCandidate &c : cands) {
        if (c.idScore < minScore) continue;
        if (c.idScore >= threshold) tied.insert(c.speciesId);
        else break;                                    // best-first: the rest are lower
    }
    return tied.size();
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
        // Same scaled tie margin as classify-candidates: the top is unique when the
        // runner-up is not strictly above best * tieMarginFactor.
        const bool topIsUnique = (cands.size() == 1)
            || (secondScore <= topScore * tieMarginFactor(topScore, tieRatio));

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
    par.minUniqueRatio = 0.0f;
    par.minUniqueCount = 10;
    par.dropEmptied = 0;
    par.filterMethod = 0;
    par.minStrongScore = 0.7f;
    par.minStrongReads = 3;
    par.minUniqueReads = 2;
    par.tieRatio = 0.99;   // matches setClassifyDefaults (classify-candidates)
    par.minScore = 0.0f;   // matches setClassifyDefaults (classify-candidates)
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

    // How to handle a read whose candidate species are all removed:
    //   dropEmptied == 1 : leave it empty -> unclassified (blunt, loses higher-rank recall)
    //   dropEmptied == 0 : keep it at its higher-rank LCA if it had >=2 tied species
    //                      (rescue), but discard it if it was uniquely mapped to a single
    //                      removed species (no honest fallback -> unclassified).
    const bool dropEmptied = par.dropEmptied != 0;
    const float tieRatio = par.tieRatio;
    const float minScore = par.minScore;

    std::atomic<size_t> keptCandidateCnt{0};
    std::atomic<size_t> removedCandidateCnt{0};
    std::atomic<size_t> resolvedQueryCnt{0};  // >=1 kept candidate -> resolves among kept
    std::atomic<size_t> rescuedQueryCnt{0};   // all removed but kept at higher-rank LCA
    std::atomic<size_t> discardedQueryCnt{0}; // all removed, no fallback -> unclassified

    if (entryCount > 0) {
#ifdef OPENMP
#pragma omp parallel default(none) shared(reader, writer, entryCount, keptSpecies, \
        keptCandidateCnt, removedCandidateCnt, resolvedQueryCnt, rescuedQueryCnt, \
        discardedQueryCnt, dropEmptied, tieRatio, minScore)
#endif
        {
            int threadIdx = 0;
#ifdef OPENMP
            threadIdx = omp_get_thread_num();
#endif
            CandidateDBEntry entry;
            size_t keptLocal = 0, removedLocal = 0;
            size_t resolvedLocal = 0, rescuedLocal = 0, discardedLocal = 0;

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

                bool anyKept = false;
                for (const SpeciesCandidate &c : entry.candidates) {
                    if (keptSpecies.count(c.speciesId) != 0) { anyKept = true; break; }
                }

                if (anyKept) {
                    // Case 1: keep the surviving candidates; the read resolves among them.
                    query.speciesCandidates.reserve(entry.candidates.size());
                    for (SpeciesCandidate &c : entry.candidates) {
                        if (keptSpecies.count(c.speciesId) != 0) {
                            query.speciesCandidates.push_back(std::move(c));
                            ++keptLocal;
                        } else {
                            ++removedLocal;
                        }
                    }
                    ++resolvedLocal;
                } else {
                    // Every candidate was removed.
                    removedLocal += entry.candidates.size();
                    bool rescue = false;
                    if (!dropEmptied && !entry.candidates.empty()) {
                        // Case 2 only when >=2 tied species would form a real higher-rank
                        // LCA; a unique hit to one removed species (Case 3) is discarded.
                        rescue = tieSetDistinctSpecies(entry.candidates, tieRatio, minScore) >= 2;
                    }
                    if (rescue) {
                        query.speciesCandidates = std::move(entry.candidates);
                        keptLocal += query.speciesCandidates.size();
                        removedLocal -= query.speciesCandidates.size(); // written after all, not dropped
                        ++rescuedLocal;
                    } else {
                        ++discardedLocal; // Case 3 (or --drop-emptied): unclassified
                    }
                }

                writer.writeQuery(entry.queryId, query, threadIdx);
            }

            keptCandidateCnt += keptLocal;
            removedCandidateCnt += removedLocal;
            resolvedQueryCnt += resolvedLocal;
            rescuedQueryCnt += rescuedLocal;
            discardedQueryCnt += discardedLocal;
        }
    }

    writer.close();
    reader.close();

    std::cout << "Candidates kept        : " << keptCandidateCnt.load() << std::endl;
    std::cout << "Candidates removed     : " << removedCandidateCnt.load() << std::endl;
    std::cout << "Reads resolved to kept : " << resolvedQueryCnt.load() << std::endl;
    std::cout << "Reads kept at LCA      : " << rescuedQueryCnt.load()
              << (dropEmptied ? " (disabled by --drop-emptied)" : " (>=2 tied removed species)") << std::endl;
    std::cout << "Reads discarded        : " << discardedQueryCnt.load()
              << " (uniquely mapped to a removed species)" << std::endl;
    std::cout << "Filtered candidate DB written to: " << outputDb << std::endl;
    return 0;
}
