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

} // namespace

// filter-candidates
//
// Reads a species-candidate DB (produced by "create-candidates") and writes a NEW
// candidate DB with low-quality species removed. Two independent filters:
//   1. --min-avg-score : drop species whose mean idScore across candidate reads is
//      below the threshold.
//   2. --min-adj-evenness : when the candidate DB stores k-mer positions, estimate
//      per-species genome coverage and drop species whose adjusted evenness is below
//      the threshold (default 0.5). --cov-use-all-hits controls whether coverage is
//      aggregated from every candidate hit (default) or only the top hit per read.
// The input DB is left untouched.
int filterCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.minAvgScore = 0.0f;
    par.minAdjEvenness = 0.5f;
    par.covUseAllHits = 1;
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string inputDb = par.filenames[0];
    const std::string dbDir = par.filenames[1];
    const std::string outputDb = par.filenames[2];
    const float minAvgScore = par.minAvgScore;
    const float minAdjEvenness = par.minAdjEvenness;
    const bool useAllHits = par.covUseAllHits != 0;

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

    const std::unordered_map<TaxID, uint64_t> sp2genomeSize = loadGenomeSizes(dbDir);

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
    std::cout << "Min. average score     : " << minAvgScore
              << (minAvgScore <= 0.0f ? " (score filter disabled)" : "") << std::endl;
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

    {
        CandidateDBEntry entry;
        for (size_t i = 0; i < entryCount; ++i) {
            if (!reader.getByIndex(i, entry, 0)) {
                continue;
            }
            for (size_t idx = 0; idx < entry.candidates.size(); ++idx) {
                const SpeciesCandidate &candidate = entry.candidates[idx];
                // Average score uses every candidate occurrence.
                speciesScoreSum[candidate.speciesId] += candidate.idScore;
                speciesScoreCount[candidate.speciesId] += 1;

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
    }

    // --- Decide which species survive ---
    struct SpeciesDecision {
        TaxID spId;
        double avgScore;
        uint64_t readCnt;
        double coverage;
        double adjEvenness;
        bool coverageEvaluated;
        bool removeScore;
        bool removeEvenness;
    };

    std::unordered_set<TaxID> keptSpecies;
    keptSpecies.reserve(speciesScoreCount.size());
    std::vector<SpeciesDecision> decisions;
    decisions.reserve(speciesScoreCount.size());
    size_t removedByScore = 0;
    size_t removedByEvenness = 0;
    size_t coverageEvaluated = 0;
    size_t coverageMissingGenomeSize = 0;

    for (const auto &kv : speciesScoreCount) {
        const TaxID spId = kv.first;
        const double avgScore = speciesScoreSum[spId] / static_cast<double>(kv.second);

        bool removeScore = (minAvgScore > 0.0f) && (avgScore < static_cast<double>(minAvgScore));

        bool removeEvenness = false;
        bool covEval = false;
        double coverage = 0.0;
        double adjEvenness = 0.0;
        uint64_t readCnt = 0;
        const auto binIt = sp2bins.find(spId);
        if (binIt != sp2bins.end()) {
            const auto gsIt = sp2genomeSize.find(spId);
            if (gsIt != sp2genomeSize.end() && gsIt->second > 0) {
                readCnt = sp2readCnt[spId];
                const CovMetric metric = computeCoverageMetric(
                    binIt->second, readCnt, sp2readLen[spId], gsIt->second);
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

        if (removeScore) ++removedByScore;
        else if (removeEvenness) ++removedByEvenness;

        if (!removeScore && !removeEvenness) {
            keptSpecies.insert(spId);
        }
        decisions.push_back({spId, avgScore, readCnt, coverage, adjEvenness,
                             covEval, removeScore, removeEvenness});
    }

    std::cout << "Distinct species       : " << speciesScoreCount.size() << std::endl;
    std::cout << "Coverage evaluated for : " << coverageEvaluated << " species" << std::endl;
    if (coverageMissingGenomeSize > 0) {
        std::cout << "Coverage skipped (no genome size) : " << coverageMissingGenomeSize << " species" << std::endl;
    }
    std::cout << "Species removed (score)    : " << removedByScore << std::endl;
    std::cout << "Species removed (evenness) : " << removedByEvenness << std::endl;
    std::cout << "Species kept               : " << keptSpecies.size() << std::endl;

    // Per-species report (capped to keep the log readable).
    constexpr size_t MAX_REPORT_ROWS = 50;
    std::sort(decisions.begin(), decisions.end(),
              [](const SpeciesDecision &a, const SpeciesDecision &b) {
                  return a.readCnt > b.readCnt;
              });
    std::cout << "species\treads\tavgScore\tcoverage\tadjEvenness\tdecision" << std::endl;
    for (size_t r = 0; r < decisions.size() && r < MAX_REPORT_ROWS; ++r) {
        const SpeciesDecision &d = decisions[r];
        std::cout << d.spId << '\t' << d.readCnt << '\t' << d.avgScore << '\t';
        if (d.coverageEvaluated) {
            std::cout << d.coverage << '\t' << d.adjEvenness << '\t';
        } else {
            std::cout << "-\t-\t";
        }
        if (d.removeScore) std::cout << "removed(score)";
        else if (d.removeEvenness) std::cout << "removed(evenness)";
        else std::cout << "kept";
        std::cout << std::endl;
    }
    if (decisions.size() > MAX_REPORT_ROWS) {
        std::cout << "... (" << (decisions.size() - MAX_REPORT_ROWS) << " more species)" << std::endl;
    }

    // --- Pass 2: write the filtered candidate DB ---
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
