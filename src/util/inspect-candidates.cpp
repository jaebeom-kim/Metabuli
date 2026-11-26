#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "TaxonomyWrapper.h"
#include "CandidateDBReader.h"
#include "DBReader.h"

#include <algorithm>
#include <bitset>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <dirent.h>
#include <fstream>
#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>

// inspect-candidates
//
// Compute per-species metrics from a species-candidate DB (from
// "create-candidates" / "filter-candidates") and dump them to a TSV, one row
// per species. A read-level inspector is "view-candidates"; this is the
// species-level counterpart for exploring filter thresholds before running
// "filter-candidates". Metrics use the same per-read accounting as
// filter-candidates' score gate: each read credits its top candidate plus
// every candidate within the tie margin (tieMarginFactor with --tie-ratio).
//
// Metric groups are selected by flags; with no flag, all are computed:
//   --avg-score    : credited reads and their mean idScore
//   --unique-ratio : unique-top reads (no tie) and their fraction of credited reads
//   --macro-cov    : genome-window occupancy of one locus per credited read
//                    (coverageLocus of the candidate's k-mer position bins):
//                    n (reads with a locus), k (distinct windows, locus >> 8),
//                    W (genomeWindowCount), macroCov = k/W, ratio = k/E with
//                    E = W(1-(1-1/W)^n), and log10 of the uniform-placement
//                    tail bound P(<= k windows) <= C(W,k) * (k/W)^n.
//                    Needs a --kmer-pos DB and genome sizes.

namespace {

// Load speciesId -> genome size from species2genomeSize.tsv (col0 = speciesId,
// col2 = genome size). The file sits next to each k-mer list: at the top level
// for classic DBs, under kmer-lists/<name>/ for sub-DB layout DBs (matching
// Classifier::parseSp2GenomeSize, which reads it per sub-DB).
std::unordered_map<TaxID, uint64_t> loadGenomeSizesForInspect(const std::string &dbDir) {
    std::unordered_map<TaxID, uint64_t> sp2genomeSize;

    std::vector<std::string> dirs = {dbDir};
    const std::string listsDir = dbDir + "/" + SUBDB_DIR_NAME;
    DIR *dp = opendir(listsDir.c_str());
    if (dp != nullptr) {
        struct dirent *entry;
        while ((entry = readdir(dp)) != nullptr) {
            const std::string name = entry->d_name;
            if (name == "." || name == "..") {
                continue;
            }
            dirs.push_back(listsDir + "/" + name);
        }
        closedir(dp);
    }

    bool anyOpened = false;
    for (const std::string &dir : dirs) {
        const std::string fileName = dir + "/species2genomeSize.tsv";
        std::ifstream infile(fileName);
        if (!infile.is_open()) {
            continue;
        }
        anyOpened = true;
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
    }
    if (!anyOpened) {
        std::cout << "Warning: could not open species2genomeSize.tsv in " << dbDir
                  << " (or its " << SUBDB_DIR_NAME << " subdirectories);"
                  << " --macro-cov columns will be empty." << std::endl;
    }
    return sp2genomeSize;
}

// log10 of the uniform-placement tail bound P(n loci occupy <= k of W windows)
// <= C(W,k) * (k/W)^n, clamped at 0 (the bound exceeds 1 when the occupancy is
// unremarkable). Exact enough at small n (where the normal approximation of
// the occupancy count is unusable) and conservative at large n.
double breadthLog10P(int W, int k, uint32_t n) {
    const double logC = (std::lgamma(W + 1.0) - std::lgamma(k + 1.0) - std::lgamma(W - k + 1.0))
                        / std::log(10.0);
    const double log10P = logC + static_cast<double>(n) * std::log10(static_cast<double>(k) / static_cast<double>(W));
    return log10P > 0.0 ? 0.0 : log10P;
}

struct SpeciesMetrics {
    uint64_t reads = 0;        // credited reads (top hit + within-margin ties)
    double scoreSum = 0.0;     // idScore sum over credited reads
    uint64_t uniqueTop = 0;    // reads where this species is the sole top (no tie)
    uint32_t lociReads = 0;    // credited reads with a locus (n)
    std::bitset<256> windows;  // distinct windows (locus >> 8) they occupy (k)
};

} // namespace

int inspectCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.tieRatio = 0.99;   // matches filter-candidates / classify-candidates
    par.inspectMacroCov = 0;
    par.inspectAvgScore = 0;
    par.inspectUniqueRatio = 0;
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string inputDb = par.filenames[0];
    const std::string dbDir = par.filenames[1];
    const std::string outputTsv = par.filenames[2];

    if (!FileUtil::fileExists(inputDb.c_str())) {
        std::cerr << "Error: candidate DB " << inputDb << " is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::fileExists((inputDb + ".index").c_str())) {
        std::cerr << "Error: candidate DB index " << inputDb << ".index is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::directoryExists(dbDir.c_str())) {
        std::cerr << "Error: DB directory " << dbDir << " is not found." << std::endl;
        return 1;
    }

    // No metric flag means "all metrics".
    const bool noFlag = par.inspectMacroCov == 0 && par.inspectAvgScore == 0 && par.inspectUniqueRatio == 0;
    const bool doMacroCov = noFlag || par.inspectMacroCov != 0;
    const bool doAvgScore = noFlag || par.inspectAvgScore != 0;
    const bool doUniqueRatio = noFlag || par.inspectUniqueRatio != 0;
    const float tieRatio = par.tieRatio;

    // Candidate DBs store INTERNAL taxonomy IDs; load the source DB's taxonomy
    // to convert them back to original (external) IDs and names.
    TaxonomyWrapper *taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    std::unordered_map<TaxID, uint64_t> sp2genomeSize;
    if (doMacroCov) {
        sp2genomeSize = loadGenomeSizesForInspect(dbDir);
    }

    const std::string outParent = FileUtil::dirName(outputTsv);
    if (!outParent.empty() && !FileUtil::directoryExists(outParent.c_str())) {
        FileUtil::makeDir(outParent.c_str());
    }
    std::ofstream out(outputTsv);
    if (!out.is_open()) {
        std::cerr << "Error: could not open output file " << outputTsv << std::endl;
        delete taxonomy;
        return 1;
    }

    const int threadCount = par.threads <= 0 ? 1 : par.threads;
    CandidateDBReader reader(inputDb, threadCount);
    if (!reader.open(DBReader<unsigned int>::SORT_BY_ID)) {
        std::cerr << "Error: failed to open candidate DB " << inputDb << std::endl;
        delete taxonomy;
        return 1;
    }
    const size_t entryCount = reader.size();

    // Per-read accounting, identical to filter-candidates' score gate: the
    // first candidate is the top hit (lists are stored best-first), and every
    // candidate above bestScore * tieMarginFactor is credited alongside it.
    std::unordered_map<TaxID, SpeciesMetrics> stats;
    size_t readsWithCandidates = 0;
    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (!reader.getByIndex(i, entry, 0)) {
            continue;
        }
        if (entry.candidates.empty()) {
            continue;
        }
        ++readsWithCandidates;

        const float bestScore = entry.candidates[0].idScore;
        const float tieThreshold = bestScore * tieMarginFactor(bestScore, tieRatio);
        size_t withinMargin = 0;
        for (const SpeciesCandidate &cand : entry.candidates) {
            if (withinMargin > 0 && cand.idScore <= tieThreshold) {
                break; // best-first: the rest are below the margin too
            }
            SpeciesMetrics &st = stats[cand.speciesId];
            st.reads += 1;
            st.scoreSum += cand.idScore;
            if (doMacroCov) {
                const uint16_t locus = coverageLocus(cand.kmerPositions);
                if (locus != 0) {
                    st.lociReads += 1;
                    st.windows.set(locus >> 8);
                }
            }
            ++withinMargin;
        }
        if (withinMargin == 1) {
            stats[entry.candidates[0].speciesId].uniqueTop += 1;
        }
    }
    reader.close();

    // Header
    out << "species\tname\treads";
    if (doAvgScore) out << "\tavgScore";
    if (doUniqueRatio) out << "\tuniqueTop\tuniqueRatio";
    if (doMacroCov) out << "\tlociReads\twindows\tgenomeWindows\texpectedWindows\tmacroCov\tratio\tlog10P";
    out << "\n";

    // Rows, most-credited species first (taxId tiebreak for determinism).
    std::vector<std::pair<TaxID, const SpeciesMetrics *>> rows;
    rows.reserve(stats.size());
    for (const auto &kv : stats) {
        rows.emplace_back(kv.first, &kv.second);
    }
    std::sort(rows.begin(), rows.end(),
              [](const std::pair<TaxID, const SpeciesMetrics *> &a,
                 const std::pair<TaxID, const SpeciesMetrics *> &b) {
                  if (a.second->reads != b.second->reads) {
                      return a.second->reads > b.second->reads;
                  }
                  return a.first < b.first;
              });

    for (const auto &row : rows) {
        const TaxID internalId = row.first;
        const SpeciesMetrics &st = *row.second;
        const TaxonNode *node = taxonomy->taxonNode(internalId);
        out << taxonomy->getOriginalTaxID(internalId) << '\t'
            << (node != nullptr ? taxonomy->getString(node->nameIdx) : "-") << '\t'
            << st.reads;
        if (doAvgScore) {
            out << '\t' << (st.reads > 0 ? st.scoreSum / static_cast<double>(st.reads) : 0.0);
        }
        if (doUniqueRatio) {
            out << '\t' << st.uniqueTop << '\t'
                << (st.reads > 0 ? static_cast<double>(st.uniqueTop) / static_cast<double>(st.reads) : 0.0);
        }
        if (doMacroCov) {
            const auto gsIt = sp2genomeSize.find(internalId);
            const uint32_t n = st.lociReads;
            const int k = static_cast<int>(st.windows.count());
            if (gsIt == sp2genomeSize.end() || gsIt->second == 0 || n == 0) {
                // No genome size or no loci: occupancy is not evaluable.
                out << '\t' << n << '\t' << k << "\t-\t-\t-\t-\t-";
            } else {
                const int windowCnt = genomeWindowCount(gsIt->second);
                const double W = static_cast<double>(windowCnt);
                const double E = W * (1.0 - std::pow(1.0 - 1.0 / W, static_cast<double>(n)));
                out << '\t' << n << '\t' << k << '\t' << windowCnt
                    << '\t' << E
                    << '\t' << static_cast<double>(k) / W
                    << '\t' << (E > 0.0 ? static_cast<double>(k) / E : 0.0)
                    << '\t' << (windowCnt >= 2 ? breadthLog10P(windowCnt, k, n) : 0.0);
            }
        }
        out << '\n';
    }
    out.close();
    delete taxonomy;

    std::cout << "Candidate entries      : " << entryCount << std::endl;
    std::cout << "Reads with candidates  : " << readsWithCandidates << std::endl;
    std::cout << "Distinct species       : " << rows.size() << std::endl;
    std::cout << "Tie ratio              : " << tieRatio << std::endl;
    std::cout << "Written to             : " << outputTsv << std::endl;
    return 0;
}
