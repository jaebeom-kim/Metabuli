#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "TaxonomyWrapper.h"
#include "CandidateDBReader.h"
#include "DBReader.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

// view-candidates
//
// Development helper: dump a species-candidate DB (from "create-candidates" /
// "filter-candidates") to a human-readable TSV so intermediate state can be
// inspected. One row per (query, candidate); queries with no candidates get a
// single row with empty candidate columns. Works for both v1 (no positions) and
// v2 (with k-mer positions) candidate DBs.
int viewCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
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

    // Candidate DBs store INTERNAL taxonomy IDs; load the source DB's taxonomy so
    // speciesId and the taxCnt taxIDs can be converted back to original (external) IDs.
    TaxonomyWrapper *taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    const std::string outParent = FileUtil::dirName(outputTsv);
    if (!outParent.empty() && !FileUtil::directoryExists(outParent.c_str())) {
        FileUtil::makeDir(outParent.c_str());
    }

    std::ofstream out(outputTsv);
    if (!out.is_open()) {
        std::cerr << "Error: could not open output file " << outputTsv << std::endl;
        return 1;
    }

    const int threadCount = par.threads <= 0 ? 1 : par.threads;
    CandidateDBReader reader(inputDb, threadCount);
    if (!reader.open(DBReader<unsigned int>::SORT_BY_ID)) {
        std::cerr << "Error: failed to open candidate DB " << inputDb << std::endl;
        return 1;
    }
    const size_t entryCount = reader.size();

    out << "queryId\tqueryName\tqueryLength\tnumCandidates\tcandRank\tspeciesId"
           "\tidScore\tsubScore\tlogE\ttaxCnt\tnumPositions\tpositions\n";

    size_t candidateRows = 0;
    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (!reader.getByIndex(i, entry, 0)) {
            std::cerr << "Warning: could not read entry at index " << i << std::endl;
            continue;
        }

        const size_t numCandidates = entry.candidates.size();
        if (numCandidates == 0) {
            out << entry.queryId << '\t' << entry.queryName << '\t' << entry.queryLength
                << '\t' << 0 << "\t-\t-\t-\t-\t-\t-\t0\t\n";
            continue;
        }

        for (size_t c = 0; c < numCandidates; ++c) {
            const SpeciesCandidate &cand = entry.candidates[c];

            std::ostringstream taxCntStr;
            for (size_t t = 0; t < cand.taxCnt.size(); ++t) {
                if (t != 0) taxCntStr << ';';
                taxCntStr << taxonomy->getOriginalTaxID(cand.taxCnt[t].first) << ':' << cand.taxCnt[t].second;
            }

            std::ostringstream posStr;
            for (size_t p = 0; p < cand.kmerPositions.size(); ++p) {
                if (p != 0) posStr << ',';
                posStr << cand.kmerPositions[p];
            }

            out << entry.queryId << '\t' << entry.queryName << '\t' << entry.queryLength
                << '\t' << numCandidates << '\t' << c << '\t' << taxonomy->getOriginalTaxID(cand.speciesId)
                << '\t' << cand.idScore << '\t' << cand.subScore << '\t' << cand.logE
                << '\t' << taxCntStr.str()
                << '\t' << cand.kmerPositions.size() << '\t' << posStr.str() << '\n';
            ++candidateRows;
        }
    }

    reader.close();
    out.close();
    delete taxonomy;

    std::cout << "Entries    : " << entryCount << std::endl;
    std::cout << "Candidates : " << candidateRows << std::endl;
    std::cout << "Written to : " << outputTsv << std::endl;
    return 0;
}
