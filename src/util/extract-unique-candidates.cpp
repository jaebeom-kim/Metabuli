#include "LocalParameters.h"
#include "Parameters.h"
#include "FileUtil.h"
#include "common.h"
#include "TaxonomyWrapper.h"
#include "CandidateDBReader.h"
#include "CandidateDBWriter.h"
#include "DBReader.h"

#include <iostream>
#include <string>
#include <utility>

// extract-unique-candidates
//
// Read a species-candidate DB (from "create-candidates" / "filter-candidates")
// and write a NEW candidate DB holding only the entries (reads) that are
// UNIQUELY mapped to a user-specified species: that species is the sole top hit
// within the --tie-ratio margin. This mirrors the "uniqueTop" definition used by
// filter-candidates (a read counts as unique for a species when exactly one
// surviving species sits within best*tieRatio of the top). Useful for pulling
// out exactly the reads that give a species its unique support.
//
// The output is a candidate DB; run "view-candidates" on it to inspect the reads
// as TSV, or "classify-candidates" to classify just this subset.
int extractUniqueCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.targetTaxId = 0;
    par.tieRatio = 0.99;   // matches filter-candidates / classify-candidates
    par.minScore = 0.0f;
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string inputDb = par.filenames[0];
    const std::string dbDir = par.filenames[1];
    const std::string outputDb = par.filenames[2];

    if (par.targetTaxId == 0) {
        std::cerr << "Error: please specify the target species with --tax-id." << std::endl;
        return 1;
    }
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
    if (inputDb == outputDb) {
        std::cerr << "Error: output DB must differ from the input DB." << std::endl;
        return 1;
    }

    // Candidate DBs store INTERNAL taxonomy IDs; map the user's external species ID.
    TaxonomyWrapper *taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    const TaxID externalTaxId = par.targetTaxId;
    const TaxID targetInternal = taxonomy->getInternalTaxID(externalTaxId);
    if (targetInternal == -1) {
        std::cerr << "Error: taxon ID " << externalTaxId << " not found in the taxonomy." << std::endl;
        delete taxonomy;
        return 1;
    }

    const std::string outParent = FileUtil::dirName(outputDb);
    if (!outParent.empty() && !FileUtil::directoryExists(outParent.c_str())) {
        FileUtil::makeDir(outParent.c_str());
    }

    const float tieRatio = par.tieRatio;
    const float minScore = par.minScore;
    const int threadCount = par.threads <= 0 ? 1 : par.threads;

    CandidateDBReader reader(inputDb, threadCount);
    if (!reader.open(DBReader<unsigned int>::SORT_BY_ID)) {
        std::cerr << "Error: failed to open candidate DB " << inputDb << std::endl;
        delete taxonomy;
        return 1;
    }
    const size_t entryCount = reader.size();

    CandidateDBWriter writer(outputDb, static_cast<unsigned int>(threadCount), 0);
    writer.open();

    size_t entriesWithTarget = 0;   // entries where the target is a candidate at all
    size_t extractedCnt = 0;        // entries uniquely mapped to the target
    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (!reader.getByIndex(i, entry, 0)) {
            continue;
        }

        // Is the target a candidate on this read at all (informational)?
        bool targetPresent = false;
        for (const SpeciesCandidate &cand : entry.candidates) {
            if (cand.speciesId == targetInternal) {
                targetPresent = true;
                break;
            }
        }
        if (targetPresent) {
            ++entriesWithTarget;
        }

        // Best-first walk (candidates are stored sorted by score descending):
        // the top hit is the first candidate at/above --min-score; count how many
        // distinct species sit within best*tieRatio of it. withinMargin == 1 means
        // the top hit has no tie, i.e. the read is uniquely mapped to that species.
        float bestScore = -1.0f;
        TaxID topSp = 0;
        size_t withinMargin = 0;
        for (const SpeciesCandidate &cand : entry.candidates) {
            if (cand.idScore < minScore) {
                continue;
            }
            if (bestScore < 0.0f) {
                bestScore = cand.idScore;
                topSp = cand.speciesId;
            } else if (cand.idScore <= bestScore * tieRatio) {
                break; // outside the tie margin; all later (sorted) candidates are too
            }
            ++withinMargin;
        }

        const bool uniqueToTarget = (withinMargin == 1) && (topSp == targetInternal);
        if (!uniqueToTarget) {
            continue;
        }

        // Re-serialize the matching entry unchanged into the output candidate DB.
        Query query;
        query.name = entry.queryName;
        query.queryLength = static_cast<int>(entry.queryLength);
        query.queryLength2 = 0;
        query.speciesCandidates = std::move(entry.candidates);
        writer.writeQuery(entry.queryId, query, 0);
        ++extractedCnt;
    }

    writer.close();
    reader.close();

    std::cout << "Target species (ext/int) : " << externalTaxId << " / " << targetInternal << std::endl;
    std::cout << "Tie ratio                : " << tieRatio << std::endl;
    std::cout << "Min. score               : " << minScore << std::endl;
    std::cout << "Candidate entries        : " << entryCount << std::endl;
    std::cout << "Entries with target      : " << entriesWithTarget << std::endl;
    std::cout << "Uniquely mapped extracted: " << extractedCnt << std::endl;
    std::cout << "Written to               : " << outputDb << std::endl;

    delete taxonomy;
    return 0;
}
