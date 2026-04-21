#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "Reporter.h"
#include "Classifier.h"

std::unordered_set<TaxID> getSpeciesAndLowerTaxIds(
    const std::vector<TaxID>& speciesTaxIds,
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren);

void setRefineReportDefault(LocalParameters & par) {
    par.minAvgScore = 0;
    par.minCladeCount = 0;
    par.minCladeProportion = 0;

}

int refine_report(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    std::string resultFileName = par.filenames[0];
    std::string dbDir = par.filenames[1];
    std::string outFileName = par.filenames[2];

    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, dbDir + "/taxonomy/");

    Reporter reporter(par, taxonomy, outFileName);

    // Open the read-by-read result file and read the classification results
    std::ifstream resultFile(resultFileName);
    if (!resultFile.is_open()) {
        std::cerr << "Error: could not open result file " << resultFileName << std::endl;
        return 1;
    }

    std::unordered_map<TaxID, TaxID> ex2inTaxId;
    taxonomy->getExternal2internalTaxID(ex2inTaxId);

    // Get the set of species or below tax IDs
    std::vector<TaxID> speciesTaxIds;
    int maxTaxId = taxonomy->getMaxTaxID();
    for (int i = 0; i <= maxTaxId; i++) {
        const char* rank = taxonomy->getString(taxonomy->taxonNodes[i].rankIdx);
        if (strcmp(rank, "species") == 0) {
            speciesTaxIds.push_back(taxonomy->taxonNodes[i].taxId);
        }
    }
    const std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    std::unordered_set<TaxID> speciesAndLowerTaxIds = getSpeciesAndLowerTaxIds(speciesTaxIds, parentToChildren);


    std::unordered_map<TaxID, unsigned int> taxonCounts;
    std::unordered_map<TaxID, double> taxon2avgScore;
    size_t totalReads = 0;
    size_t totalClassifiedReads = 0;
    std::string line;
    int rankCol = -1;
    int taxidCol = -1;
    int scoreCol = -1;
    std::vector<std::string> items;
    while (std::getline(resultFile, line)) {
        if (line[0] == '#') {
            // Header line, find the column index for rank
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
            for (size_t i = 0; i < columns.size(); i++) {
                if (columns[i] == "rank") { rankCol = i; }
                if (columns[i] == "taxID") { taxidCol = i; }
                if (columns[i] == "score") { scoreCol = i; }
            }
            continue;
        }

        if (line.empty()) {
            continue;
        }

        totalReads++;

        items = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
        if (rankCol == -1 || static_cast<size_t>(rankCol) >= items.size()) {
            std::cerr << "Error: rank column not found in result file header." << std::endl;
            return 1;
        }

        TaxID taxID = ex2inTaxId[std::stoi(items[taxidCol])];
        if (taxID != 0) {
            totalClassifiedReads++;
        }
        taxonCounts[taxID]++;

        // Check if species or below, and if so, add the score to taxon2avgScore
        // if (par.minAvgScore > 0 && speciesAndLowerTaxIds.find(taxID) != speciesAndLowerTaxIds.end()) {
        //     double score = std::stod(items[scoreCol]);
        //     TaxID spTaxId = taxonomy->getTaxIdAtRank(taxID, "species");
        //     sp2avgScore[spTaxId] += score;
        // }

        if (par.minAvgScore > 0 && taxID != 0) {
            double score = std::stod(items[scoreCol]);

            if (speciesAndLowerTaxIds.find(taxID) != speciesAndLowerTaxIds.end()) {
                // If the taxID is species or below, add the score to the corresponding species
                TaxID spTaxId = taxonomy->getTaxIdAtRank(taxID, "species");
                taxon2avgScore[spTaxId] += score;
            } else {
                // If the taxID is above species, use it directly to add the score
                taxon2avgScore[taxID] += score;
            }
        }
    }

    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);

    // Calculate average score for each species
    if (par.minAvgScore > 0) {        
        // for (const auto& kv : sp2avgScore) {
        //     sp2avgScore[kv.first] = kv.second / cladeCounts[kv.first].cladeCount;
        // }

        Classifier::rollUpScore(taxon2avgScore, parentToChildren, 1);
        for (const auto& kv : taxon2avgScore) {
            taxon2avgScore[kv.first] = kv.second / cladeCounts[kv.first].cladeCount;
        }
        std::unordered_map<TaxID, double> taxon2avgScore_copy = taxon2avgScore;
        
        taxon2avgScore.clear();
        taxonCounts.clear();
        cladeCounts.clear();
        totalClassifiedReads = 0;

        // Filter out species that do not meet the minimum average score threshold
        resultFile.clear();
        resultFile.seekg(0, std::ios::beg);
        while (std::getline(resultFile, line)) {
            if (line[0] == '#') {
                continue;
            }

            if (line.empty()) {
                continue;
            }

            items = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
            TaxID taxID = ex2inTaxId[std::stoi(items[taxidCol])];

            if (taxID != 0) {
                double score = std::stod(items[scoreCol]);
                TaxID mappableTaxId = taxID;
                if (speciesAndLowerTaxIds.find(taxID) != speciesAndLowerTaxIds.end()) {
                    mappableTaxId = taxonomy->getTaxIdAtRank(taxID, "species");
                }

                if (taxon2avgScore_copy[mappableTaxId] >= par.minAvgScore) {
                    // If the taxID (or its corresponding species) meets the minimum average score, keep it
                    taxon2avgScore[mappableTaxId] += score;
                    taxonCounts[taxID]++;
                    totalClassifiedReads++;
                } else {
                    // Otherwise, treat it as unclassified (taxID = 0)
                    taxonCounts[0]++;
                }
            }
        }

        cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);

        Classifier::rollUpScore(taxon2avgScore, parentToChildren, 1);
        // Recalculate average score for each species after filtering by minimum average score
        for (const auto& kv : taxon2avgScore) {
            taxon2avgScore[kv.first] = kv.second / cladeCounts[kv.first].cladeCount;
        }
        
    } 

    if (par.minCladeCount > 0 || par.minCladeProportion > 0) {
        // Get species to filter out based on clade count and proportion
        // std::unordered_set<TaxID> rejectedSpeciesIds;
        std::vector<TaxID> rejectedSpeciesIds;
        for (size_t i = 0; i < speciesTaxIds.size(); i++) {
            TaxID spTaxId = speciesTaxIds[i];
            auto it = cladeCounts.find(spTaxId);
            if (it == cladeCounts.end()) {
                continue;
            }
            unsigned int cladeCount = it->second.cladeCount;
            double cladeProportion = cladeCount / double(totalClassifiedReads);
            if ((par.minCladeCount > 0 && cladeCount < static_cast<unsigned int>(par.minCladeCount)) ||
                (par.minCladeProportion > 0 && cladeProportion < par.minCladeProportion)) {
                rejectedSpeciesIds.push_back(spTaxId);
            }
        }

        std::unordered_set<TaxID> taxIdToClean = getSpeciesAndLowerTaxIds(rejectedSpeciesIds, parentToChildren);
        for (TaxID taxId : taxIdToClean) {
            taxonCounts.erase(taxId);
        }

        cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);
    }    

    reporter.writeReportFile(totalReads, cladeCounts, taxon2avgScore, ReportType::Default);

    return 0;
}


std::unordered_set<TaxID> getSpeciesAndLowerTaxIds(
    const std::vector<TaxID>& inputTaxIds,
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren) 
{
    std::unordered_set<TaxID> result;
    std::vector<TaxID> stack = inputTaxIds;

    // 1. Add the input species tax IDs to the result set
    for (TaxID taxId : stack) {
        result.insert(taxId);
    }

    // 2. Traverse downwards to grab all subspecies, strains, and isolates
    while (!stack.empty()) {
        TaxID currentTaxId = stack.back();
        stack.pop_back();

        // Check if this TaxID has any children in the map
        auto it = parentToChildren.find(currentTaxId);
        if (it != parentToChildren.end()) {
            for (TaxID childId : it->second) {
                // Insert the child into the result.
                // Checking if it already exists prevents infinite loops 
                // just in case the taxonomy data has a corrupted cycle.
                if (result.insert(childId).second) { 
                    stack.push_back(childId); // Push child to explore its descendants later
                }
            }
        }
    }

    return result;
}

