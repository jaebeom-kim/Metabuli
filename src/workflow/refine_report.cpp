#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "Reporter.h"
#include "Classifier.h"

std::unordered_set<TaxID> getSelfAndChildren(
    const std::vector<TaxID>& taxIds,
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren);

void setRefineReportDefault(LocalParameters & par) {
    par.minAvgScore = 0;
    par.minCladeCount = 0;
    par.minCladeProportion = 0;
    par.outFilteredResults = "";
    par.maxEValue = 0;
    par.excludeTaxid = "";
    par.minScore = 0;

}

int refine_report(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    std::string resultFileName = par.filenames[0];
    std::string dbDir = par.filenames[1];
    std::string outFileName = par.filenames[2];

    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, dbDir + "/taxonomy/");
    std::unordered_map<TaxID, TaxID> ex2inTaxId;
    taxonomy->getExternal2internalTaxID(ex2inTaxId);
    const std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();

    Reporter reporter(par, taxonomy, outFileName);

    std::vector<TaxID> denylist;
    std::unordered_set<TaxID> allDenyTaxIds;
    if (!par.excludeTaxid.empty()) {
        std::vector<std::string> contams = Util::split(par.excludeTaxid, ",");
        for (const std::string &contam : contams) {
            // Check if the tax ID exists in the taxonomy
            TaxID externalTaxId = std::stoi(contam);
            if (ex2inTaxId.find(externalTaxId) == ex2inTaxId.end()) {
                std::cerr << "Warning: tax ID " << externalTaxId << " in the exclude list does not exist in the taxonomy and will be ignored." << std::endl;
                continue;
            }
            denylist.push_back(ex2inTaxId[externalTaxId]);
        }
        allDenyTaxIds = getSelfAndChildren(denylist, parentToChildren);
    }
    
    // Get the set of species or below tax IDs
    std::vector<TaxID> speciesTaxIds;
    int maxTaxId = taxonomy->getMaxTaxID();
    for (int i = 0; i <= maxTaxId; i++) {
        const char* rank = taxonomy->getString(taxonomy->taxonNodes[i].rankIdx);
        if (strcmp(rank, "species") == 0) {
            speciesTaxIds.push_back(taxonomy->taxonNodes[i].taxId);
        }
    }
    
    std::unordered_set<TaxID> speciesAndLowerTaxIds = getSelfAndChildren(speciesTaxIds, parentToChildren);


    unordered_map<TaxID, TaxonCounts> cladeCounts;
    std::unordered_map<TaxID, unsigned int> taxonCounts;
    std::unordered_map<TaxID, double> taxon2avgScore;
    std::unordered_map<TaxID, double> taxon2avgScore_copy;
    size_t totalReads = 0;
    size_t totalClassifiedReads = 0;
    std::string line;
    int taxidCol = -1;
    int scoreCol = -1;
    int eValueCol = -1;
    int lineageCol = -1;
    std::vector<std::string> items;

    // Open the read-by-read result file and read the classification results
    std::ifstream resultFile(resultFileName);
    if (!resultFile.is_open()) {
        std::cerr << "Error: could not open result file " << resultFileName << std::endl;
        return 1;
    }

    

    // Calculate clade counts and average scores for each taxon
    if (par.minAvgScore > 0) {
        while (std::getline(resultFile, line)) {
            if (line[0] == '#') {
                std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
                for (size_t i = 0; i < columns.size(); i++) {
                    if (columns[i] == "taxID") { taxidCol = i; }
                    if (columns[i] == "score") { scoreCol = i; }
                }
                continue;
            }

            if (line.empty()) continue;
        
            if (line[0] == '0') {
                taxonCounts[0]++;
                continue;
            }
        
            items = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
            double score = std::stod(items[scoreCol]);
            TaxID mappableTaxId = ex2inTaxId[std::stoi(items[taxidCol])];
            taxonCounts[mappableTaxId]++;
        
            if (speciesAndLowerTaxIds.find(mappableTaxId) != speciesAndLowerTaxIds.end()) {
                mappableTaxId = taxonomy->getTaxIdAtRank(mappableTaxId, "species");
            }

            taxon2avgScore[mappableTaxId] += score;            
        }

        cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);
        Classifier::rollUpScore(taxon2avgScore, parentToChildren, 1);
        for (const auto& kv : taxon2avgScore) {
            auto it = cladeCounts.find(kv.first);
            if (it != cladeCounts.end() && it->second.cladeCount > 0) {
                taxon2avgScore[kv.first] = kv.second / static_cast<double>(it->second.cladeCount);
            } else {
                taxon2avgScore[kv.first] = 0.0;
            }
            // taxon2avgScore[kv.first] = kv.second / cladeCounts[kv.first].cladeCount;
        }
        taxon2avgScore_copy = taxon2avgScore;
    }
   
    
    std::vector<uint64_t> filteredQueryIndices;

    taxon2avgScore.clear();
    taxonCounts.clear();
    cladeCounts.clear();
    
    // Re-read the result file and apply filtering criteria
    resultFile.clear();
    resultFile.seekg(0, std::ios::beg);
    uint64_t queryIndex = 0;
    while (std::getline(resultFile, line)) {
        if (line[0] == '#') {
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
            for (size_t i = 0; i < columns.size(); i++) {
                if (columns[i] == "taxID") { taxidCol = i; }
                if (columns[i] == "score") { scoreCol = i; }
                if (columns[i] == "e_value") { eValueCol = i; }
                if (columns[i] == "lineage") { lineageCol = i; }
            }
            queryIndex++;
            continue;
        }

        if (line.empty()) {
            queryIndex++;
            continue;
        }

        totalReads++;

        items = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
        TaxID taxID = ex2inTaxId[std::stoi(items[taxidCol])];

        if (taxID != 0) {
            if (par.maxEValue > 0 && eValueCol != -1 && std::stod(items[eValueCol]) > par.maxEValue) {
                filteredQueryIndices.push_back(queryIndex);
                taxonCounts[0]++;
                queryIndex++;
                continue;
            }

            if (par.minScore > 0 && scoreCol != -1 && std::stod(items[scoreCol]) < par.minScore) {
                filteredQueryIndices.push_back(queryIndex);
                taxonCounts[0]++;
                queryIndex++;
                continue;
            }

            if (par.minAvgScore > 0 && scoreCol != -1) {
                TaxID mappableTaxId = taxID;
                if (speciesAndLowerTaxIds.find(taxID) != speciesAndLowerTaxIds.end()) {
                    mappableTaxId = taxonomy->getTaxIdAtRank(taxID, "species");
                }
                if (taxon2avgScore_copy[mappableTaxId] < par.minAvgScore) {
                    filteredQueryIndices.push_back(queryIndex);
                    taxonCounts[0]++;
                    queryIndex++;
                    continue;
                }
            }

            if (denylist.size() > 0) {
                if (allDenyTaxIds.find(taxID) != allDenyTaxIds.end()) {
                    filteredQueryIndices.push_back(queryIndex);
                    taxonCounts[0]++;
                    queryIndex++;
                    continue;
                }
            }

            // If it passes all filters, keep the classification and add the score to taxon2avgScore
            double score = std::stod(items[scoreCol]);    
            TaxID mappableTaxId = taxID;
            if (speciesAndLowerTaxIds.find(taxID) != speciesAndLowerTaxIds.end()) {
                mappableTaxId = taxonomy->getTaxIdAtRank(taxID, "species");
            }
            taxon2avgScore[mappableTaxId] += score;
            taxonCounts[taxID]++;
            totalClassifiedReads++;
        }

        queryIndex++;
    }

    // Recalculate clade counts and average score after filtering
    cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);
    Classifier::rollUpScore(taxon2avgScore, parentToChildren, 1);
    for (const auto& kv : taxon2avgScore) {
        auto it = cladeCounts.find(kv.first);
        if (it != cladeCounts.end() && it->second.cladeCount > 0) {
            taxon2avgScore[kv.first] = kv.second / static_cast<double>(it->second.cladeCount);
        } else {
            taxon2avgScore[kv.first] = 0.0;
        }
        // taxon2avgScore[kv.first] = kv.second / cladeCounts[kv.first].cladeCount;
    }

    // Filter out taxa based on clade count and proportion, and get the set of tax IDs to clean (set to unclassified)
    std::unordered_set<TaxID> taxIdToClean;
    if (par.minCladeCount > 0 || par.minCladeProportion > 0) {
        // Get species to filter out based on clade count and proportion
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

        taxIdToClean = getSelfAndChildren(rejectedSpeciesIds, parentToChildren);
        for (TaxID taxId : taxIdToClean) {
            taxonCounts.erase(taxId);
        }

        cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);
    }    

    reporter.writeReportFile(totalReads, cladeCounts, taxon2avgScore, ReportType::Default);

    if (!par.outFilteredResults.empty()) {
        std::ofstream filteredResultsFile(par.outFilteredResults);
        if (!filteredResultsFile.is_open()) {
            std::cerr << "Error: could not open output file for filtered results: " << par.outFilteredResults << std::endl;
            return 1;
        }

        // Re-read the result file and write the filtered results to the output file
        resultFile.clear();
        resultFile.seekg(0, std::ios::beg);
        queryIndex = 0;
        uint64_t filteredIndex = 0;
        while (std::getline(resultFile, line)) {
            // Print as is for header lines, empty lines, and unclassified lines (taxID 0)
            if (line[0] == '#' || line.empty() || line[0] == '0') {
                filteredResultsFile << line << "\n";
                queryIndex++;
                continue;
            }

            // Print filtered queries as unclassified (taxID 0)
            if (filteredIndex < filteredQueryIndices.size() && queryIndex == filteredQueryIndices[filteredIndex]) {
                std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
                filteredResultsFile 
                    << "0\t"                    // 0 for unclassified
                    << columns[1]               // query name
                    << "\t0\t"                  // taxID 0 for unclassified         
                    << columns[3] << "\t"       // query length
                    << "0\t"                    // score
                    << "-\t"                    // eValue
                    << "-\t";                   // rank
                if (lineageCol != -1) {
                    filteredResultsFile << "-\t";
                }
                filteredResultsFile << "-\t\n"; // taxID:match_count

                filteredIndex++;
                queryIndex++;
                continue;
            } 


            // Print low clade count/proportion queries as unclassified (taxID 0)
            if (taxIdToClean.size() > 0) {
                std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
                TaxID taxID = ex2inTaxId[std::stoi(columns[taxidCol])];
                if (taxIdToClean.find(taxID) != taxIdToClean.end()) {
                    filteredResultsFile 
                        << "0\t"                    // 0 for unclassified
                        << columns[1]               // query name
                        << "\t0\t"                  // taxID 0 for unclassified         
                        << columns[3] << "\t"       // query length
                        << "0\t"                    // score
                        << "-\t"                    // eValue
                        << "-\t";                   // rank
                    if (lineageCol != -1) {
                        filteredResultsFile << "-\t";
                    }
                    filteredResultsFile << "-\t\n"; // taxID:match_count
                    queryIndex++;
                    continue;
                }
                
            }

            // Print the rest of the lines as is
            filteredResultsFile << line << "\n";
            queryIndex++;
        }
    }

    return 0;
}


std::unordered_set<TaxID> getSelfAndChildren(
    const std::vector<TaxID>& inputTaxIds,
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren) 
{
    std::unordered_set<TaxID> result;
    std::vector<TaxID> stack = inputTaxIds;

    // 1. Add the input tax IDs to the result set
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

