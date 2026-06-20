#include "ClassificationResultRefiner.h"

#include "Classifier.h"
#include "Debug.h"
#include "FileUtil.h"
#include "Reporter.h"

#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>

#ifdef OPENMP
#include "omp.h"
#endif

using namespace std;

ClassificationResultRefiner::ClassificationResultRefiner(const LocalParameters &par) : par(par) {}

ClassificationResultRefiner::ParsedClassification ClassificationResultRefiner::parseFields(
    const std::vector<std::string> &fields,
    bool hasEvalue,
    bool hasLineage)
{
    ParsedClassification result;

    if (fields.size() < 6) {
        return result;
    }

    result.isClassified = (fields[0] == "1");
    result.readId = fields[1];
    result.taxonomyId = std::stoi(fields[2]);
    result.effectiveReadLength = std::stoi(fields[3]);
    result.dnaIdentityScore = std::stof(fields[4]);

    int rankIdx = hasEvalue ? 6 : 5;
    int evalueIdx = hasEvalue ? 5 : -1;
    if (evalueIdx == 5) {
        if (fields[evalueIdx] == "-") {
            result.evalue = -1.0f;
        } else {
            result.evalue = std::stod(fields[evalueIdx]);
        }
    }
    result.classificationRank = fields[rankIdx];

    if (hasLineage) {
        result.fullLineage = fields[rankIdx + 1];
        result.taxIdKmerCounts = fields[rankIdx + 2];
    } else {
        result.fullLineage = "-";
        result.taxIdKmerCounts = fields[rankIdx + 1];
    }

    return result;
}

bool ClassificationResultRefiner::isInAnyClade(
    TaxonomyWrapper *taxonomy,
    const std::vector<TaxID> &cladeTaxIds,
    TaxID taxId)
{
    for (TaxID cladeTaxId : cladeTaxIds) {
        if (taxonomy->IsAncestor(cladeTaxId, taxId)) {
            return true;
        }
    }
    return false;
}

std::unordered_set<TaxID> ClassificationResultRefiner::getSelfAndChildren(
    const std::vector<TaxID> &taxIds,
    const std::unordered_map<TaxID, std::vector<TaxID>> &parentToChildren)
{
    std::unordered_set<TaxID> result;
    std::vector<TaxID> stack = taxIds;

    for (TaxID taxId : stack) {
        result.insert(taxId);
    }

    while (!stack.empty()) {
        TaxID currentTaxId = stack.back();
        stack.pop_back();

        auto it = parentToChildren.find(currentTaxId);
        if (it != parentToChildren.end()) {
            for (TaxID childId : it->second) {
                if (result.insert(childId).second) {
                    stack.push_back(childId);
                }
            }
        }
    }

    return result;
}

int ClassificationResultRefiner::runRefineReport(
    const std::string &resultFileName,
    const std::string &dbDir,
    const std::string &outFileName) const
{
    std::unique_ptr<TaxonomyWrapper> taxonomy(loadTaxonomy(dbDir, dbDir + "/taxonomy/"));
    std::unordered_map<TaxID, TaxID> ex2inTaxId;
    taxonomy->getExternal2internalTaxID(ex2inTaxId);
    const std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();

    Reporter reporter(par, taxonomy.get(), outFileName);

    std::vector<TaxID> denylist;
    std::unordered_set<TaxID> allDenyTaxIds;
    if (!par.excludeTaxid.empty()) {
        std::vector<std::string> contams = Util::split(par.excludeTaxid, ",");
        for (const std::string &contam : contams) {
            TaxID externalTaxId = std::stoi(contam);
            if (ex2inTaxId.find(externalTaxId) == ex2inTaxId.end()) {
                std::cerr << "Warning: tax ID " << externalTaxId << " in the exclude list does not exist in the taxonomy and will be ignored." << std::endl;
                continue;
            }
            denylist.push_back(ex2inTaxId[externalTaxId]);
        }
        allDenyTaxIds = getSelfAndChildren(denylist, parentToChildren);
    }

    std::vector<TaxID> speciesTaxIds;
    int maxTaxId = taxonomy->getMaxTaxID();
    for (int i = 0; i <= maxTaxId; i++) {
        const char *rank = taxonomy->getString(taxonomy->taxonNodes[i].rankIdx);
        if (strcmp(rank, "species") == 0) {
            speciesTaxIds.push_back(taxonomy->taxonNodes[i].taxId);
        }
    }

    std::unordered_set<TaxID> speciesAndLowerTaxIds = getSelfAndChildren(speciesTaxIds, parentToChildren);

    unordered_map<TaxID, TaxonCounts> cladeCounts;
    std::unordered_map<TaxID, unsigned int> taxonCounts;
    std::unordered_map<TaxID, double> taxon2avgScore;
    std::unordered_map<TaxID, double> taxon2avgScoreCopy;
    size_t totalReads = 0;
    size_t totalClassifiedReads = 0;
    std::string line;
    int taxidCol = -1;
    int scoreCol = -1;
    int eValueCol = -1;
    int lineageCol = -1;
    std::vector<std::string> items;

    std::ifstream resultFile(resultFileName);
    if (!resultFile.is_open()) {
        std::cerr << "Error: could not open result file " << resultFileName << std::endl;
        return 1;
    }

    if (par.minAvgScore > 0) {
        while (std::getline(resultFile, line)) {
            if (line.empty()) {
                continue;
            }
            if (line[0] == '#') {
                std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
                for (size_t i = 0; i < columns.size(); i++) {
                    if (columns[i] == "taxID") { taxidCol = static_cast<int>(i); }
                    if (columns[i] == "score") { scoreCol = static_cast<int>(i); }
                }
                continue;
            }

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
        for (const auto &kv : taxon2avgScore) {
            auto it = cladeCounts.find(kv.first);
            if (it != cladeCounts.end() && it->second.cladeCount > 0) {
                taxon2avgScore[kv.first] = kv.second / static_cast<double>(it->second.cladeCount);
            } else {
                taxon2avgScore[kv.first] = 0.0;
            }
        }
        taxon2avgScoreCopy = taxon2avgScore;
    }

    std::vector<uint64_t> filteredQueryIndices;

    taxon2avgScore.clear();
    taxonCounts.clear();
    cladeCounts.clear();

    resultFile.clear();
    resultFile.seekg(0, std::ios::beg);
    uint64_t queryIndex = 0;
    while (std::getline(resultFile, line)) {
        if (line.empty()) {
            queryIndex++;
            continue;
        }
        if (line[0] == '#') {
            std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
            for (size_t i = 0; i < columns.size(); i++) {
                if (columns[i] == "taxID") { taxidCol = static_cast<int>(i); }
                if (columns[i] == "score") { scoreCol = static_cast<int>(i); }
                if (columns[i] == "e_value") { eValueCol = static_cast<int>(i); }
                if (columns[i] == "lineage") { lineageCol = static_cast<int>(i); }
            }
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
                if (taxon2avgScoreCopy[mappableTaxId] < par.minAvgScore) {
                    filteredQueryIndices.push_back(queryIndex);
                    taxonCounts[0]++;
                    queryIndex++;
                    continue;
                }
            }

            if (!denylist.empty() && allDenyTaxIds.find(taxID) != allDenyTaxIds.end()) {
                filteredQueryIndices.push_back(queryIndex);
                taxonCounts[0]++;
                queryIndex++;
                continue;
            }

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

    cladeCounts = taxonomy->getCladeCounts(taxonCounts, parentToChildren);
    Classifier::rollUpScore(taxon2avgScore, parentToChildren, 1);
    for (const auto &kv : taxon2avgScore) {
        auto it = cladeCounts.find(kv.first);
        if (it != cladeCounts.end() && it->second.cladeCount > 0) {
            taxon2avgScore[kv.first] = kv.second / static_cast<double>(it->second.cladeCount);
        } else {
            taxon2avgScore[kv.first] = 0.0;
        }
    }

    std::unordered_set<TaxID> taxIdToClean;
    if (par.minCladeCount > 0 || par.minCladeProportion > 0) {
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

        resultFile.clear();
        resultFile.seekg(0, std::ios::beg);
        queryIndex = 0;
        uint64_t filteredIndex = 0;
        while (std::getline(resultFile, line)) {
            if (line.empty() || line[0] == '#' || line[0] == '0') {
                filteredResultsFile << line << "\n";
                queryIndex++;
                continue;
            }

            if (filteredIndex < filteredQueryIndices.size() && queryIndex == filteredQueryIndices[filteredIndex]) {
                std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
                filteredResultsFile
                    << "0\t"
                    << columns[1]
                    << "\t0\t"
                    << columns[3] << "\t"
                    << "0\t"
                    << "-\t"
                    << "-\t";
                if (lineageCol != -1) {
                    filteredResultsFile << "-\t";
                }
                filteredResultsFile << "-\t\n";

                filteredIndex++;
                queryIndex++;
                continue;
            }

            if (!taxIdToClean.empty()) {
                std::vector<std::string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 10);
                TaxID taxID = ex2inTaxId[std::stoi(columns[taxidCol])];
                if (taxIdToClean.find(taxID) != taxIdToClean.end()) {
                    filteredResultsFile
                        << "0\t"
                        << columns[1]
                        << "\t0\t"
                        << columns[3] << "\t"
                        << "0\t"
                        << "-\t"
                        << "-\t";
                    if (lineageCol != -1) {
                        filteredResultsFile << "-\t";
                    }
                    filteredResultsFile << "-\t\n";
                    queryIndex++;
                    continue;
                }
            }

            filteredResultsFile << line << "\n";
            queryIndex++;
        }
    }

    return 0;
}

int ClassificationResultRefiner::runRefineResult(const std::string &classifiedFile) const
{
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    const string &dbDir = par.filenames[1];
    const string &taxonomyDir = par.taxonomyPath;
    std::unique_ptr<TaxonomyWrapper> taxonomy(loadTaxonomy(dbDir, taxonomyDir));
    unordered_map<TaxID, TaxID> extern2intern;
    taxonomy->getExternal2internalTaxID(extern2intern);
    bool useEvalue = (par.maxEValue > 0);

    string reportFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined_report.tsv";
    Reporter reporter(par, taxonomy.get(), reportFileName);

    string refinedFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined.tsv";
    cout << "Write refined classification results to: " << endl;
    cout << refinedFileName << endl;

    if (FileUtil::fileExists(refinedFileName.c_str())) {
        Debug(Debug::INFO) << refinedFileName << " is already exists.\n";
        return 0;
    }
    ofstream refinedFile(refinedFileName.c_str());

    if (!refinedFile.is_open()) {
        Debug(Debug::ERROR) << "Could not open " << refinedFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }

    refinedFile.close();

    vector<string> contams;
    vector<TaxID> contamsTaxIds;
    if (!par.excludeTaxid.empty()) {
        contams = Util::split(par.excludeTaxid, ",");
        for (const string &contam : contams) {
            contamsTaxIds.push_back(extern2intern[stoi(contam)]);
        }
    }

    vector<string> targets;
    vector<TaxID> targetsTaxIds;
    if (!par.selectTaxid.empty()) {
        targets = Util::split(par.selectTaxid, ",");
        for (const string &target : targets) {
            targetsTaxIds.push_back(extern2intern[stoi(target)]);
            const TaxID targetId = extern2intern[stoi(target)];

            if (isInAnyClade(taxonomy.get(), contamsTaxIds, targetId)) {
                Debug(Debug::ERROR) << "Excluded taxid is selected : " << target << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    }

    vector<string> columns;
    vector<int> columnsIdx;
    if (!par.selectColumns.empty()) {
        columns = Util::split(par.selectColumns, ",");
        for (const auto &column : columns) {
            columnsIdx.push_back(stoi(column));
        }
    }

    std::string criterionRank = par.rank;
    int createUpperRanksFile = par.higherRankFile;

    if (!criterionRank.empty() && taxonomy->findRankIndex(criterionRank) == -1) {
        Debug(Debug::ERROR) << "Invalid criterion rank: " << criterionRank << ". Rank not found in NcbiRanks.\n";
        EXIT(EXIT_FAILURE);
    }

    std::string upperRankFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined_higherRanks.tsv";

    if (createUpperRanksFile == 2) {
        cout << "Write higher rank reads to: " << endl;
        cout << upperRankFileName << endl;
        ofstream upperRankFile(upperRankFileName.c_str());
        if (!upperRankFile.is_open()) {
            Debug(Debug::ERROR) << "Could not open " << upperRankFileName << " for writing\n";
            EXIT(EXIT_FAILURE);
        }
        upperRankFile.close();
    }

    ifstream file(classifiedFile);
    ofstream refinedFileAppend(refinedFileName.c_str(), ios::app);
    if (!file.is_open() || !refinedFileAppend.is_open()) {
        Debug(Debug::ERROR) << "Could not open input or output file\n";
        EXIT(EXIT_FAILURE);
    }

    std::string firstLine;
    std::getline(file, firstLine);
    bool hasEvalue = (firstLine.find("e_value") != std::string::npos);
    bool hasLineage = (firstLine.find("lineage") != std::string::npos);

    if (firstLine[0] != '#') {
        file.seekg(0);
    }

    ofstream upperRankFileAppend;
    if (createUpperRanksFile == 2) {
        upperRankFileAppend.open(upperRankFileName.c_str(), ios::app);
        if (!upperRankFileAppend.is_open()) {
            Debug(Debug::ERROR) << "Could not open " << upperRankFileName << " for writing\n";
            EXIT(EXIT_FAILURE);
        }
    }

    const size_t CHUNK_SIZE = 1000;

    size_t totalSeqCnt = 0;
    std::unordered_map<TaxID, unsigned int> taxcntSum;
    std::vector<std::string> resultChunk;
    std::vector<std::string> upperRanks;
    const bool printLineage = par.printLineage;
    const bool hasSelectedColumns = !par.selectColumns.empty();
    const bool removeUnclassified = par.removeUnclassified;
    const float minScore = par.minScore;
    const double maxEValue = par.maxEValue;

#pragma omp parallel default(none) shared(file, refinedFileAppend, cout, taxonomy, \
    extern2intern, contamsTaxIds, targetsTaxIds, columnsIdx, totalSeqCnt, \
    resultChunk, taxcntSum, upperRanks, upperRankFileAppend, \
    createUpperRanksFile, criterionRank, hasEvalue, hasLineage, useEvalue, \
    printLineage, hasSelectedColumns, removeUnclassified, minScore, maxEValue)
    {
#pragma omp single
        {
#ifdef OPENMP
            resultChunk.resize(10000 * omp_get_num_threads());
            upperRanks.resize(10000 * omp_get_num_threads());
#else
            resultChunk.resize(10000);
            upperRanks.resize(10000);
#endif
            vector<string> chunk;
            string line;
            size_t chunkCnt = 0;

            while (true) {
                chunk.clear();
                chunkCnt++;
                for (size_t i = 0; i < CHUNK_SIZE && getline(file, line); ++i) {
                    totalSeqCnt++;
                    chunk.push_back(line);
                }

                if (chunkCnt > (size_t) 10 * omp_get_num_threads()) {
#pragma omp taskwait
                    for (const string &line : resultChunk) {
                        refinedFileAppend << line;
                    }

                    resultChunk.clear();
                    resultChunk.resize(10000 * omp_get_num_threads());

                    if (createUpperRanksFile == 2) {
                        for (const string &line : upperRanks) {
                            if (!line.empty()) {
                                upperRankFileAppend << line;
                            }
                        }

                        upperRanks.clear();
                        upperRanks.resize(10000 * omp_get_num_threads());
                    }

                    chunkCnt = 1;
                }

                if (chunk.empty()) {
#pragma omp taskwait
                    for (const string &line : resultChunk) {
                        refinedFileAppend << line;
                    }
                    resultChunk.clear();

                    if (createUpperRanksFile == 2) {
                        for (const string &line : upperRanks) {
                            if (!line.empty()) {
                                upperRankFileAppend << line;
                            }
                        }

                        upperRanks.clear();
                    }
                    break;
                }

#pragma omp task firstprivate(chunk, chunkCnt)
                {
                    std::unordered_map<TaxID, unsigned int> taxCounts;
                    for (size_t localIdx = 0; localIdx < chunk.size(); ++localIdx) {
                        const string &localLine = chunk[localIdx];
                        vector<string> fields = Util::split(localLine, "\t");

                        size_t matchCountIdx = 6;
                        matchCountIdx += hasEvalue ? 1 : 0;
                        matchCountIdx += hasLineage ? 1 : 0;

                        while (fields.size() <= matchCountIdx) {
                            fields.push_back("-");
                        }

                        ParsedClassification data = parseFields(fields, hasEvalue, hasLineage);
                        data.taxonomyId = extern2intern[data.taxonomyId];

                        if (printLineage && !hasLineage) {
                            if (data.isClassified) {
                                data.fullLineage = taxonomy->taxLineage2(taxonomy->taxonNode(data.taxonomyId));
                            } else {
                                data.fullLineage = "-";
                            }
                            fields.push_back(data.fullLineage);
                        }

                        std::vector<int> printCols;
                        if (hasSelectedColumns) {
                            printCols = columnsIdx;
                        } else {
                            for (size_t i = 0; i < fields.size(); ++i) {
                                printCols.push_back(static_cast<int>(i));
                            }
                        }

                        if (!(removeUnclassified == true && data.isClassified == false) &&
                            !(contamsTaxIds.size() > 0 && isInAnyClade(taxonomy.get(), contamsTaxIds, data.taxonomyId)) &&
                            !(targetsTaxIds.size() > 0 && !isInAnyClade(taxonomy.get(), targetsTaxIds, data.taxonomyId)) &&
                            !(data.isClassified == true && minScore > data.dnaIdentityScore) &&
                            !(data.isClassified == true && hasEvalue && useEvalue && data.evalue > maxEValue)) {

                            stringstream ss;
                            stringstream tt;
                            if (!criterionRank.empty()) {
                                int criterion = taxonomy->findRankIndex(criterionRank);
                                if (taxonomy->findRankIndex(data.classificationRank) <= criterion) {
                                    data.taxonomyId = taxonomy->getTaxIdAtRank(data.taxonomyId, criterionRank);
                                    fields[2] = std::to_string(data.taxonomyId);
                                    data.classificationRank = criterionRank;
                                    fields[hasEvalue ? 6 : 5] = data.classificationRank;
                                } else {
                                    if (createUpperRanksFile == 0) { continue; }
                                    if (createUpperRanksFile == 2) {
                                        for (size_t i = 0; i < printCols.size(); i++) {
                                            tt << fields[printCols[i]] << "\t";
                                        }
                                        tt << endl;

                                        size_t globalIdx = (chunkCnt - 1) * CHUNK_SIZE + localIdx;
                                        upperRanks[globalIdx] = tt.str();
                                        continue;
                                    }
                                }
                            }

                            for (size_t i = 0; i < printCols.size(); i++) {
                                cout << fields[printCols[i]] << "\t";
                                if (printCols[i] < static_cast<int>(fields.size())) {
                                    ss << fields[printCols[i]] << "\t";
                                }
                            }
                            if (printLineage && !hasLineage) {
                                cout << data.fullLineage << "\t";
                                ss << data.fullLineage << "\t";
                            }

                            cout << endl;
                            ss << endl;
                            ++taxCounts[data.taxonomyId];

                            size_t globalIdx = (chunkCnt - 1) * CHUNK_SIZE + localIdx;
                            resultChunk[globalIdx] = ss.str();
                        }
                    }

#pragma omp critical
                    {
                        for (const auto &it : taxCounts) {
                            taxcntSum[it.first] += it.second;
                        }
                    }
                }
            }
        }
    }

    refinedFileAppend.close();
    if (createUpperRanksFile == 2) {
        upperRankFileAppend.close();
    }

    if (par.report) {
        cout << "Write report to: " << endl;
        cout << reportFileName << endl;

        std::string kronaFileName = classifiedFile.substr(0, classifiedFile.find_last_of('.')) + "_refined_krona.html";

        std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
        unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxcntSum, parentToChildren);
        reporter.writeReportFile(totalSeqCnt + 1, cladeCounts, ReportType::Default, kronaFileName);
    }

    return 0;
}
