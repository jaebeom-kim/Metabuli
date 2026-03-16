#include "Reporter.h"
#include "taxonomyreport.cpp"

Reporter::Reporter(const LocalParameters &par, TaxonomyWrapper *taxonomy, const std::string &customReportFileName) : par(par), taxonomy(taxonomy) {
    if (!customReportFileName.empty()){
        reportFileName = customReportFileName;
    } else {
        if (par.targetTaxId != 0) {return;}
        if (par.contamList == "") { // classify module
            if (par.seqMode == 2) {
                outDir = par.filenames[3];
                jobId = par.filenames[4];
            } else {
                outDir = par.filenames[2];
                jobId = par.filenames[3];
            }
            // Output file names
            reportFileName = outDir + + "/" + jobId + "_report.tsv";
            readClassificationFileName = outDir + "/" + jobId + "_classifications.tsv";
            if (par.em) {
                reportFileName_em            = outDir + "/" + jobId + "_EM_report.tsv";
                reportFileName_em_reclassify = outDir + "/" + jobId + "_EM+reclassify_report.tsv";
                reclassifyFileName           = outDir + "/" + jobId + "_EM+reclassify_results.tsv";
                mappingResFileName           = outDir + "/" + jobId + "_mapping_results.txt";
                mappingResBuffer = new WriteBuffer<MappingRes>(mappingResFileName, 1000000);
            }
        }
    }    
}

void Reporter::openReadClassificationFile() {
    readClassificationFile.open(readClassificationFileName);
}

void Reporter::writeReadClassification(const vector<Query> & queryList, bool classifiedOnly) {
    if (isFirstTime) {
        readClassificationFile << "#is_classified\tname\ttaxID\tquery_length\tscore\trank";
        if (par.printLineage) {
            readClassificationFile << "\tlineage";
        }
        readClassificationFile << "\ttaxID:match_count\n";
        isFirstTime = false;
    }
    for (size_t i = 0; i < queryList.size(); i++) {
        if (classifiedOnly && !(queryList[i].classification == 0)) {
            continue;
        }
        if (queryList[i].name.empty()) {
            break;
        }
        if (queryList[i].classification != 0) {
            readClassificationFile 
                << "1\t" 
                << queryList[i].name << "\t"
                << taxonomy->getOriginalTaxID(queryList[i].classification) << "\t"
                << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                << queryList[i].idScore << "\t"
                << queryList[i].subScore << "\t"
                << queryList[i].eValue << "\t"
                << taxonomy->getString(taxonomy->taxonNode(queryList[i].classification)->rankIdx) << "\t";
            
            if (par.printLineage) {
                readClassificationFile << taxonomy->taxLineage2(taxonomy->taxonNode(queryList[i].classification)) << "\t";
            }
            
            for (auto it = queryList[i].taxCnt.begin(); it != queryList[i].taxCnt.end(); ++it) {
                readClassificationFile << taxonomy->getOriginalTaxID(it->first) << ":" << it->second << " ";
            }
            readClassificationFile << "\n";
        } else {
            readClassificationFile 
                << "0\t" 
                << queryList[i].name << "\t"
                << taxonomy->getOriginalTaxID(queryList[i].classification) << "\t"
                << queryList[i].queryLength + queryList[i].queryLength2 << "\t"
                << queryList[i].idScore << "\t"
                << "-" << "\t" // subScore
                << "-" << "\t" // eValue
                << "-" << "\t";
            
            if (par.printLineage) {
                readClassificationFile << "-\t";
            }
            readClassificationFile << "-\t\n";
        }
    }
}

void Reporter::closeReadClassificationFile() {
    readClassificationFile.close();
}

void Reporter::kronaReport(FILE *FP, const TaxonomyWrapper &taxDB, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts, unsigned long totalReads, TaxID taxID, int depth) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "<node name=\"unclassified\"><magnitude><val>%d</val></magnitude></node>", cladeCount);
        }
        kronaReport(FP, taxDB, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxDB.taxonNode(taxID);
        std::string escapedName = escapeAttribute(taxDB.getString(taxon->nameIdx));
        fprintf(FP, "<node name=\"%s\"><magnitude><val>%d</val></magnitude>", escapedName.c_str(), cladeCount);
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                kronaReport(FP, taxDB, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
        fprintf(FP, "</node>");
    }
}

void Reporter::writeReportFile(
    int numOfQuery, 
    unordered_map<TaxID, unsigned int> &taxCnt, 
    ReportType reportType,
    string kronaFileName) 
{
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
    FILE *fp = nullptr;
    if (reportType == ReportType::Default) {
        fp = fopen(reportFileName.c_str(), "w");
    } else if (reportType == ReportType::EM) {
        fp = fopen(reportFileName_em.c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        fp = fopen(reportFileName_em_reclassify.c_str(), "w");
    }
    fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n");
    writeReport(fp, cladeCounts, numOfQuery);
    fclose(fp);

    // Write Krona chart
    if (jobId.empty()) { return; }
    
    FILE *kronaFile = nullptr;

    if (reportType == ReportType::Default) {
        if (!kronaFileName.empty()) {
            kronaFile = fopen(kronaFileName.c_str(), "w");
        } else {
            kronaFile = fopen((outDir + "/" + jobId + "_krona.html").c_str(), "w");
        }
    } else if (reportType == ReportType::EM) {
        kronaFile = fopen((outDir + "/" + jobId + "_EM_krona.html").c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        kronaFile = fopen((outDir + "/" + jobId + "_EM+reclassify_krona.html").c_str(), "w");
    }
    if (kronaFile == nullptr) {
        Debug(Debug::ERROR) << "Could not open Krona file for writing: " << kronaFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), kronaFile);
    fprintf(kronaFile, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", (size_t) numOfQuery);
    kronaReport(kronaFile, *taxonomy, cladeCounts, numOfQuery);
    fprintf(kronaFile, "</node></krona></div></body></html>");
    fclose(kronaFile);
}

void Reporter::writeReport(FILE *FP, const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
                             unsigned long totalReads, TaxID taxID, int depth) {
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        writeReport(FP, cladeCounts, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);
        fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t%s%s\n",
                100 * cladeCount / double(totalReads), cladeCount, taxCount,
                taxonomy->getString(taxon->rankIdx), taxonomy->getOriginalTaxID(taxID), std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                writeReport(FP, cladeCounts, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

// void Reporter::writeEMreportFile(
//     unordered_map<TaxID, double> &taxProbs) 
// {
//     std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
//     unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeProbs(taxCnt, parentToChildren);
//     FILE *fp = fopen(reportFileName_em.c_str(), "w");
//     fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tname\n");
//     writeReport(fp, cladeCounts, numOfQuery);
//     fclose(fp);

//     // Write Krona chart
//     if (jobId.empty()) {
//         return;
//     }
//     FILE *kronaFile = nullptr;
//     if (!em) { 
//         if (!kronaFileName.empty()) {
//             kronaFile = fopen(kronaFileName.c_str(), "w");
//         } else {
//             kronaFile = fopen((outDir + "/" + jobId + "_krona.html").c_str(), "w");
//         }
//     } else {
//         kronaFile = fopen((outDir + "/" + jobId + "_EM_krona.html").c_str(), "w");
//     }
//     fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), kronaFile);
//     fprintf(kronaFile, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", (size_t) numOfQuery);
//     kronaReport(kronaFile, *taxonomy, cladeCounts, numOfQuery);
//     fprintf(kronaFile, "</node></krona></div></body></html>");
//     fclose(kronaFile);
// }


unsigned int Reporter::cladeCountVal(const std::unordered_map<TaxID, TaxonCounts> &map, TaxID key) {
    typename std::unordered_map<TaxID, TaxonCounts>::const_iterator it = map.find(key);
    if (it == map.end()) {
        return 0;
    } else {
        return it->second.cladeCount;
    }
}

void Reporter::getReadsClassifiedToClade(TaxID cladeId,
                                         const string &readClassificationFileName,
                                         vector<size_t> &readIdxs) {
    FILE *results = fopen(readClassificationFileName.c_str(), "r");
    if (!results) {
        perror("Failed to open read-by-read classification file");
        return;
    }
    char line[4096];
    size_t idx = 0;
    if (cladeId == -1) {
        while (fgets(line, sizeof(line), results)) {
            if (line[0] == '#') {
                continue;
            }
            if (line[0] == '0') { // unclassified
                readIdxs.push_back(idx);
            }
            idx++;
        }
    } else if (taxonomy->hasInternalTaxID()) {
        unordered_map<TaxID, TaxID> extern2intern;
        taxonomy->getExternal2internalTaxID(extern2intern);
        while (fgets(line, sizeof(line), results)) {
            if (line[0] == '#') {
                continue;
            }
            int taxId;
            if (sscanf(line, "%*s %*s %d", &taxId) == 1) {            
                if (taxonomy->IsAncestor(cladeId, extern2intern[taxId])) {
                    readIdxs.push_back(idx);
                }
            }
            idx++;
        }
    } else {
        while (fgets(line, sizeof(line), results)) {
            if (line[0] == '#') {
                continue;
            }
            int taxId;
            if (sscanf(line, "%*s %*s %d", &taxId) == 1) {            
                if (taxonomy->IsAncestor(cladeId, taxId)) {
                    readIdxs.push_back(idx);
                }
            }
            idx++;
        }
    }
    fclose(results);
}

void Reporter::printSpecifiedReads(const vector<size_t> & readIdxs,
                                   const string & readFileName,
                                   string & outFileName) {
    // Check FASTA or FASTQ
    KSeqWrapper* tempKseq = KSeqFactory(readFileName.c_str());
    tempKseq->ReadEntry();
    bool isFasta = tempKseq->entry.qual.l == 0;

    if (isFasta && par.extractMode == 2) {
        Debug(Debug::ERROR) << "Cannot convert FASTA to FASTQ\n";
        EXIT(EXIT_FAILURE);
    }

    delete tempKseq;

    bool printFasta;
    if (isFasta || par.extractMode == 1) {
        printFasta = true;
        outFileName += ".fna";
    } else {
        printFasta = false;
        outFileName += ".fq";
    }
    
    KSeqWrapper* kseq = KSeqFactory(readFileName.c_str());
    FILE *outFile = fopen(outFileName.c_str(), "w");
    if (!outFile) {
        perror("Failed to open file");
        return;
    }

    size_t readCnt = 0;
    size_t idx = 0;

    if (printFasta) {
        while (kseq->ReadEntry()) {
            if (readCnt == readIdxs[idx]) {
                fprintf(outFile, ">%s", kseq->entry.name.s);
                if (kseq->entry.comment.l > 0) {
                    fprintf(outFile, " %s\n", kseq->entry.comment.s);
                } else {
                    fprintf(outFile, "\n");
                }
                fprintf(outFile, "%s\n", kseq->entry.sequence.s);
                idx++;
                if (idx == readIdxs.size()) {
                    break;
                }
            }
            readCnt++;
        }
    } else {
        while (kseq->ReadEntry()) {
            if (readCnt == readIdxs[idx]) {
                fprintf(outFile, "@%s", kseq->entry.name.s);
                if (kseq->entry.comment.l > 0) {
                    fprintf(outFile, " %s\n", kseq->entry.comment.s);
                } else {
                    fprintf(outFile, "\n");
                }
                fprintf(outFile, "%s\n", kseq->entry.sequence.s);
                fprintf(outFile, "+%s", kseq->entry.name.s);
                if (kseq->entry.comment.l > 0) {
                    fprintf(outFile, " %s\n", kseq->entry.comment.s);
                } else {
                    fprintf(outFile, "\n");
                }
                fprintf(outFile, "%s\n", kseq->entry.qual.s);
                idx++;
                if (idx == readIdxs.size()) {
                    break;
                }
            }
            readCnt++;
        }
    }
    delete kseq;
}

void Reporter::writeReclassifyResults(const std::vector<Classification> & results)
{   
    ofstream emResultFile(reclassifyFileName, std::ios::out | std::ios::trunc);
    if (!emResultFile.is_open()) {
        cerr << "Error: Could not open EM results file " << reclassifyFileName << endl;
        return;
    }

    emResultFile << "#is_classified\tname\ttaxID\tquery_length\tscore\trank";
    if (par.printLineage) {
        emResultFile << "\tlineage";
    }
    emResultFile << endl;

    for (size_t i = 0; i < results.size(); ++i) {
        const Classification &result = results[i];
        if (result.taxId != 0) {
            emResultFile << (result.taxId != 0) << "\t" 
                         << result.name << "\t"
                         << taxonomy->getOriginalTaxID(result.taxId) << "\t"
                         << result.length << "\t"
                         << result.score << "\t"
                         << taxonomy->getString(taxonomy->taxonNode(result.taxId)->rankIdx);
            if (par.printLineage) {
                emResultFile << "\t" << taxonomy->taxLineage2(taxonomy->taxonNode(result.taxId));
            }
        } else {
            emResultFile << (result.taxId != 0) << "\t" 
                         << result.name << "\t"
                         << taxonomy->getOriginalTaxID(result.taxId) << "\t"
                         << result.length << "\t"
                         << result.score << "\t"
                         << "-";
            if (par.printLineage) {
                emResultFile << "\t-";
            }        
        }
        emResultFile << "\n";
    }

    emResultFile.close();
    cout << "EM results written to " << reclassifyFileName << endl;
}

void Reporter::writeReportFile(
    int numOfQuery, 
    unordered_map<TaxID, unsigned int> &taxCnt, 
    const unordered_map<TaxID, double> &species2adjustedEvenness, // NEW PARAMETER
    ReportType reportType,
    string kronaFileName) 
{
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
    FILE *fp = nullptr;
    
    if (reportType == ReportType::Default) {
        fp = fopen(reportFileName.c_str(), "w");
    } else if (reportType == ReportType::EM) {
        fp = fopen(reportFileName_em.c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        fp = fopen(reportFileName_em_reclassify.c_str(), "w");
    }
    
    // UPDATED HEADER: Added 'evenness' column before 'name'
    fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tevenness\tname\n");
    writeReport(fp, cladeCounts, species2adjustedEvenness, numOfQuery); // PASS MAP
    fclose(fp);

    // Write Krona chart (Krona HTML structure left unmodified to prevent rendering breaks)
    if (jobId.empty()) { return; }
    
    FILE *kronaFile = nullptr;

    if (reportType == ReportType::Default) {
        if (!kronaFileName.empty()) {
            kronaFile = fopen(kronaFileName.c_str(), "w");
        } else {
            kronaFile = fopen((outDir + "/" + jobId + "_krona.html").c_str(), "w");
        }
    } else if (reportType == ReportType::EM) {
        kronaFile = fopen((outDir + "/" + jobId + "_EM_krona.html").c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        kronaFile = fopen((outDir + "/" + jobId + "_EM+reclassify_krona.html").c_str(), "w");
    }
    if (kronaFile == nullptr) {
        Debug(Debug::ERROR) << "Could not open Krona file for writing: " << kronaFileName << "\n";
        EXIT(EXIT_FAILURE);
    }
    fwrite(krona_prelude_html, krona_prelude_html_len, sizeof(char), kronaFile);
    fprintf(kronaFile, "<node name=\"all\"><magnitude><val>%zu</val></magnitude>", (size_t) numOfQuery);
    kronaReport(kronaFile, *taxonomy, cladeCounts, numOfQuery);
    fprintf(kronaFile, "</node></krona></div></body></html>");
    fclose(kronaFile);
}


void Reporter::writeReport(
    FILE *FP, 
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
    const std::unordered_map<TaxID, double> &species2adjustedEvenness, 
    unsigned long totalReads, 
    TaxID taxID, 
    int depth) 
{
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    
    if (taxID == 0) {
        if (cladeCount > 0) {
            // Unclassified naturally gets a '-' for evenness
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\t-\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        writeReport(FP, cladeCounts, species2adjustedEvenness, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);
        
        // NEW: Check if this specific taxID has a calculated evenness score
        auto evIt = species2adjustedEvenness.find(taxID);
        if (evIt != species2adjustedEvenness.end()) {
            // Found: Print the score formatted to 4 decimal places
            fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t%.4f\t%s%s\n",
                    100 * cladeCount / double(totalReads), cladeCount, taxCount,
                    taxonomy->getString(taxon->rankIdx), taxonomy->getOriginalTaxID(taxID),
                    evIt->second, // Evenness score
                    std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));
        } else {
            // Not Found (Higher rank/not a species): Print '-' placeholder
            fprintf(FP, "%.4f\t%i\t%i\t%s\t%i\t-\t%s%s\n",
                    100 * cladeCount / double(totalReads), cladeCount, taxCount,
                    taxonomy->getString(taxon->rankIdx), taxonomy->getOriginalTaxID(taxID), 
                    std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));
        }

        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { 
            return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); 
        });
        
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                // Pass the map down the recursive tree
                writeReport(FP, cladeCounts, species2adjustedEvenness, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}


void Reporter::filterClassificationFile(
    const std::string& inputFilePath, 
    const std::string& outputFilePath, 
    const std::unordered_map<TaxID, double>& species2adjustedEvenness, 
    double cutoff) 
{
    std::ifstream inFile(inputFilePath);
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open input file " << inputFilePath << "\n";
        return;
    }

    std::ofstream outFile(outputFilePath);
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open output file " << outputFilePath << "\n";
        return;
    }

    std::unordered_map<TaxID, TaxID> ex2inTaxId;
    taxonomy->getExternal2internalTaxID(ex2inTaxId);
    std::string line;
    while (std::getline(inFile, line)) {
        // 1. Pass the header or empty lines directly to the new file
        if (line.empty() || line[0] == '#') {
            outFile << line << "\n";
            continue;
        }

        // 2. Split the line by tabs to parse the columns
        std::vector<std::string> columns;
        size_t start = 0, end = 0;
        while ((end = line.find('\t', start)) != std::string::npos) {
            columns.push_back(line.substr(start, end - start));
            start = end + 1;
        }
        columns.push_back(line.substr(start)); // Grab the final column

        // Safety check to ensure the line has at least the basic columns
        if (columns.size() < 5) {
            outFile << line << "\n";
            continue;
        }

        // 3. Check if the read is currently classified
        if (columns[0] == "1") {
            TaxID taxID = ex2inTaxId[std::stoull(columns[2])]; // Convert the string taxID to your integer type

            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxID, "species");

            auto evIt = species2adjustedEvenness.find(speciesTaxID);
            
            // 4. If the species is in the map AND its score is below the cutoff, rewrite it
            if (evIt != species2adjustedEvenness.end() && evIt->second < cutoff) {
                // Reconstruct the line exactly as your original "unclassified" else-block did
                outFile << "0\t"          // is_classified
                        << columns[1] << "\t" // name
                        << "0\t"          // taxID (forced to 0)
                        << columns[3] << "\t" // query_length
                        << columns[4] << "\t" // idScore
                        << "-\t"          // subScore
                        << "-\t"          // eValue
                        << "-\t";         // rank

                // Check if lineage was printed (original format has > 9 columns if lineage exists)
                if (columns.size() > 9) {
                    outFile << "-\t";     // empty lineage
                }
                
                outFile << "-\n";         // empty taxID:match_count
                continue;                 // Skip writing the original line
            }
        }

        // 5. If it wasn't filtered, write the original line exactly as it was
        outFile << line << "\n";
    }

    inFile.close();
    outFile.close();
}