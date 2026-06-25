#include "Reporter.h"
#include "taxonomyreport.cpp"
#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <type_traits>

namespace {
inline void appendUnsignedNumber(std::string &out, unsigned long long value) {
    char buffer[32];
    char *ptr = buffer + sizeof(buffer);
    do {
        *--ptr = static_cast<char>('0' + (value % 10));
        value /= 10;
    } while (value != 0);
    out.append(ptr, static_cast<size_t>(buffer + sizeof(buffer) - ptr));
}

template <typename T>
inline void appendNumber(std::string &out, T value, std::true_type) {
    if (value < 0) {
        out += '-';
        appendUnsignedNumber(out, static_cast<unsigned long long>(-(value + 1)) + 1);
    } else {
        appendUnsignedNumber(out, static_cast<unsigned long long>(value));
    }
}

template <typename T>
inline void appendNumber(std::string &out, T value, std::false_type) {
    appendUnsignedNumber(out, static_cast<unsigned long long>(value));
}

template <typename T>
inline void appendNumber(std::string &out, T value) {
    appendNumber(out, value, typename std::is_signed<T>::type());
}

inline void appendFloat(std::string &out, double value) {
    char buffer[64];
    int len = std::snprintf(buffer, sizeof(buffer), "%.6g", value);
    if (len > 0) {
        out.append(buffer, static_cast<size_t>(len));
    }
}

inline void writeMetricOrDash(FILE *fp, double value) {
    if (value >= 0.0) {
        fprintf(fp, "\t%.4f", value);
    } else {
        fprintf(fp, "\t-");
    }
}

inline bool parseTaxId(const std::string &value, TaxID &taxId) {
    errno = 0;
    char *parseEnd = nullptr;
    const long long parsed = std::strtoll(value.c_str(), &parseEnd, 10);
    if (errno == ERANGE
        || parseEnd == value.c_str()
        || *parseEnd != '\0'
        || parsed < std::numeric_limits<TaxID>::min()
        || parsed > std::numeric_limits<TaxID>::max()) {
        return false;
    }

    taxId = static_cast<TaxID>(parsed);
    return true;
}

bool isReadInClade(const char *line,
                   TaxID cladeId,
                   TaxonomyWrapper *taxonomy,
                   const unordered_map<TaxID, TaxID> *extern2intern) {
    if (cladeId == -1) {
        return line[0] == '0';
    }

    if (line[0] == '0') {
        return false;
    }

    int taxId;
    if (sscanf(line, "%*s %*s %d", &taxId) != 1) {
        return false;
    }

    TaxID taxonomyTaxId = taxId;
    if (extern2intern != nullptr) {
        auto it = extern2intern->find(taxId);
        if (it == extern2intern->end()) {
            return false;
        }
        taxonomyTaxId = it->second;
    }

    return taxonomy->IsAncestor(cladeId, taxonomyTaxId);
}

void getReadsByCladeMembership(TaxID cladeId,
                               const string &readClassificationFileName,
                               vector<size_t> &readIdxs,
                               TaxonomyWrapper *taxonomy,
                               bool keepMatches) {
    FILE *results = fopen(readClassificationFileName.c_str(), "r");
    if (!results) {
        perror("Failed to open read-by-read classification file");
        return;
    }

    unordered_map<TaxID, TaxID> extern2intern;
    const unordered_map<TaxID, TaxID> *extern2internPtr = nullptr;
    if (cladeId != -1 && taxonomy->hasInternalTaxID()) {
        taxonomy->getExternal2internalTaxID(extern2intern);
        extern2internPtr = &extern2intern;
    }

    char line[4096];
    size_t idx = 0;
    while (fgets(line, sizeof(line), results)) {
        if (line[0] == '#') {
            continue;
        }

        const bool inClade = isReadInClade(line, cladeId, taxonomy, extern2internPtr);
        if (inClade == keepMatches) {
            readIdxs.push_back(idx);
        }
        idx++;
    }

    fclose(results);
}
}

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
    if (readClassificationFileBuffer.empty()) {
        readClassificationFileBuffer.resize(1 << 20);
    }
    readClassificationFile.rdbuf()->pubsetbuf(readClassificationFileBuffer.data(), static_cast<std::streamsize>(readClassificationFileBuffer.size()));
    readClassificationFile.open(readClassificationFileName);
}

void Reporter::writeReadClassification(const vector<Query> & queryList, bool classifiedOnly) {
    std::string &out = readClassificationOutputBuffer;
    out.clear();
    out.reserve(std::max<size_t>(4096, queryList.size() * 96));

    if (isFirstTime) {
        out += "#is_classified\tname\ttaxID\tquery_length\tscore\te_value\trank";
        if (par.printLineage) {
            out += "\tlineage";
        }
        out += "\ttaxID:match_count\n";
        isFirstTime = false;
    }

    auto originalTaxId = [this](TaxID taxId) {
        auto it = originalTaxIdCache.find(taxId);
        if (it != originalTaxIdCache.end()) {
            return it->second;
        }
        TaxID original = taxonomy->getOriginalTaxID(taxId);
        originalTaxIdCache.emplace(taxId, original);
        return original;
    };

    auto rankName = [this](TaxID taxId) -> const string& {
        auto it = rankCache.find(taxId);
        if (it != rankCache.end()) {
            return it->second;
        }
        const TaxonNode *node = taxonomy->taxonNode(taxId);
        auto inserted = rankCache.emplace(taxId, taxonomy->getString(node->rankIdx));
        return inserted.first->second;
    };

    auto lineage = [this](TaxID taxId) -> const string& {
        auto it = lineageCache.find(taxId);
        if (it != lineageCache.end()) {
            return it->second;
        }
        const TaxonNode *node = taxonomy->taxonNode(taxId);
        auto inserted = lineageCache.emplace(taxId, taxonomy->taxLineage2(node));
        return inserted.first->second;
    };

    for (size_t i = 0; i < queryList.size(); i++) {
        if (classifiedOnly && !(queryList[i].classification == 0)) {
            continue;
        }
        if (queryList[i].name.empty()) {
            break;
        }
        const Query &query = queryList[i];
        if (query.classification != 0) {
            out += "1\t";
            out += query.name;
            out += '\t';
            appendNumber(out, originalTaxId(query.classification));
            out += '\t';
            appendNumber(out, query.queryLength + query.queryLength2);
            out += '\t';
            appendFloat(out, query.idScore);
            out += '\t';
            appendFloat(out, query.eValue);
            out += '\t';
            out += rankName(query.classification);
            out += '\t';

            if (par.printLineage) {
                out += lineage(query.classification);
                out += '\t';
            }

            for (auto it = query.taxCnt.begin(); it != query.taxCnt.end(); ++it) {
                appendNumber(out, originalTaxId(it->first));
                out += ':';
                appendNumber(out, it->second);
                out += ' ';
            }
            out += '\n';
        } else {
            out += "0\t";
            out += query.name;
            out += '\t';
            appendNumber(out, originalTaxId(query.classification));
            out += '\t';
            appendNumber(out, query.queryLength + query.queryLength2);
            out += '\t';
            appendFloat(out, query.idScore);
            out += "\t-\t-";
            out += '\t';

            if (par.printLineage) {
                out += "-\t";
            }
            out += "-\t\n";
        }
    }

    readClassificationFile.write(out.data(), static_cast<std::streamsize>(out.size()));
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
    unordered_map<TaxID, TaxonCounts> cladeCounts,
    const unordered_map<TaxID, double> &taxon2avgScore,
    ReportType reportType,
    string kronaFileName) 
{
    // std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    // unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
    FILE *fp = nullptr;
    if (reportType == ReportType::Default) {
        fp = fopen(reportFileName.c_str(), "w");
    } else if (reportType == ReportType::EM) {
        fp = fopen(reportFileName_em.c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        fp = fopen(reportFileName_em_reclassify.c_str(), "w");
    }
    fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\tavg_score\trank\ttaxID\tname\n");
    writeReport(fp, cladeCounts, taxon2avgScore, numOfQuery);
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

void Reporter::writeReportFile(
    int numOfQuery, 
    unordered_map<TaxID, TaxonCounts> cladeCounts,
    ReportType reportType,
    string kronaFileName) 
{
    // std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    // unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
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

void Reporter::writeReport(
    FILE *FP, 
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts, 
    const unordered_map<TaxID, double> &taxon2avgScore,
    unsigned long totalReads, 
    TaxID taxID, 
    int depth) 
{
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        writeReport(FP, cladeCounts, taxon2avgScore, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);

        auto scoreIt = taxon2avgScore.find(taxID);
        if (scoreIt != taxon2avgScore.end()) {
            fprintf(FP, "%.4f\t%i\t%i\t%.4f\t%s\t%i\t%s%s\n",
                100 * cladeCount / double(totalReads),
                cladeCount, 
                taxCount,
                scoreIt->second,
                taxonomy->getString(taxon->rankIdx), 
                taxonomy->getOriginalTaxID(taxID), 
                std::string(2 * depth, ' ').c_str(), 
                taxonomy->getString(taxon->nameIdx));
        } else {
            fprintf(FP, "%.4f\t%i\t%i\t-\t%s\t%i\t%s%s\n",
                100 * cladeCount / double(totalReads),
                cladeCount, 
                taxCount,
                taxonomy->getString(taxon->rankIdx), 
                taxonomy->getOriginalTaxID(taxID), 
                std::string(2 * depth, ' ').c_str(), 
                taxonomy->getString(taxon->nameIdx));
        }
        
        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); });
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                writeReport(FP, cladeCounts, taxon2avgScore, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

void Reporter::writeReport(
    FILE *FP, 
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
    unsigned long totalReads,
    TaxID taxID,
    int depth)
{
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
    getReadsByCladeMembership(cladeId, readClassificationFileName, readIdxs, taxonomy, true);
}

void Reporter::getReadsNotClassifiedToClade(TaxID cladeId,
                                            const string &readClassificationFileName,
                                            vector<size_t> &readIdxs) {
    getReadsByCladeMembership(cladeId, readClassificationFileName, readIdxs, taxonomy, false);
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
        delete kseq;
        return;
    }

    if (readIdxs.empty()) {
        fclose(outFile);
        delete kseq;
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
    fclose(outFile);
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
    ReportType reportType,
    string kronaFileName)
{
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
    writeReportFile(numOfQuery, cladeCounts, reportType, kronaFileName);
}

void Reporter::writeReportFile(
    int numOfQuery,
    unordered_map<TaxID, unsigned int> &taxCnt,
    unordered_map<TaxID, CovMetric> &species2covMetrics,
    const unordered_map<TaxID, double> &taxon2avgScore,
    ReportType reportType,
    string kronaFileName)
{
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCnt, parentToChildren);
    writeReportFile(numOfQuery, cladeCounts, species2covMetrics, taxon2avgScore, reportType, kronaFileName);
}

void Reporter::writeReportFile(
    int numOfQuery,
    unordered_map<TaxID, TaxonCounts> cladeCounts,
    unordered_map<TaxID, CovMetric> &species2covMetrics,
    const unordered_map<TaxID, double> &taxon2avgScore,
    ReportType reportType,
    string kronaFileName)
{
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    rollUpCoverageMetrics(parentToChildren, cladeCounts, species2covMetrics, 1);
    FILE *fp = nullptr;
    
    if (reportType == ReportType::Default) {
        fp = fopen(reportFileName.c_str(), "w");
    } else if (reportType == ReportType::EM) {
        fp = fopen(reportFileName_em.c_str(), "w");
    } else if (reportType == ReportType::EM_RECLASSIFY) {
        fp = fopen(reportFileName_em_reclassify.c_str(), "w");
    }
    
    fprintf(fp, "#clade_proportion\tclade_count\ttaxon_count\trank\ttaxID\tevenness\tcoverage\tmacro_coverage\tadjusted_evenness\tunified_score\tavg_score\tname\n");
    writeReport(fp, cladeCounts, species2covMetrics, taxon2avgScore, numOfQuery);
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
    const std::unordered_map<TaxID, CovMetric> &species2covMetrics,
    const unordered_map<TaxID, double> &taxon2avgScore,
    unsigned long totalReads,
    TaxID taxID,
    int depth)
{
    std::unordered_map<TaxID, TaxonCounts>::const_iterator it = cladeCounts.find(taxID);
    unsigned int cladeCount = it == cladeCounts.end() ? 0 : it->second.cladeCount;
    unsigned int taxCount = it == cladeCounts.end() ? 0 : it->second.taxCount;
    
    if (taxID == 0) {
        if (cladeCount > 0) {
            fprintf(FP, "%.4f\t%i\t%i\tno rank\t0\t-\t-\t-\t-\t-\t-\tunclassified\n",
                    100 * cladeCount / double(totalReads),
                    cladeCount, taxCount);
        }
        writeReport(FP, cladeCounts, species2covMetrics, taxon2avgScore, totalReads, 1);
    } else {
        if (cladeCount == 0) {
            return;
        }
        const TaxonNode *taxon = taxonomy->taxonNode(taxID);
        
        fprintf(FP, "%.4f\t%i\t%i\t%s\t%i",
                100 * cladeCount / double(totalReads), cladeCount, taxCount,
                taxonomy->getString(taxon->rankIdx), taxonomy->getOriginalTaxID(taxID));

        auto evIt = species2covMetrics.find(taxID);
        if (evIt != species2covMetrics.end()) {
            writeMetricOrDash(FP, evIt->second.evenness);
            writeMetricOrDash(FP, evIt->second.coverage);
            writeMetricOrDash(FP, evIt->second.macroCoverage);
            writeMetricOrDash(FP, evIt->second.adjustedEvenness);
            writeMetricOrDash(FP, evIt->second.unifiedScore);
        } else {
            fprintf(FP, "\t-\t-\t-\t-\t-");
        }

        auto scoreIt = taxon2avgScore.find(taxID);
        if (scoreIt != taxon2avgScore.end()) {
            writeMetricOrDash(FP, scoreIt->second);
        } else {
            fprintf(FP, "\t-");
        }
        fprintf(FP, "\t%s%s\n", std::string(2 * depth, ' ').c_str(), taxonomy->getString(taxon->nameIdx));

        std::vector<TaxID> children = it->second.children;
        SORT_SERIAL(children.begin(), children.end(), [&](int a, int b) { 
            return cladeCountVal(cladeCounts, a) > cladeCountVal(cladeCounts, b); 
        });
        
        for (size_t i = 0; i < children.size(); ++i) {
            TaxID childTaxId = children[i];
            if (cladeCounts.count(childTaxId)) {
                writeReport(FP, cladeCounts, species2covMetrics, taxon2avgScore, totalReads, childTaxId, depth + 1);
            } else {
                break;
            }
        }
    }
}

void Reporter::filterClassificationFile(
    const std::string& inputFilePath,
    const std::string& outputFilePath,
    const std::unordered_map<TaxID, double> &taxon2avgScore,
    double cutoff)
{
    static const std::unordered_map<TaxID, CovMetric> noCoverageMetrics;
    filterClassificationFile(inputFilePath, outputFilePath, noCoverageMetrics, taxon2avgScore, cutoff);
}


void Reporter::filterClassificationFile(
    const std::string& inputFilePath,
    const std::string& outputFilePath,
    const std::unordered_map<TaxID, CovMetric> &sp2covMetric,
    const std::unordered_map<TaxID, double> &taxon2avgScore,
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
    bool hasLineageColumn = false;
    while (std::getline(inFile, line)) {
        // 1. Pass the header or empty lines directly to the new file
        if (line.empty()) {
            outFile << line << "\n";
            continue;
        }
        if (line[0] == '#') {
            std::vector<std::string> headerColumns = TaxonomyWrapper::splitByDelimiter(line, "\t", 16);
            for (const std::string &column : headerColumns) {
                if (column == "lineage") {
                    hasLineageColumn = true;
                    break;
                }
            }
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
            TaxID externalTaxID = 0;
            if (!parseTaxId(columns[2], externalTaxID)) {
                outFile << line << "\n";
                continue;
            }

            auto taxIdIt = ex2inTaxId.find(externalTaxID);
            if (taxIdIt == ex2inTaxId.end()) {
                outFile << line << "\n";
                continue;
            }

            TaxID taxID = taxIdIt->second;
            TaxID filterTaxID = taxID;

            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxID, "species");
            if (speciesTaxID != 0) {
                filterTaxID = speciesTaxID;
            }

            auto evIt = sp2covMetric.find(filterTaxID);
            auto avgIt = taxon2avgScore.find(filterTaxID);
            
            // 4. If the species is in the map AND its score is below the cutoff, rewrite it
            if ((evIt != sp2covMetric.end() && evIt->second.adjustedEvenness < 0.6) ||
                (avgIt != taxon2avgScore.end() && avgIt->second < cutoff)) {
                // Reconstruct the line exactly as your original "unclassified" else-block did
                outFile << "0\t"          // is_classified
                        << columns[1] << "\t" // name
                        << "0\t"          // taxID (forced to 0)
                        << columns[3] << "\t" // query_length
                        << columns[4] << "\t" // idScore
                        << "-\t"          // eValue
                        << "-";           // rank

                if (hasLineageColumn) {
                    outFile << "\t-\t";   // empty lineage
                } else {
                    outFile << "\t";
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

void Reporter::rollUpCoverageMetrics(
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren,
    const std::unordered_map<TaxID, TaxonCounts>& cladeCounts,
    std::unordered_map<TaxID, CovMetric>& allMetrics, // Starts with species, gets filled with all ranks
    TaxID currentTaxID) 
{
// NEW BASE CASE: If this node already has a metric calculated (i.e., it's a Species),
    // it is the base of our roll-up. Stop recursing down to subspecies and return immediately.
    if (allMetrics.find(currentTaxID) != allMetrics.end()) {
        return;
    }

    // Secondary Base Case: True leaf node (no children at all) but no pre-calculated metric
    auto childrenIt = parentToChildren.find(currentTaxID);
    if (childrenIt == parentToChildren.end() || childrenIt->second.empty()) {
        return; 
    }

    double weightedEvenness = 0.0;
    double weightedCoverage = 0.0;
    double weightedAdjEvenness = 0.0;
    double weightedUnifiedScore = 0.0;
    double weightedMacroCoverage = 0.0;
    uint64_t totalChildReads = 0;

    // 1. Recursively process all children FIRST (Post-order traversal)
    for (TaxID childID : childrenIt->second) {
        rollUpCoverageMetrics(parentToChildren, cladeCounts, allMetrics, childID);
        
        // 2. Gather data for the weighted average
        auto countIt = cladeCounts.find(childID);
        auto metricIt = allMetrics.find(childID);

        if (countIt != cladeCounts.end() && metricIt != allMetrics.end()) {
            uint64_t reads = countIt->second.cladeCount;
            
            // Only include children that actually have metrics and reads
            if (reads > 0) {
                weightedEvenness += metricIt->second.evenness * reads;
                weightedCoverage += metricIt->second.coverage * reads;
                weightedAdjEvenness += metricIt->second.adjustedEvenness * reads;
                weightedUnifiedScore += metricIt->second.unifiedScore * reads;
                weightedMacroCoverage += metricIt->second.macroCoverage * reads; // Assuming macro coverage is the same as coverage for weighting

                totalChildReads += reads;
            }
        }
    }

    // 3. Calculate the final weighted metrics for THIS parent node
    if (totalChildReads > 0) {
        CovMetric parentMetric;
        parentMetric.evenness = weightedEvenness / static_cast<double>(totalChildReads);
        parentMetric.coverage = weightedCoverage / static_cast<double>(totalChildReads);
        parentMetric.adjustedEvenness = weightedAdjEvenness / static_cast<double>(totalChildReads);
        parentMetric.unifiedScore = weightedUnifiedScore / static_cast<double>(totalChildReads);
        parentMetric.macroCoverage = weightedMacroCoverage / static_cast<double>(totalChildReads);
        
        // Save it to the map
        allMetrics[currentTaxID] = parentMetric;
    } 
}
