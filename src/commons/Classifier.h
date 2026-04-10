#ifndef METABULI_CLASSIFIER_H
#define METABULI_CLASSIFIER_H

#include "BitManipulateMacros.h"
#include "Mmap.h"
#include <fstream>
#include "Kmer.h"
#include "SeqIterator.h"
#include "printBinary.h"
#include "common.h"
#include "NcbiTaxonomy.h"
#include "Debug.h"
#include "IndexCreator.h"
#include <cstdio>
#include <ctime>
#include <vector>
#include <algorithm>
#include <regex>
#include "FastSort.h"
#include "KSeqWrapper.h"
#include "LocalParameters.h"
#include <set>
#include <cmath>
#include "Match.h"
#include <unordered_set>
#include "LocalUtil.h"
#include "QueryIndexer.h"
#include "KmerMatcher.h"
#include "KmerExtractor.h"
#include "Taxonomer.h"
#include "Reporter.h"
#include <cstdint>
#include "TaxonomyWrapper.h"

#define BufferSize 16'777'216 //16 * 1024 * 1024 // 16 M
using namespace std;




class Classifier {
protected:
    // Parameters
    LocalParameters & par;
    string dbDir;
    size_t matchPerKmer;
    int kmerFormat;

    // Agents
    GeneticCode * geneticCode = nullptr;
    MetamerPattern * metamerPattern = nullptr;
    QueryIndexer * queryIndexer = nullptr;
    KmerExtractor * kmerExtractor = nullptr;
    KmerMatcher * kmerMatcher = nullptr;
    Reporter * reporter = nullptr;
    TaxonomyWrapper * taxonomy = nullptr;

    unordered_map<TaxID, unsigned int> taxCounts;

    size_t mappingResListSize = 0;

    // EM algorithm
    unordered_map<TaxID, unsigned int> emTaxCounts;
    unordered_map<TaxID, unsigned int> reclassifyTaxCounts;
    unordered_map<TaxID, double> taxProbs;
    MappingRes * mappingResList = nullptr;
    std::vector<Classification> emResults;
    std::unordered_set<TaxID> topSpeciesSet;

    // Coverage
    unordered_map<TaxID, uint64_t> sp2genomeSize;
    unordered_map<TaxID, uint64_t> sp2totalReadLength;
    unordered_map<TaxID, vector<uint8_t>> sp2coverage_global;
    unordered_map<TaxID, double> sp2scoreSum_global;
    unordered_map<TaxID, CovMetric> sp2covMetric;
    std::vector<double> C_LOG2_C;

    CovMetric calCovMetrics(
        const std::vector<uint8_t>& bins, 
        double k_mers_per_read, 
        int readCnt, 
        uint64_t totalReadLength,
        uint64_t genomeSize = 65536);

    void countUniqueKmerPerSpecies(vector<uint32_t> & sp2uniqKmerCnt);

    void loadMappings(const string & mappingResFileName);

    void loadOriginalResults(const string & classificationFileName, size_t seqNum);

    void preciseModePreset(LocalParameters & par);

    void parseSp2GenomeSize() {
        std::string fileName = dbDir + "/species2genomeSize.tsv";
        ifstream infile(fileName);
        if (!infile.is_open()) {
            cerr << "Error: Could not open species to genome size file " << fileName << endl;
            return;
        }
        string line;
        while (getline(infile, line)) {
            if (line.empty()) continue;
            vector<string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 3);
            sp2genomeSize[static_cast<TaxID>(stoul(columns[0]))] = stoull(columns[2]);
        }
        infile.close();
    }

    unordered_map<TaxID, TaxonCounts> getCladeCounts() {
        return taxonomy->getCladeCounts(taxCounts, taxonomy->getParentToChildren());
    }

public:
    void classifyReads();
    void classifyReadsWithPos();

    uint64_t calculateBufferSize(
        uint64_t queryListSize,
        uint64_t matchPerKmer,
        int mode);

    template <typename MatchType>
    void assignTaxonomy(const MatchType *matchList,
                        size_t numOfMatches,
                        std::vector<Query> & queryList,
                        const LocalParameters &par);

    explicit Classifier(LocalParameters & par);

    virtual ~Classifier();

    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }

    void getTopSpecies(const std::vector<Query> & queries) {
        for (const auto & query : queries) {
            if (query.classification != 0) {
                topSpeciesSet.insert(query.topSpeciesId);
            }
        }
    }

    void em(size_t totalQueryCnt);

    void reclassify(
        const std::vector<std::pair<size_t, size_t>> & queryRanges,
        const MappingRes * mappingResList,
        const vector<double> & sp2lengthFactor,
        size_t totalQueryCnt);

};


#endif //METABULI_CLASSIFIER_H
