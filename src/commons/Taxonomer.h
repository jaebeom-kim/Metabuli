#ifndef METABULI_TAXONOMER_H
#define METABULI_TAXONOMER_H

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <type_traits>

#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "Match.h"
#include "common.h"
#include "BitManipulateMacros.h"
#include "MetamerPattern.h"
#include "EvalueComputation.h"

using namespace std;

struct TaxonScore {
    TaxID taxId;
    MatchScore score;
    int hammingDist;
    bool LCA;
    std::pair<size_t, size_t> matchPathRange;
    TaxonScore(TaxID taxId, MatchScore score, int hammingDist, bool LCA) :
            taxId(taxId), score(score), hammingDist(hammingDist), LCA(LCA) {}
    TaxonScore() : taxId(0), score(), hammingDist(0), LCA(false) {}
};

template <typename MatchType>
class Taxonomer {
private:
    const LocalParameters & par;
    TaxonomyWrapper * taxonomy;
    const MetamerPattern *metamerPattern = nullptr;
    SubstitutionMatrix * substitutionMatrix = nullptr;
    EvalueComputation *evaluer = nullptr;

    // spaced k-mer
    int bitPerCodon;
    int bitPerAA;
    int spaceNum;
    int kmerLen;
    int windowSize;
    uint32_t windowMask;

    // Parameters from user
    int accessionLevel;
    int eukaryotaTaxId;
    std::vector<TaxID> priorityTaxa;

    // Internal
    int denominator;
    int maxCodonShift;
    int dnaShift;
    // int smerLength;
    double logMaxEValue;
    bool useEvalueFilter = false;

    // vector<const Match *> speciesMatches;

    // chooseBestTaxon
    unordered_map<TaxID, unsigned int> taxCnt;

    // getBestSpeciesMatches
    vector<MatchPath<MatchType>> matchPaths;
    vector<MatchPath<MatchType>> combinedMatchPaths;
    vector<TaxID> maxSpecies;

    // getMatchPaths
    vector<bool> connectedToNext;
    vector<MatchPath<MatchType>> localMatchPaths;

    // lowerRankClassification
    unordered_map<TaxID, TaxonCounts> cladeCnt;

    // filterRedundantMatches
    const MatchType **bestMatchForQuotient;
    TaxID *bestMatchTaxIdForQuotient;
    uint8_t *minHammingForQuotient;
    size_t arraySize_filterRedundantMatches;


    // Output
    unordered_map<TaxID, unsigned int> taxCounts;

    void ensureArraySize(size_t newSize);

    void printSpeciesMatches (
       const MatchType *matchList,
       const std::pair<size_t, size_t> & bestSpeciesRange
    );

    
    TaxonScore getBestSpeciesMatches(
        std::pair<size_t, size_t> & bestSpeciesRange,
        const MatchType *matchList,
        size_t end,
        size_t offset,
        Query & query);

    MatchPath<MatchType> makeMatchPath(
        const MatchType * match
    );

    // void makeMatchPath(
    //     const MatchType * match,
    //     const MatchType * matchList,
    //     size_t matchNum,
    //     vector<MatchPath<MatchType>> & matchPaths,
    //     TaxID speciesId);

    void getMatchPaths_lookbackDP(
        const MatchType * matchList,
        size_t matchNum,
        vector<MatchPath<MatchType>> & matchPaths,
        TaxID speciesId);

    void getSpacedMatchPaths_lookbackDP(
        const MatchType * matchList,
        size_t matchNum,
        vector<MatchPath<MatchType>> & matchPaths,
        TaxID speciesId);

    void getMatchPaths(
        const MatchType * matchList,
        size_t matchNum,
        vector<MatchPath<MatchType>> & matchPaths,
        TaxID speciesId);
    
    void getSpacedMatchPaths(
        const MatchType * matchList,
        size_t matchNum,
        vector<MatchPath<MatchType>> & matchPaths,
        TaxID speciesId); 

    void makeSpacedMatchPath(
        const MatchType * match,
        size_t index
    );

    MatchScore combineMatchPaths(
        vector<MatchPath<MatchType>> & matchPaths,
        size_t matchPathStart,
        vector<MatchPath<MatchType>> & combMatchPaths,
        size_t combMatchPathStart,
        int queryLength);
        
    bool isMatchPathOverlapped(const MatchPath<MatchType> & matchPath1, const MatchPath<MatchType> & matchPath2);

    bool trimMatchPath(MatchPath<MatchType> & path1, const MatchPath<MatchType> & path2, int overlapLength);
    bool trimSpacedMatchPath(MatchPath<MatchType> & path1, const MatchPath<MatchType> & path2, int overlapLength);
    void sortMatchPath(std::vector<MatchPath<MatchType>> & matchPaths, size_t i);

   
public:

    unordered_map<TaxID, vector<uint8_t>> sp2coverage;
    unordered_map<TaxID, std::bitset<65536>> sp2coveredBins;
    unordered_map<TaxID, double> sp2scoreSum;

    Taxonomer(
        const LocalParameters & par, 
        TaxonomyWrapper * taxonomy, 
        const MetamerPattern *metamerPattern);
        
    ~Taxonomer();

    void chooseBestTaxon(uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         const MatchType *matchList,
                         vector<Query> & queryList);      

    void filterRedundantMatches(const MatchType *matchList,
                                const std::pair<size_t, size_t> & bestSpeciesRange,
                                unordered_map<TaxID, unsigned int> & taxCnt,
                                int queryLength);

    TaxID lowerRankClassification(const unordered_map<TaxID, unsigned int> & taxCnt, TaxID speciesID, int queryLength);

    void getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> & taxCnt,
                               unordered_map<TaxID, TaxonCounts> & cladeCnt,
                               TaxID spciesID);

    TaxID BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root, unsigned int maxCnt);

    // Getters
    unordered_map<TaxID, unsigned int> & getTaxCounts() { return taxCounts; }

};


#endif //METABULI_TAXONOMER_H