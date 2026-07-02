#ifndef METABULI_TAXONOMER_H
#define METABULI_TAXONOMER_H

#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>
#include <type_traits>

#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "CandidateDBReader.h"
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

struct SpeciesScoreCandidate {
    TaxID taxId;
    MatchScore score;
    std::pair<size_t, size_t> speciesRange;
    std::pair<size_t, size_t> matchPathRange;

    SpeciesScoreCandidate(TaxID taxId,
                          MatchScore score,
                          std::pair<size_t, size_t> speciesRange,
                          std::pair<size_t, size_t> matchPathRange)
        : taxId(taxId),
          score(score),
          speciesRange(speciesRange),
          matchPathRange(matchPathRange) {}
};

template <typename MatchType>
class Taxonomer {
private:
    const LocalParameters & par;
    TaxonomyWrapper * taxonomy;
    const MetamerPattern *metamerPattern = nullptr;
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
    std::vector<uint8_t> ownedPriorityTaxonLookup;
    const std::vector<uint8_t> *priorityTaxonLookup = &ownedPriorityTaxonLookup;

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
    vector<SpeciesScoreCandidate> sp2score;
    vector<size_t> tiedIndices;
    vector<TaxID> tempPrioritySpecies;
    vector<size_t> tempPriorityIndices;

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
    vector<size_t> touchedQuotients;


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

    void addTerminalMatchPath(
        size_t pathIdx,
        const MatchType * matchList,
        vector<MatchPath<MatchType>> & filteredMatchPaths);

   
public:

    unordered_map<TaxID, vector<uint8_t>> sp2coverage;
    unordered_map<TaxID, uint64_t> sp2totalReadLength;

    Taxonomer(
        const LocalParameters & par, 
        TaxonomyWrapper * taxonomy, 
        const MetamerPattern *metamerPattern,
        const std::vector<uint8_t> *priorityTaxonLookup = nullptr);
        
    ~Taxonomer();

    static std::vector<uint8_t> makePriorityTaxonLookup(
        const LocalParameters &par,
        TaxonomyWrapper *taxonomy);

    std::unordered_map<TaxID, double> sp2scoreSum;

    void chooseBestTaxon(uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         const MatchType *matchList,
                         vector<Query> & queryList);      

    void chooseBestTaxonFromCandidates(
        const CandidateDBEntry &candidateEntry,
        Query &query);

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
