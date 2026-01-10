#ifndef METABULI_TAXONOMER_H
#define METABULI_TAXONOMER_H
#include "TaxonomyWrapper.h"
#include "LocalParameters.h"
#include "Match.h"
#include "common.h"
#include "BitManipulateMacros.h"
#include "MetamerPattern.h"
#include <unordered_map>
#include <unordered_set>
#include <algorithm>
#include <cstdint>

using namespace std;

struct TaxonScore {
    TaxID taxId;
    float score;
    int hammingDist;
    bool LCA;
    TaxonScore(TaxID taxId, float score, int hammingDist, bool LCA) :
            taxId(taxId), score(score), hammingDist(hammingDist), LCA(LCA) {}
    TaxonScore() : taxId(0), score(0.0f), hammingDist(0), LCA(false) {}
};

struct MatchPath {
    MatchPath(int start, int end, float score, int hammingDist, int depth, const Match2 * startMatch, const Match2 * endMatch) :
         start(start), end(end), score(score), hammingDist(hammingDist), depth(depth), startMatch(startMatch), endMatch(endMatch) {}
    MatchPath() : start(0), end(0), score(0.f), hammingDist(0), depth(0), startMatch(nullptr), endMatch(nullptr) {}


    MatchPath(const Match2 * startMatch, float score, int hammingDist) 
        : start(startMatch->qKmer.qInfo.pos),
          end(startMatch->qKmer.qInfo.pos + 23),
          score(score),
          hammingDist(hammingDist),
          depth(1),
          startMatch(startMatch),
          endMatch(startMatch) {}
    
    int start;                // query coordinate
    int end;                  // query coordinate
    float score;
    int hammingDist;
    int depth;
    const Match2 * startMatch;
    const Match2 * endMatch;

    void printMatchPath() {
        std::cout << start << " " << end << " " << score << " " << hammingDist << " " << depth << std::endl;
    }
};

class Taxonomer {
private:
    const LocalParameters & par;
    TaxonomyWrapper * taxonomy;
    const MetamerPattern *metamerPattern = nullptr;
    const SubstitutionMatrix * substitutionMatrix = nullptr;
    

    // spaced k-mer
    int unmaskedPos[9];
    int spaceNum;

    // Parameters from user
    int maxGap;
    int accessionLevel;
    int minSSMatch;
    size_t minConsCnt;
    size_t minConsCntEuk;
    int eukaryotaTaxId;
    float tieRatio;

    // Internal
    int denominator;
    int bitsPerCodon;
    int totalDnaBits;
    uint32_t lastCodonMask;
    int maxCodonShift;
    int dnaShift;
    int smerLength;
    int minSubSpeciesMatch;

    // vector<const Match *> speciesMatches;

    // chooseBestTaxon
    unordered_map<TaxID, unsigned int> taxCnt;

    // getBestSpeciesMatches
    vector<MatchPath> matchPaths;
    vector<MatchPath> combinedMatchPaths;
    vector<TaxID> maxSpecies;

    // getMatchPaths
    vector<bool> connectedToNext;
    vector<MatchPath> localMatchPaths;

    // lowerRankClassification
    unordered_map<TaxID, TaxonCounts> cladeCnt;

    // filterRedundantMatches
    const Match2 **bestMatchForQuotient;
    TaxID *bestMatchTaxIdForQuotient;
    uint8_t *minHammingForQuotient;
    size_t arraySize_filterRedundantMatches;


    // Output
    unordered_map<TaxID, unsigned int> taxCounts;

    void ensureArraySize(size_t newSize);

    void printSpeciesMatches (
       const Match2 *matchList,
       const std::pair<size_t, size_t> & bestSpeciesRange
    );

    TaxonScore getBestSpeciesMatches(
        std::pair<size_t, size_t> & bestSpeciesRange,
        const Match2 *matchList,
        size_t end,
        size_t offset,
        Query & query);

    void getMatchPaths(
        const Match2 * matchList,
        size_t start,
        size_t end,
        vector<MatchPath> & matchPaths,
        TaxID speciesId); 
    
    // void getMatchPaths2(
    //     const Match2 * matchList,
    //     size_t start,
    //     size_t end,
    //     vector<MatchPath> & matchPaths,
    //     TaxID speciesId);

    float combineMatchPaths(
        vector<MatchPath> & matchPaths,
        size_t matchPathStart,
        vector<MatchPath> & combMatchPaths,
        size_t combMatchPathStart,
        int readLength);
        
    bool isMatchPathOverlapped(const MatchPath & matchPath1, const MatchPath & matchPath2);
    void trimMatchPath(MatchPath & path1, const MatchPath & path2, int overlapLength);

public:
    Taxonomer(
        const LocalParameters & par, 
        TaxonomyWrapper * taxonomy, 
        const MetamerPattern *metamerPattern);

    ~Taxonomer();

    // void assignTaxonomy(const Match *matchList,
    //                     size_t numOfMatches,
    //                     std::vector<Query> & queryList,
    //                     const LocalParameters &par);

    void chooseBestTaxon(uint32_t currentQuery,
                         size_t offset,
                         size_t end,
                         const Match2 *matchList,
                         vector<Query> & queryList,
                         const LocalParameters &par);      

    void filterRedundantMatches(const Match2 *matchList,
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

    bool compareMatchPaths(const MatchPath& a, const MatchPath& b) const {
        if (a.score != b.score)
            return a.score < b.score;
        if (a.hammingDist != b.hammingDist)
            return a.hammingDist < b.hammingDist;
        return a.start < b.start;
    }
};


#endif //METABULI_TAXONOMER_H