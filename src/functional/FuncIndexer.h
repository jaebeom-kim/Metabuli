#ifndef METABULI_FUNCINDEXER_H
#define METABULI_FUNCINDEXER_H

#include "IndexCreator.h"
#include "Kmer.h"
#include "KmerMatcher.h"
#include "LocalParameters.h"
#include "Match.h"
#include "common.h"
#include <cstdint>
#include <unordered_map>

struct QueryCDS {
    QueryCDS(const std::string &cdsId, int cdsIdx, int cdsLength, const std::string &proteinId, uint32_t protIdx, float score)
        : cdsId(cdsId), cdsIdx(cdsIdx), cdsLength(cdsLength), proteinId(proteinId), protIdx(protIdx), score(score) {}
    QueryCDS(const std::string &cdsId, int cdsIdx, int cdsLength)
        : cdsId(cdsId), cdsIdx(cdsIdx), cdsLength(cdsLength) {}
    QueryCDS() = default;
    std::string cdsId;
    int cdsIdx;
    int cdsLength;
    std::string proteinId;
    uint32_t protIdx;
    float score;
    
};

struct ProtScore {
    ProtScore(uint32_t protIdx, float score) : protIdx(protIdx), score(score) {}
    ProtScore(uint32_t protIdx) : protIdx(protIdx), score(0) {}
    ProtScore() = default;
    uint32_t protIdx;
    float score;
};

// struct cds2protScore {
//     cds2protScore(const std::string &cdsId, uint32_t protIdx, float score)
//         : cdsId(cdsId), protIdx(protIdx), score(score) {}
//     cds2protScore() = default;
//     std::string cdsId;
//     uint32_t protIdx;
//     float score;
// };

class FuncIndexer : public IndexCreator {
protected:
    const LocalParameters &par;
    string protDBFileName;
    string protDbSplitFileName;
    string protIdMapFileName;

    KmerMatcher *kmerMatcher;

    std::unordered_map<uint32_t, string> protIdMap;
    std::vector<QueryCDS> queryCdsList;
    std::unordered_map<uint32_t, int> cds2lengthMap;
    std::unordered_map<uint32_t, std::vector<ProtScore>> cds2protScoreMap;

    void loadProtIdMap();

    size_t fillTargetKmerBuffer(Buffer<TargetMetamerF> &kmerBuffer,
                                bool * tempChecker,
                                size_t &processedSplitCnt);
    
    int getProteinId(Buffer<TargetMetamerF> &kmerBuffer);

    void labelCdsWithProtId(Buffer<ProtMatch> & protMatch, unordered_map<uint32_t, uint32_t> & cdsId2protId);

    void labelKmerWithProtId(Buffer<TargetMetamerF> &kmerBuffer, const unordered_map<uint32_t, uint32_t> & cdsId2protId);

    uint32_t chooseBestProtein(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    ProtScore scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    bool sortTargetMetamerF(const TargetMetamerF &a, const TargetMetamerF &b);
public:
    FuncIndexer(const LocalParameters &par);
    ~FuncIndexer();
    void createIndex();

};
#endif //METABULI_PROTEINDBINDEXER_H