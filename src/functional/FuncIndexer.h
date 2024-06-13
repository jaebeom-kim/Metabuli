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
    QueryCDS(uint32_t cds_protein_id, int cdsIdx, int cdsLength, const string & assigend_uniref_id, uint32_t assigend_uniref_idx, float score)
        : cds_protein_id(cds_protein_id), cdsIdx(cdsIdx), cdsLength(cdsLength), assigend_uniref_id(assigend_uniref_id), assigend_uniref_idx(assigend_uniref_idx), score(score) {}
    QueryCDS(uint32_t cds_protein_id, int cdsIdx, int cdsLength)
        : cds_protein_id(cds_protein_id), cdsIdx(cdsIdx), cdsLength(cdsLength) {}
    QueryCDS() = default;
    uint32_t cds_protein_id;
    int cdsIdx;
    int cdsLength;
    string assigend_uniref_id;
    uint32_t assigend_uniref_idx;
    float score;
};

struct ProtScore {
    ProtScore(uint32_t targetProtIdx, float score) : targetProtIdx(targetProtIdx), score(score) {}
    ProtScore(uint32_t targetProtIdx) : targetProtIdx(targetProtIdx), score(0) {}
    ProtScore() = default;
    uint32_t targetProtIdx;
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

    uint32_t lastProtIdx;
    std::unordered_map<string, vector<CDSinfo>> cdsInfoMap;
    std::unordered_map<uint32_t, string> uniRefIdMap;
    std::unordered_map<uint32_t, string> protIdMap;
    std::vector<QueryCDS> queryCdsList;
    std::unordered_map<uint32_t, int> cds2lengthMap;
    std::unordered_map<uint32_t, std::vector<ProtScore>> query2targetProtScMap;
    std::unordered_map<uint32_t, TaxID> protIdx2taxId;
    std::unordered_map<uint32_t, uint32_t> protIdx2unirefId;


    void loadProtIdMap();

    size_t fillTargetKmerBuffer(Buffer<ExtractedMetamer> &kmerBuffer,
                                bool * tempChecker,
                                size_t &processedSplitCnt,
                                uint32_t & nonCdsIdx);
    
    int getProteinId(Buffer<ExtractedMetamer> &kmerBuffer);

    int getProteinId(Buffer<TargetMetamerF> &kmerBuffer);

    void mapQueryProt2TargetProt(Buffer<ProtMatch> & protMatch);

    void labelExtractedMetamerWithUniRefId(Buffer<ExtractedMetamer> &kmerBuffer);

    void labelKmerWithProtId(Buffer<TargetMetamerF> &kmerBuffer, const unordered_map<uint32_t, uint32_t> & cdsId2protId);

    void reduceRedundancy(Buffer<TargetMetamerF> &kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void reduceRedundancy(Buffer<ExtractedMetamer> &kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    uint32_t chooseBestProtein(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    ProtScore scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    static bool sortTargetMetamerF(const TargetMetamerF &a, const TargetMetamerF &b);

    static bool sortExtractedMetamer(const ExtractedMetamer &a, const ExtractedMetamer &b);

    void loadCdsInfo(const string & cdsInfoFileList);
public:
    FuncIndexer(const LocalParameters &par);
    ~FuncIndexer();
    void createIndex();

};
#endif //METABULI_PROTEINDBINDEXER_H