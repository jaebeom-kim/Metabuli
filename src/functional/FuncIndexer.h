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
#include <unordered_set>

struct QueryCDS {
    QueryCDS(uint32_t codingRegionId, int cdsLength, const string & assigend_uniref_id, uint32_t assigend_uniref_idx, float score)
        : codingRegionId(codingRegionId), cdsLength(cdsLength), assigend_uniref_id(assigend_uniref_id), assigend_uniref_idx(assigend_uniref_idx), score(score) {}
    QueryCDS(uint32_t codingRegionId, int cdsLength)
        : codingRegionId(codingRegionId), cdsLength(cdsLength) {}
    QueryCDS() = default;
    uint32_t codingRegionId;
    int cdsLength;
    string assigend_uniref_id;
    uint32_t assigend_uniref_idx;
    float score;
};

struct QueryCodingRegionInfo {
    QueryCodingRegionInfo(uint32_t codingRegionId, int cdsLength)
        : codingRegionId(codingRegionId), cdsLength(cdsLength) {}
    QueryCodingRegionInfo() = default;
    uint32_t codingRegionId;
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

class FuncIndexer : public IndexCreator {
protected:
    // Inputs
    string protDBFileName;
    string protDbSplitFileName;
    string protIdMapFileName;
    string unirefIdx2taxIdFileName;
    string ncbi2gtdbFileName;

    // Outputs
    string protIdx2taxIdFileName;
    string protIdx2unirefIdFileName;
    string protIdx2taxIdFileName_debug;
    string protIdx2unirefIdFileName_debug;
    // string diffIdxFileName;
    // string infoFileName;
    // string deltaIdxFileName;

    KmerMatcher *kmerMatcher;
    uint32_t lastProtIdx;
    std::unordered_map<TaxID, TaxID> taxId2speciesId;
    std::unordered_set<uint32_t> codingRegionIdSet;
    std::unordered_map<string, vector<CDSinfo>> cdsInfoMap;
    std::unordered_map<uint32_t, string> uniRefIdMap;
    std::unordered_map<uint32_t, string> protIdMap;
    std::unordered_map<uint32_t, QueryCodingRegionInfo> queryCodingRegionMap;
    std::unordered_map<uint32_t, int> cds2lengthMap;
    std::unordered_map<uint32_t, std::vector<ProtScore>> query2targetProtScMap;
    std::unordered_map<uint32_t, TaxID> regionId2taxId;
    std::unordered_map<uint32_t, uint32_t> regionId2unirefId;
    std::unordered_map<uint32_t, int> unirefIdx2taxId;
    std::unordered_map<TaxID, TaxID> ncbi2gtdb;
    std::unordered_map<TaxID, std::unordered_set<TaxID>> spTaxId2UniRefTaxIds;
    std::unordered_set<TaxID> speciesTaxIds;
    std::unordered_set<TaxID> unirefTaxIds;



    void loadProtIdMap();

    void makeSpTaxId2UniRefTaxIds();

    size_t fillTargetKmerBuffer(Buffer<ExtractedMetamer> &kmerBuffer,
                                bool * tempChecker,
                                size_t &processedSplitCnt,
                                uint32_t & nonCdsIdx);
    
    size_t fillTargetKmerBufferUsingProdigal(Buffer<ExtractedMetamer> &kmerBuffer,
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

    void writeTargetFilesAndSplits(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void writeTargetFilesAndSplits2(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void writeTargetFiles(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void writeTargetFiles2(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void mergeDeltaIndexFiles();

    size_t getSmallestMetamer(const Metamer * lookingMetamers, size_t fileCnt); 

    void fillDeltaIndexing(const Metamer & previousMetamer,
                           const Metamer & currentMetamer,
                           FILE* handleKmerTable,
                           uint16_t * deltaIndexBuffer,
                           size_t bufferSize,
                           size_t & localBufIdx);

    void fillDeltaIndexing(const Metamer & previousMetamer,
                           const Metamer & currentMetamer,
                           FILE* handleKmerTable,
                           uint16_t * deltaIndexBuffer,
                           size_t bufferSize,
                           size_t & localBufIdx,
                           size_t & totalBufferIdx);

    void writeTargetFilesInText(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqeKmerIdx, size_t & uniqKmerCnt);

    void testDeltaIndexing(const Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqeKmerIdx, size_t uniqKmerCnt);

    uint32_t chooseBestProtein(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    ProtScore scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    static bool sortTargetMetamerF(const TargetMetamerF &a, const TargetMetamerF &b);

    static bool sortExtractedMetamer(const ExtractedMetamer &a, const ExtractedMetamer &b);

    void loadCdsInfo(const string & cdsInfoFileList);

    void generateTaxId2SpeciesIdMap();

    
public:
    FuncIndexer(const LocalParameters &par);
    ~FuncIndexer();
    void createIndex();

};
#endif //METABULI_PROTEINDBINDEXER_H