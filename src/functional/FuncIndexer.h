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
#include "KaijuWrapper.h"

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
    string unirefId2IdxFileName;
    string unirefIdx2taxIdFileName;
    string ncbi2gtdbFileName;

    // Inputs (from Kaiju)
    string kaijuFmiFileName;
    string kaijuUniRefIdsFileName;
    string kaijuUniRefTaxIdsFileName;

    // Outputs
    string regionId2taxIdFileName;
    string regionId2unirefIdFileName;
    string regionId2taxIdFileName_debug;
    string regionId2unirefIdFileName_debug;


    
    uint32_t lastProtIdx;
    std::unordered_map<TaxID, TaxID> taxId2speciesId;
    std::unordered_map<string, vector<CDSinfo>> cdsInfoMap;

    // Required to annoate with UniRef
    KaijuWrapper *kaiju;
    std::vector<string> uniRefIds;
    std::vector<TaxID> uniRefTaxIds;
    std::unordered_set<TaxID> speciesTaxIds;
    std::unordered_map<TaxID, TaxID> ncbi2gtdb;

    // Results
    std::unordered_map<uint32_t, TaxID> regionId2taxId;
    std::unordered_map<uint32_t, uint32_t> regionId2unirefId;

    // Debug
    std::unordered_map<uint32_t, string> uniRefIdx2Id;

    std::unordered_map<uint32_t, QueryCodingRegionInfo> queryCodingRegionMap;

    // void loadUniRefId2Idx();
    // void loadUniRefIdx2TaxId();

    void makeSpTaxId2UniRefTaxIds();

    uint32_t getUniRefIdx(const string & cds, TaxID taxIds);

    size_t fillTargetKmerBuffer(Buffer<ExtractedMetamer> &kmerBuffer,
                                bool * tempChecker,
                                size_t &processedSplitCnt,
                                uint32_t & nonCdsIdx);

    // size_t fillTargetKmerBufferUsingProdigal(Buffer<ExtractedMetamer> &kmerBuffer,
    //                                          bool * tempChecker,
    //                                          size_t &processedSplitCnt,
    //                                          uint32_t & nonCdsIdx);

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

    // uint32_t chooseBestProtein(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    // ProtScore scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch);

    // static bool sortTargetMetamerF(const TargetMetamerF &a, const TargetMetamerF &b);

    static bool sortExtractedMetamer(const ExtractedMetamer &a, const ExtractedMetamer &b);

    void loadCdsInfo(const string & cdsInfoFileList);

    void generateTaxId2SpeciesIdMap();

    int selectReadingFrame(unordered_set<uint64_t> & kmerSet, const string & seq);

    
public:
    FuncIndexer(const LocalParameters &par);
    ~FuncIndexer();
    void createIndex();
    // std::unordered_map<string, uint32_t> & getUniRefId2Idx() { return uniRefId2Idx; }

};
#endif //METABULI_PROTEINDBINDEXER_H