#ifndef METABULI_PROTEINDBINDEXER_H
#define METABULI_PROTEINDBINDEXER_H

#include <cstdint>
#include <cstdio>
#include <string>
#include <unordered_map>
#include "LocalParameters.h"
#include "common.h"
#include "kseq.h"
#include "SeqIterator.h"
#include "KmerBuffer.h"
#include "KSeqWrapper.h"
#include "FileUtil.h"

class ProteinDbIndexer{
private:
    // Parameters
    const LocalParameters &par;
    int threadNum;
    size_t bufferSize;
    std::string dbName;
    std::string dbDate;
    
    // Inputs
    std::string dbDir;
    std::string proteinFileName;
    std::string proteinIndexSplitFileName;

    // Internal
    uint64_t mask_getId;
    uint64_t mask_getAA;
    Buffer<uint64_t> * aaKmerBuffer;
    std::vector<SequenceBlock> sequenceBlocks;
    std::unordered_map<string, size_t> proteinId2Index;
    int flushCnt;

    // Outputs
    std::string versionFileName;
    std::string paramterFileName;

    size_t numOfFlush=0;

    void splitFasta(const std::string & fastaFile, std::vector<SequenceBlock> & blocks);
    size_t fillAAKmerBuffer(bool *checker, size_t & processedSplitCnt);
    void reduceRedundantAAKmers(size_t * uniqeKmerIdx, size_t & uniqueKmerCnt);
    void writeProteinIndexFile(const size_t * uniqeKmerIdx, size_t uniqKmerCnt, bool writeSplit=false);
    void writeProteinDeltaIndexFile(const size_t * uniqeKmerIdx, size_t uniqKmerCnt, bool writeSplit=false);
    void mergeProteinIndexFiles();
    void mergeProteinDeltaIndexFiles();
    size_t getSmallestKmer(const uint64_t lookingKmers[], size_t fileCnt); 
    void writePrtIdMap();

public:
    ProteinDbIndexer(const LocalParameters &par);
    ~ProteinDbIndexer();
    void splitProteinFile();
    void index();

};
#endif //METABULI_PROTEINDBINDEXER_H
