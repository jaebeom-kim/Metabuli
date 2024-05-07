#ifndef METABULI_PROTEINDBSEARCHER_H
#define METABULI_PROTEINDBSEARCHER_H

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

class ProteinDbSearcher{
private:
    // Parameters
    const LocalParameters &par;
    int threadNum;
    size_t bufferSize;
    std::string dbName;
    
    // Inputs
    std::string dbDir;
    std::string proteinFileName;
    std::string proteinIndexSplitFileName;

    // Internal
    uint64_t mask_getId;
    uint64_t mask_getAA;
    Buffer<uint64_t> * aaKmerBuffer;
    std::vector<SequenceBlock> sequenceBlocks;
    std::unordered_map<string, int> proteinId2Index;
    int flushCnt;

    // Outputs
    std::string versionFileName;
    std::string paramterFileName;

    size_t numOfFlush=0;

    void splitFasta(const std::string & fastaFile, std::vector<SequenceBlock> & blocks);
    size_t fillAAKmerBuffer(bool *checker, size_t & processedSplitCnt);
    void reduceRedundantAAKmers(size_t * uniqeKmerIdx, size_t & uniqueKmerCnt);
    void writeProteinIndexFile(const size_t * uniqeKmerIdx, size_t uniqKmerCnt, bool writeSplit=false);
    void mergeProteinIndexFiles();
    size_t getSmallestKmer(const uint64_t lookingKmers[], size_t fileCnt); 

public:
    ProteinDbSearcher(const LocalParameters &par);
    ~ProteinDbSearcher();
    void splitProteinFile();
    void index();

};
#endif //METABULI_PROTEINDBSEARCHER_H
