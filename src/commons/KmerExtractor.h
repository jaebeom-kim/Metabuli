#ifndef METABULI_KMEREXTRACTER_H
#define METABULI_KMEREXTRACTER_H
#include "SeqIterator.h"
#include "QueryIndexer.h"
#include "KSeqWrapper.h"
#include "common.h"
#include <unordered_map>
#include <atomic>
#include <cstdint>
#include "GeneticCode.h"
#include "MetamerPattern.h"

class KmerExtractor {
private:
    const LocalParameters &par;
    const GeneticCode * geneticCode;
    const MetamerPattern * metamerPattern;
    std::vector<std::unique_ptr<KmerScanner>> kmerScanners;
    
    // Parameters
    int spaceNum;
    int maskMode;
    float maskProb;
    int kmerLen;
    int windowSize;
    float syncmerRatio;

    // For masking reads
    ProbabilityMatrix * probMatrix;
    BaseMatrix * subMat;

    // It converts forward frame to the reverse frame
    const int frameCoversion[3][3] = { //[seq%3][frame%3]
        {3, 5, 4},  // n = 0
        {4, 3, 5},  // n = 1
        {5, 4, 3}   // n = 2
    };

    static constexpr int MAX_N = 32;
    size_t binom[MAX_N + 1][MAX_N + 1];
    



    // Extract query k-mer
    void fillQueryKmerBufferParallel(KSeqWrapper* kseq1,
                                     Buffer<Kmer> &kmerBuffer,
                                     vector<Query> & queryList,
                                     const QuerySplit & currentSplit,
                                     const LocalParameters &par);

    void fillQueryKmerBufferParallel_paired(KSeqWrapper* kseq1,
                                            KSeqWrapper* kseq2,
                                            Buffer<Kmer> &kmerBuffer,
                                            vector<Query> &queryList,
                                            const QuerySplit & currentSplit,
                                            const LocalParameters &par);

    void loadChunkOfReads(KSeqWrapper *kseq,
                          vector<Query> & queryList,
                          size_t & processedQueryNum,
                          size_t chunkSize,
                          size_t chunkEnd,
                          vector<string> & reads,
                          vector<bool> & emptyReads,
                          size_t & count,
                          bool isReverse);

    void processSequence(
        size_t count,
        size_t processedQueryNum,
        const vector<string> & reads,
        const vector<bool> & emptyReads,
        char *seq,
        char *maskedSeq,
        size_t & maxReadLength,
        Buffer<Kmer> &kmerBuffer,
        const vector<Query> & queryList,
        bool isReverse);

    void fillQueryKmerBuffer(
        const char *seq,
        int seqLen, 
        Buffer<Kmer> &kmerBuffer, 
        size_t & posToWrite, 
        uint32_t seqID, 
        uint32_t offset = 0);

    void generatePDMNeighborKmers(
        const char *seq,
        size_t seqStart, 
        int seqLen,
        Buffer<Kmer> &kmerBuffer,
        int threadID,
        size_t & posToWrite,
        int frame,
        uint32_t seqID, 
        uint32_t offset);

    int getKmerCount(
        const char *seq,
        int seqLen);    
    
        int getPDMKmerCount(
    const char *seq,
    int seqLen);

    inline size_t countMutationComb(size_t nt, size_t na, int maxDamage) {
        size_t count = 0;
        size_t N = nt + na;

        for (int d = 1; d <= maxDamage; ++d) {
            if (d > N) break;
            count += binom[N][d];
        }
        return count;
    }
public:
    explicit KmerExtractor(
        const LocalParameters & par,
        const GeneticCode *geneticCode,
        int kmerFormat);

    explicit KmerExtractor(
        const LocalParameters &par,
        const MetamerPattern * metamerPattern);

    ~KmerExtractor();
    
    void extractQueryKmers(
        Buffer<Kmer> &kmerBuffer,
        vector<Query> & queryList,
        const QuerySplit & currentSplit,
        const LocalParameters &par,
        KSeqWrapper* kseq1,
        KSeqWrapper* kseq2 = nullptr);

    bool extractQueryKmers(
        Buffer<Kmer> &kmerBuffer,
        vector<Query> & queryList,
        uint64_t & seqCnt,
        SeqEntry * savedSeq_1,
        SeqEntry * savedSeq_2,
        KSeqWrapper* kseq_1,
        KSeqWrapper* kseq_2 = nullptr);

    bool extractQueryKmers_aa2aa(
        Buffer<Kmer> &kmerBuffer,
        std::vector<ProteinQuery> & queryList,
        KSeqWrapper* kseq,
        uint64_t & processedSeqCnt,
        SeqEntry & savedSeq
    );

    void processSequenceChunk(
        Buffer<Kmer> &kmerBuffer,
        size_t & writePos,
        uint32_t queryOffset,
        const std::vector<std::string> & reads, 
        size_t seqNum,
        char *maskedSeq,
        size_t & maxReadLength,
        const std::vector<std::string> * pairedReads = nullptr
    );

    void processSequenceChunk_aa2aa(
        Buffer<Kmer> &kmerBuffer,
        size_t writePos,
        uint32_t queryOffset,
        const std::vector<std::string> & aaSeqs, 
        size_t seqNum,
        int threadID
    );

    int extractTargetKmers(
        const char *seq,
        Buffer<Kmer> &kmerBuffer,
        size_t &posToWrite,
        int seqID,
        int taxIdAtRank,
        SequenceBlock block);

    bool extractKmers(
        KSeqWrapper *kseq,
        Buffer<Kmer> &kmerBuffer,
        std::unordered_map<string, uint32_t> & accession2index,
        uint32_t & idOffset,
        SeqEntry & savedSeq);
    
    void extractKmer_dna2aa(
        const char *seq,
        int seqLen, 
        Buffer<Kmer> &kmerBuffer, 
        size_t &posToWrite,
        uint32_t seqId1, 
        uint32_t seqId2 = 0);

    bool extractUnirefKmers(
        KSeqWrapper *kseq,
        Buffer<Kmer> &kmerBuffer,
        std::unordered_map<std::string, uint32_t> & unirefName2Id,
        uint32_t & seqCnt,
        SeqEntry & savedSeq);
};

#endif //METABULI_KMEREXTRACTER_H
