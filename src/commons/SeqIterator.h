#ifndef ADKMER4_KMEREXTRACTOR_H
#define ADKMER4_KMEREXTRACTOR_H

#include <iostream>
#include <vector>
#include <queue>
#include <cstdint>
#include <algorithm>
#include <functional>

#include "common.h"
#include "Mmap.h"
#include "xxhash.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "LocalUtil.h"
#include "SyncmerScanner.h"
#include "Kmer.h"
#include "printBinary.h"
#include "kseq.h"
#include "KSeqBufferReader.h"
#include "ProdigalWrapper.h"
#ifdef OPENMP
    #include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define nuc2int(x) (x & 14u)>>1u

using namespace std;

class SeqIterator {
private:
    uint64_t * powers;
    uint32_t * mask;
    int * mask_int;
    uint32_t spaceNum;
    int spaceNum_int;
    int bitsForCodon;
    uint32_t smerMask;
    uint64_t dnaMask;
    int kmerLen = 12;
    int smerLen;
    
public:
    SeqIterator(const LocalParameters &par); 
    
    ~SeqIterator();
    
    string reverseComplement(string &read) const;
    
    void devideToCdsAndNonCds(const char *maskedSeq,
                              size_t seqLen,
                              const vector<CDSinfo> &cdsInfo, 
                              vector<string> &cds,
                              vector<string> &nonCds);

    char *reverseComplement(char *read, size_t length) const;

    void generateIntergenicKmerList(struct _gene *genes, struct _node *nodes, int numberOfGenes,
                                    vector<uint64_t> &intergenicKmerList, const char *seq);

    void getMinHashList(priority_queue<uint64_t> &sortedHashQue, const char *seq);

    bool compareMinHashList(priority_queue<uint64_t> list1, priority_queue<uint64_t> &list2, size_t length1, size_t length2);

    bool isSyncmer(const vector<int> &aaSeq, int startPos, int k, int s) {
        size_t min_smer_value = UINT64_MAX;
        int min_smer_pos = -1;
        size_t current_value = 0;
        for (int i = 0; i <= k - s; ++i) {
            if (i == 0) {
                for (int j = 0; j < s; ++j) {
                    current_value = (current_value << 5) | aaSeq[startPos + i + j];
                }
            } else {
                current_value = (current_value << 5) | aaSeq[startPos + i + s - 1];
            }
            current_value = current_value & ((1ULL << (5 * s)) - 1); // Mask to keep only the last s amino acids
            if (current_value < min_smer_value) {
                min_smer_value = current_value;
                min_smer_pos = i;
            }
        }
        return (min_smer_pos == 0 || min_smer_pos == (k - s));
    }

    static void maskLowComplexityRegions(const unsigned char * seq, unsigned char * maskedSeq, ProbabilityMatrix & probMat,
                                         float maskProb, const BaseMatrix * subMat);

};

#endif //ADKMER4_KMEREXTRACTOR_H

