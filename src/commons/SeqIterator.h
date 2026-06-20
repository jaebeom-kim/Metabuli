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
    int kmerLen;
    
public:
    SeqIterator(const LocalParameters &par); 
    
    ~SeqIterator();
        
    void devideToCdsAndNonCds(const char *maskedSeq,
                              size_t seqLen,
                              const vector<CDSinfo> &cdsInfo, 
                              vector<string> &cds,
                              vector<string> &nonCds);
    void getMinHashList(priority_queue<uint64_t> &sortedHashQue, const char *seq);

    bool compareMinHashList(priority_queue<uint64_t> list1, priority_queue<uint64_t> &list2, size_t length1, size_t length2);

    static void maskLowComplexityRegions(const unsigned char * seq, unsigned char * maskedSeq, ProbabilityMatrix & probMat,
                                         float maskProb, const BaseMatrix * subMat);

};

#endif //ADKMER4_KMEREXTRACTOR_H

