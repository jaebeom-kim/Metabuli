#ifndef ADKMER4_KMEREXTRACTOR_H
#define ADKMER4_KMEREXTRACTOR_H

#include <iostream>
#include <unordered_map>
#include <vector>
#include "Kmer.h"
#include "printBinary.h"
#include "common.h"
#include "Mmap.h"
#include "KmerBuffer.h"
#include <algorithm>
#include "kseq.h"
#include "KSeqBufferReader.h"
#include "ProdigalWrapper.h"
#include <functional>
//#include "xxh3.h"
#include "xxhash.h"
#include <queue>
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"

#ifdef OPENMP
#include <omp.h>
#endif

KSEQ_INIT(kseq_buffer_t*, kseq_buffer_reader)

#define kmerLength 8

#define nuc2int(x) (x & 14u)>>1u

using namespace std;

typedef struct PredictedBlock {
    PredictedBlock(int start, int end, int strand) : start(start), end(end), strand(strand) {}

    void printPredictedBlock() {
        cout << strand << " " << start << " " << end << endl;
    }

    int start;
    int end;
    int strand; //true for forward
} PredictedBlock;

class SeqIterator {
private:
    static const string iRCT;
    static const string atcg;
    int * aa2num;
    vector<int> aaFrames[6];
    uint64_t powers[10];
    int nuc2aa[8][8][8];
    uint64_t nuc2num[4][4][4];
    unordered_map<int, char> num2aa;
    uint32_t * mask;
    int * mask_int;
    uint32_t spaceNum;
    int spaceNum_int;
    int bitsForCodon;
    int bitsFor8Codons;


    void addDNAInfo_QueryKmer(uint64_t &kmer, const char *seq, int forOrRev, uint32_t kmerCnt, uint32_t frame,
                              int readLength);

    void addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, const PredictedBlock &block, const int &kmerCnt);
    void addDNAInfo_TargetKmer(uint64_t &kmer, const char *seq, int kmerCnt, int frame);

public:
    
  void devideToCdsAndNonCds(const char *seq, size_t seqLen,
                            const vector<CDSinfo> &cdsInfo, vector<string> &cds,
                            vector<string> &nonCds);

  void devideToCdsAndNonCds(const char *seq, size_t seqLen,
                            ProdigalWrapper * prodigal, vector<string> &cds,
                            vector<string> &nonCds);
                            
  void fillQueryKmerBuffer(const char *seq, int seqLen,
                           Buffer<QueryKmer> &kmerBuffer, size_t &posToWrite,
                           uint32_t seqID, uint32_t offset = 0);

  string reverseComplement(string &read) const;

  char *reverseComplement(const char *read, size_t length) const;

  void sixFrameTranslation(const char *seq, int seqLen);

  bool translateBlock(const char *seq, PredictedBlock block);

  bool translate(const string & cds, int frame = 0);

  bool translate(const char * cds, int frame = 0);

  bool translate(string & aa, const string & cds);

  void generateIntergenicKmerList(struct _gene *genes, struct _node *nodes,
                                  int numberOfGenes,
                                  vector<uint64_t> &intergenicKmerList,
                                  const char *seq);

  void getExtendedORFs(struct _gene *genes, struct _node *nodes,
                       vector<PredictedBlock> &blocks, size_t numOfGene,
                       size_t length, size_t &numOfBlocks,
                       vector<uint64_t> &intergenicKmerList, const char *seq);

  static void getMinHashList(priority_queue<uint64_t> &sortedHashQue, const char *seq, size_t listSize = 3000);

  bool compareMinHashList(priority_queue<uint64_t> list1,
                          priority_queue<uint64_t> &list2, size_t length1,
                          size_t length2);

  static size_t kmerNumOfSixFrameTranslation(const char *seq);

  size_t getNumOfKmerForBlock(const PredictedBlock &block);

  int fillBufferWithKmerFromBlock(const PredictedBlock &block, const char *seq,
                                  Buffer<TargetKmer> &kmerBuffer,
                                  size_t &posToWrite, int seqID,
                                  int taxIdAtRank);

  int computeMetamers(const char * seq, int frame, Buffer<TargetKmer> & kmerBuffer, size_t & posToWrite, int seqID, int taxIdAtRank);
  int computeMetamerF(const char * seq, int frame, Buffer<TargetMetamerF> & kmerBuffer, size_t & posToWrite, uint32_t seqID, int taxIdAtRank, uint32_t cdsIdx);
  int computeMetamerF(const char * seq, int frame, Buffer<ExtractedMetamer> & kmerBuffer, size_t & posToWrite, int taxIdAtRank, uint32_t protId, uint32_t unirefId);
  size_t computeAAKmer(Buffer<uint64_t> * kmerBuffer, const char * seq, size_t seqLength, size_t & posToWrite, size_t seqId); 

  static void maskLowComplexityRegions(const char *seq, char *maskedSeq,
                                       ProbabilityMatrix &probMat,
                                       float maskProb,
                                       const BaseMatrix *subMat);

  void printKmerInDNAsequence(uint64_t kmer);
  void printAAKmer(uint64_t kmer, int shits = 28);
  void printTranslation(const string & dna);

  explicit SeqIterator(const LocalParameters &par);
  ~SeqIterator();
};

#endif //ADKMER4_KMEREXTRACTOR_H

