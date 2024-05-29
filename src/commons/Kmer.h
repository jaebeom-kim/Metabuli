#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
#include "NcbiTaxonomy.h"

struct QueryKmerInfo {
    explicit QueryKmerInfo(uint32_t seqID = 0, uint32_t pos = 0, uint8_t frame = 0 ) : pos(pos), sequenceID(seqID), frame(frame) {}
    uint64_t pos : 32;
    uint64_t sequenceID : 29;
    uint64_t frame : 3; // 0, 1, 2 are forward, and 3, 4, 5 are reverse 1 byte
}; // 8 byte


typedef struct QueryKmer {
    QueryKmer(uint64_t ADkmer, uint32_t seqID, uint32_t pos, uint8_t frame) : ADkmer(ADkmer), info(seqID, pos, frame) {}
    QueryKmer():ADkmer(0), info(0,0,0){}
    uint64_t ADkmer; // 8 byte
    QueryKmerInfo info; // 8 byte
} QueryKmer; // 16 byte


struct TargetKmerInfo{
    explicit TargetKmerInfo(int seqID = 0, bool redundancy = false) : sequenceID(seqID), redundancy(redundancy) {}
    int sequenceID : 31;
    int redundancy : 1;
    bool operator == (const TargetKmerInfo & info) const{
        return (sequenceID == info.sequenceID && this->redundancy==info.redundancy);
    }
}; // 4 bytes

struct TargetKmer{
    TargetKmer(): info(0, false), taxIdAtRank(0), ADkmer(0)  { };
    TargetKmer(uint64_t ADkmer, TaxID taxIdAtRank, int seqID, bool redundacy)
        : info(seqID, redundacy), taxIdAtRank(taxIdAtRank), ADkmer(ADkmer) {}
    TargetKmerInfo info; // 4 byte
    TaxID taxIdAtRank; // 4 byte
    uint64_t ADkmer; // 8 byte
};

struct DiffIdxSplit{
    DiffIdxSplit(uint64_t ADkmer, size_t diffIdxOffset, size_t infoIdxOffset) : ADkmer(ADkmer), diffIdxOffset(diffIdxOffset), infoIdxOffset(infoIdxOffset) { }
    DiffIdxSplit(const DiffIdxSplit & copy) {ADkmer = copy.ADkmer; diffIdxOffset = copy.diffIdxOffset; infoIdxOffset=copy.infoIdxOffset;}
    DiffIdxSplit() {};
    uint64_t ADkmer;
    size_t diffIdxOffset;
    size_t infoIdxOffset;
};

struct AAKmer {
    uint64_t kmer : 36;
    uint64_t id : 28;
};
struct ProtIdxSplit{
    ProtIdxSplit(uint64_t kmer, size_t idxOffset) : kmer(kmer), idxOffset(idxOffset) { }
    ProtIdxSplit() {};
    uint64_t kmer; 
    size_t idxOffset;
    // getNextKmer(kmer, idxOffset) = next kmer
};

struct MetamerF {
    MetamerF(uint64_t metamer, uint32_t seqId, uint32_t protId) : metamer(metamer), seqId(seqId), protId(protId) {}
    MetamerF() {}
    uint64_t metamer;
    uint32_t seqId;
    uint32_t protId; 
    // The ID of originated CDS is first stored in protId to sort the metamerF list
    // After the CDS is mapped to a protein, the protein ID is stored in protId
     
};

struct TargetMetamerF {
    TargetMetamerF(uint64_t metamer, uint32_t seqId, uint32_t protId, TaxID taxId, uint32_t cdsPos) 
        : metamerF(metamer, seqId, protId), speciesId(taxId), cdsPos(cdsPos) {}
    MetamerF metamerF;
    TaxID speciesId;
    uint32_t cdsPos; // The position of the metamer in the CDS
};


#endif //ADKMER3_KMER_H
