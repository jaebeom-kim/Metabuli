#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <cstdint>
#include <iostream>
#include "NcbiTaxonomy.h"
#include <bitset>
#include <sys/types.h>

#define DNA_MASK 0Xffffffff

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



union MetamerID {
    uint64_t id;
    struct {
        uint64_t seqID : 32; // 30?
        uint64_t protID : 32;
    } ids;

    MetamerID(uint32_t seqID, uint32_t protID) : ids{seqID, protID} {}
    MetamerID(uint64_t id) : id(id) {}
    MetamerID() : id(0) {}
};

struct MetamerF { // Deprecated
    MetamerF(uint64_t metamer, uint32_t seqId, uint32_t protId) : metamer(metamer), seqId(seqId), protId(protId) {}
    MetamerF() {}
    uint64_t metamer;
    uint32_t seqId; // taxonomy ID will be stored at last
    uint32_t protId; 
    // The ID of originated CDS is first stored in protId to sort the metamerF list
    // After the CDS is mapped to a protein, the protein ID is stored in protId
};

struct Metamer {
    Metamer() : metamer(0), id(0) {}
    Metamer(uint64_t metamer, uint32_t id) : metamer(metamer), id(id) {}
    uint64_t metamer;
    uint32_t id; // it is mapped to taxonomy ID and protein ID

    Metamer substract(const Metamer & other) const { // self is equal or greater than other
        if (metamer < other.metamer) {
            std::cerr << "Metamer: substract: metamer is smaller than other.metamer" << std::endl;
        } else if (metamer == other.metamer) {
            return Metamer(0, id - other.id);
        } else { // metamer > other.metamer
            if (id > other.id) {
                return Metamer(metamer - other.metamer, id - other.id);
            } else {
                return Metamer((metamer - 1) - other.metamer, UINT32_MAX - other.id + id + 1);
            }
        }
        return Metamer(0, 0);
    }

    static std::bitset<96> substract(const Metamer & metamer1, const Metamer & metamer2) {
        // metamer 1 is the same or greater than metamer2
        if (metamer1.metamer == metamer2.metamer) {
            return std::bitset<96>(metamer1.id - metamer2.id);
        }
        if (metamer1.id >= metamer2.id) {
            std::bitset<96> result;
            result = metamer1.metamer - metamer2.metamer;
            result <<= 32;
            result |= (metamer1.id - metamer2.id);
            return result;
        }
        std::bitset<96> result;
        uint64_t diff = metamer1.metamer - metamer2.metamer - 1;
        result = diff;
        result <<= 32;
        result |= (UINT32_MAX - metamer2.id + metamer1.id + 1);
        return result;    
    }

    Metamer add(const std::bitset<96> & diff) const {
        uint64_t idSum = this->id + (diff & std::bitset<96>(0xFFFFFFFF)).to_ullong();
        uint64_t metamerSum = this->metamer + (diff >> 32).to_ullong() + (idSum >> 32);
        idSum &= 0xFFFFFFFF;
        return Metamer(metamerSum, idSum);
    }

    bool operator < (const Metamer & other) const {
        if (metamer != other.metamer) {
            return metamer < other.metamer;
        }
        return id < other.id;
    }

    bool operator == (const Metamer & other) const {
        return metamer == other.metamer && id == other.id;
    }
};

struct DeltaIdxOffset{
    DeltaIdxOffset(Metamer metamer, size_t offset) : metamer(metamer), offset(offset) { }
    DeltaIdxOffset() {};
    Metamer metamer;
    size_t offset;
};

struct ExtractedMetamer {
    ExtractedMetamer(uint64_t metamer, uint32_t id, TaxID speciesId, uint32_t cdsPos, uint32_t unirefId) 
     : metamer(metamer, id), speciesId(speciesId), cdsPos(cdsPos), unirefId(unirefId) {}
    Metamer metamer;
    TaxID speciesId;
    uint32_t cdsPos;
    uint32_t unirefId;
};

struct TargetMetamerF { // Deprecated
    TargetMetamerF(uint64_t metamer, uint32_t seqId, uint32_t protId, TaxID taxId, uint32_t cdsPos) 
        : metamerF(metamer, seqId, protId), speciesId(taxId), cdsPos(cdsPos) {}
    MetamerF metamerF;
    TaxID speciesId;
    uint32_t cdsPos; // The position of the metamer in the CDS
};


#endif //ADKMER3_KMER_H