#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include "Kmer.h"
#include <cstdint>
#include <iostream>
#include "BitManipulateMacros.h"
class Match { // 32 byte
public:
    Match(){}
    Match(Kmer qKmer, Kmer tKmer): qKmer(qKmer), tKmer(tKmer) { }

    Kmer qKmer;
    Kmer tKmer;

    void printMatch() const {
        std::cout << qKmer.qInfo.sequenceID << " " << qKmer.qInfo.pos << " " << qKmer.qInfo.frame << " "
        << tKmer.tInfo.taxId << " " << tKmer.tInfo.speciesId << "\n";
    }

    static bool compare(const Match &a, const Match &b) {
        if (a.qKmer.qInfo.sequenceID != b.qKmer.qInfo.sequenceID)
            return a.qKmer.qInfo.sequenceID < b.qKmer.qInfo.sequenceID;
        
        if (a.tKmer.tInfo.speciesId != b.tKmer.tInfo.speciesId)
            return a.tKmer.tInfo.speciesId < b.tKmer.tInfo.speciesId;

        if (a.qKmer.qInfo.frame != b.qKmer.qInfo.frame)
            return a.qKmer.qInfo.frame < b.qKmer.qInfo.frame;

        if (a.qKmer.qInfo.pos != b.qKmer.qInfo.pos)
            return a.qKmer.qInfo.pos < b.qKmer.qInfo.pos;

        if (a.tKmer.tInfo.taxId != b.tKmer.tInfo.taxId)
            return a.tKmer.tInfo.taxId < b.tKmer.tInfo.taxId;

        if (a.qKmer.value != b.qKmer.value)
            return a.qKmer.value < b.qKmer.value;
        
        return a.tKmer.value < b.tKmer.value;
    }

};


struct Match_AA {
    uint32_t queryId;
    uint32_t targetId;
    uint32_t pos;     // For developing purpose only
    uint64_t kmer;    // For developing purpose only

    Match_AA(uint32_t queryId, uint32_t targetId) : queryId(queryId), targetId(targetId) { }

    Match_AA(uint32_t queryId, uint32_t targetId, uint64_t kmer) 
        : queryId(queryId), targetId(targetId), kmer(kmer) { }
    
    Match_AA(uint32_t queryId, uint32_t targetId, uint32_t pos, uint64_t kmer) 
        : queryId(queryId), targetId(targetId), pos(pos), kmer(kmer) { }

    static bool compare(const Match_AA &a, const Match_AA &b) {
        if (a.queryId != b.queryId)
            return a.queryId < b.queryId;
        if (a.pos != b.pos)
            return a.pos < b.pos;
        return a.targetId < b.targetId;
    }
};


struct MatchBlock {
    MatchBlock(size_t start, size_t end, int id) : start(start), end(end), id(id) {}
    MatchBlock() : start(0), end(0), id(0) {}
    size_t start;
    size_t end;
    uint32_t id;
};


#endif //ADCLASSIFIER2_MATCH_H