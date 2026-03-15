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

    ~Match() {}

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

class MatchWithPos {
    public:
    Kmer qKmer;
    Kmer tKmer;
    uint16_t posId;
    MatchWithPos(Kmer qKmer, Kmer tKmer, uint16_t posId) : qKmer(qKmer), tKmer(tKmer), posId(posId) {}

    void printMatch() const {
        std::cout << qKmer.qInfo.sequenceID << " " << qKmer.qInfo.pos << " " << qKmer.qInfo.frame << " "
        << tKmer.tInfo.taxId << " " << tKmer.tInfo.speciesId << " " << posId << "\n";
    }  

    static bool compare(const MatchWithPos &a, const MatchWithPos &b) {
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


struct MatchScore {
    float idScore;
    float subScore;
    float logE;
    double logP; // 2 * sigma(logPi)

    MatchScore() : idScore(0.0f), subScore(0.0f), logE(0.0f), logP(0.0) {}
    MatchScore(float idScore, float subScore) : idScore(idScore), subScore(subScore), logE(0.0f), logP(0.0) {}
    MatchScore(float idScore, float subScore, double logP) : idScore(idScore), subScore(subScore), logE(0.0f), logP(logP) {}

    void print() const {
        std::cout << "ID Score: " << idScore << ", Substitution Score: " << subScore << ", logE: " << logE << ", logP: " << logP << std::endl;
    }

    bool isLargerThan(const MatchScore & other, int mode) const {
        if (mode == 0) {
            return idScore > other.idScore;
        } else if (mode == 1) {
            return subScore > other.subScore;
        } else {
            return (idScore + subScore) > (other.idScore + other.subScore);
        }
    }

    MatchScore operator*(float factor) const {
        return MatchScore(idScore * factor, subScore * factor);
    }

    MatchScore & operator+=(const MatchScore & other) {
        idScore += other.idScore;
        subScore += other.subScore;
        logP += other.logP;
        return *this;
    }

    MatchScore & operator-=(const MatchScore & other) {
        idScore -= other.idScore;
        subScore -= other.subScore;
        logP -= other.logP;
        return *this;
    }

};

inline MatchScore operator+(MatchScore lhs, const MatchScore& rhs) {
    lhs += rhs;
    return lhs;
}


template <typename MatchType>
struct MatchPath {
    MatchPath() : start(0), end(0), score(), hammingDist(0), depth(0),  historyMask(0), startMatch(nullptr), endMatch(nullptr), prevMatchIdx(-1) {}

    MatchPath(int start, int end, MatchScore score, int hammingDist, int depth, const MatchType * startMatch, const MatchType * endMatch) :
         start(start), end(end), score(score), hammingDist(hammingDist), depth(depth),  historyMask(0), startMatch(startMatch), endMatch(endMatch), prevMatchIdx(-1) {}

    
    MatchPath(const MatchType * startMath, int windowSizeNt) 
        : start(startMath->qKmer.qInfo.pos),
          end(startMath->qKmer.qInfo.pos + windowSizeNt - 1), 
          score(),
          hammingDist(0),
          depth(1),
          historyMask(0),
          startMatch(startMath),
          endMatch(startMath),
          prevMatchIdx(-1) {}
    
    MatchPath(const MatchType * startMatch, MatchScore score, int hammingDist, int kmerLenNt) 
        : start(startMatch->qKmer.qInfo.pos),
          end(startMatch->qKmer.qInfo.pos + kmerLenNt - 1),
          score(score),
          hammingDist(hammingDist),
          depth(1),
          historyMask(0),
          startMatch(startMatch),
          endMatch(startMatch),
          prevMatchIdx(-1) {}
    
    int start;                // query coordinate
    int end;                  // query coordinate
    MatchScore score;
    int hammingDist;
    int depth;
    uint32_t historyMask;         // For spaced pattern with history

    const MatchType * startMatch; // Start match of the path
    const MatchType * endMatch;   // It is always the current match
    int prevMatchIdx;
    std::vector<const MatchType*> chain;

    void printMatchPath() {
        std::cout << start << " " << end << " " << score.idScore << " " << score.subScore << " " << hammingDist << " " << depth << std::endl;
    }
};

#endif //ADCLASSIFIER2_MATCH_H