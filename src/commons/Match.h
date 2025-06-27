#ifndef ADCLASSIFIER2_MATCH_H
#define ADCLASSIFIER2_MATCH_H

#include "Kmer.h"
#include <cstdint>
#include <iostream>
#include "BitManipulateMacros.h"

struct Match { // 24 byte
    Match(){}
    Match(QueryKmerInfo qInfo,
          int targetId,
          TaxID speciesId,
          uint32_t dnaEncoding,
          uint16_t eachHamming,
          uint8_t hamming):
          qInfo(qInfo), targetId(targetId), speciesId(speciesId), dnaEncoding(dnaEncoding),
          rightEndHamming(eachHamming), hamming(hamming) { }

    QueryKmerInfo qInfo;      // 8 // Query K-mer information
    TaxID targetId;           // 4 // axonomy id infact
    TaxID speciesId;          // 4 // Used to group matches by species
    uint32_t dnaEncoding;     // 4 // Used to check if two matches are consecutive
    uint16_t rightEndHamming; // 2 // Used to calculate score
    uint8_t hamming;          // 1 // Used to filter redundant matches

    void printMatch() const {
        std::cout << qInfo.sequenceID << " " << qInfo.pos << " " << qInfo.frame << " "
        << targetId << " " << speciesId << " " << rightEndHamming << " " << (int)hamming << " " << getScore() << std::endl;
    }

    float getScore(float score = 0.0f, int cnt = 0) const { 
        int currentHamming = GET_2_BITS(rightEndHamming >> (cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        if (cnt == 7) {
            return score;
        } else {
        return getScore(score, cnt + 1);    
        }
    }

    bool isConsecutive_DNA(const Match * match2) const {
        // match1 87654321 -> 08765432
        // match2 98765432 -> 08765432
        return (dnaEncoding >> 3) == (match2->dnaEncoding & 0x1FFFFF);
    }

    float getRightPartScore(const int range, float score = 0.0f, int cnt = 0) const {
        if (cnt == range) {
            return score;
        }
        int currentHamming = GET_2_BITS(rightEndHamming >> (cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        return getRightPartScore(range, score, cnt + 1);    
    }

    float getLeftPartScore(const int range, float score = 0.0f, int cnt = 0) const {
        if (cnt == range) {
            return score;
        }
        int currentHamming = GET_2_BITS(rightEndHamming >> (14 - cnt * 2));
        if (currentHamming == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * currentHamming;
        }
        return getLeftPartScore(range, score, cnt + 1);    
    }

    int getRightPartHammingDist(const int range) const {
        int sum = 0;
        for (int i = 0; i < range; i++) {
            sum += GET_2_BITS(rightEndHamming >> (i * 2));
        }
        return sum;
    }

    int getLeftPartHammingDist(const int range) const {
        int sum = 0;
        for (int i = 0; i < range; i++) {
            sum += GET_2_BITS(rightEndHamming >> (14 - i * 2));
        }
        return sum;
    }
};

#endif //ADCLASSIFIER2_MATCH_H