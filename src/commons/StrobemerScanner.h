#ifndef METABULI_STROBEMER_SCANNER_H
#define METABULI_STROBEMER_SCANNER_H

#include <algorithm>
#include <cstdlib>
#include <limits>
#include <vector>

#include "KmerScanner.h"


class StrobemerScanner : public MetamerScanner {
protected:
    int strobeNum;
    int strobeLen;
    int windowStart;
    int windowEnd;
    int maxSpan;
    int syncmerSmerLen;
    uint64_t syncmerSmerMask;

    struct StrobeValue {
        uint64_t aa = 0;
        uint64_t dna = 0;
        uint64_t hash = 0;
        bool valid = false;
    };

    std::vector<int> strobeStarts;
    std::vector<StrobeValue> strobes;

    static int checkedKmerSize(const GeneticCode &geneticCode, int strobeNum, int strobeLen) {
        if (strobeNum < 2) {
            std::cerr << "Error: The number of strobes must be at least 2." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (strobeLen < 1) {
            std::cerr << "Error: Strobe length must be positive." << std::endl;
            exit(EXIT_FAILURE);
        }

        const int totalStrobeLen = strobeNum * strobeLen;
        if ((geneticCode.bitPerAA + geneticCode.bitPerCodon) * totalStrobeLen > 64) {
            std::cerr << "Error: Strobemer is too long to fit in 64 bits." << std::endl;
            exit(EXIT_FAILURE);
        }

        return totalStrobeLen;
    }

    static uint64_t mix64(uint64_t x) {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    StrobeValue extractStrobe(int aaStart) const {
        StrobeValue result;
        if (aaStart < 0 || aaStart + strobeLen > aaLen) {
            return result;
        }

        for (int i = 0; i < strobeLen; ++i) {
            int ci;
            int aa;
            int codon;
            if (isForward) {
                ci = seqStart + (aaStart + i) * 3;
                aa = geneticCode.getAA(
                    atcg[seq[ci]],
                    atcg[seq[ci + 1]],
                    atcg[seq[ci + 2]]);
                codon = geneticCode.getCodon(
                    atcg[seq[ci]],
                    atcg[seq[ci + 1]],
                    atcg[seq[ci + 2]]);
            } else {
                ci = seqEnd - (aaStart + i) * 3;
                aa = geneticCode.getAA(
                    iRCT[atcg[seq[ci]]],
                    iRCT[atcg[seq[ci - 1]]],
                    iRCT[atcg[seq[ci - 2]]]);
                codon = geneticCode.getCodon(
                    iRCT[atcg[seq[ci]]],
                    iRCT[atcg[seq[ci - 1]]],
                    iRCT[atcg[seq[ci - 2]]]);
            }

            if (aa < 0) {
                return result;
            }

            result.aa = (result.aa << bitsPerAA) | static_cast<uint64_t>(aa);
            result.dna = (result.dna << bitsPerCodon) | static_cast<uint64_t>(codon);
        }

        result.hash = mix64(result.aa);
        result.valid = true;
        return result;
    }

    int chooseNextStrobeStart(uint64_t chainHash, int previousStart) const {
        uint64_t bestHash = std::numeric_limits<uint64_t>::max();
        int bestStart = -1;

        const int previousLast = previousStart + strobeLen - 1;
        const int begin = previousLast + windowStart;
        const int end = previousLast + windowEnd;
        for (int candidateStart = begin; candidateStart <= end; ++candidateStart) {
            StrobeValue candidate = extractStrobe(candidateStart);
            if (!candidate.valid) {
                continue;
            }

            uint64_t randstrobeHash = mix64(chainHash ^ candidate.hash ^ static_cast<uint64_t>(candidateStart));
            if (randstrobeHash < bestHash) {
                bestHash = randstrobeHash;
                bestStart = candidateStart;
            }
        }

        return bestStart;
    }

    bool buildRandstrobe(int firstStart) {
        strobeStarts[0] = firstStart;
        strobes[0] = extractStrobe(firstStart);
        if (!strobes[0].valid) {
            return false;
        }

        uint64_t chainHash = strobes[0].hash;
        for (int i = 1; i < strobeNum; ++i) {
            int nextStart = chooseNextStrobeStart(chainHash, strobeStarts[i - 1]);
            if (nextStart < 0) {
                return false;
            }

            strobeStarts[i] = nextStart;
            strobes[i] = extractStrobe(nextStart);
            chainHash = mix64(chainHash ^ strobes[i].hash ^ static_cast<uint64_t>(nextStart));
        }

        return true;
    }

    uint64_t packRandstrobe() const {
        uint64_t packedAA = 0;
        uint64_t packedDNA = 0;
        for (int i = 0; i < strobeNum; ++i) {
            packedAA = (packedAA << (bitsPerAA * strobeLen)) | strobes[i].aa;
            packedDNA = (packedDNA << (bitsPerCodon * strobeLen)) | strobes[i].dna;
        }
        return ((packedAA & aaMask) << dnaBits) | (packedDNA & dnaMask);
    }

    bool passesClosedSyncmer(uint64_t packedAA) const {
        if (syncmerSmerLen <= 0) {
            return true;
        }

        uint64_t minSmer = std::numeric_limits<uint64_t>::max();
        int minPos = -1;
        for (int i = 0; i <= kmerSize - syncmerSmerLen; ++i) {
            const uint64_t smer = (packedAA >> (i * bitsPerAA)) & syncmerSmerMask;
            if (smer < minSmer) {
                minSmer = smer;
                minPos = i;
            }
        }

        return minPos == 0 || minPos == (kmerSize - syncmerSmerLen);
    }

    uint8_t packDefaultOffsets() const {
        if (strobeNum != 3 || windowStart != 2 || windowEnd != 4) {
            return 0;
        }

        const int baseDelta = strobeLen - 1 + windowStart;
        const int delta1 = strobeStarts[1] - strobeStarts[0];
        const int delta2 = strobeStarts[2] - strobeStarts[1];
        return static_cast<uint8_t>(((delta1 - baseDelta) & 0x3) |
                                    (((delta2 - baseDelta) & 0x3) << 2));
    }

    void validateParameters() const {
        if (windowStart < 0 || windowEnd < windowStart) {
            std::cerr << "Error: Invalid strobemer extraction window." << std::endl;
            exit(EXIT_FAILURE);
        }
    }

public:
    StrobemerScanner(
        const GeneticCode &geneticCode,
        int strobeNum,
        int strobeLen,
        int windowStart,
        int windowEnd,
        int syncmerSmerLen = 0)
        : MetamerScanner(geneticCode, checkedKmerSize(geneticCode, strobeNum, strobeLen)),
          strobeNum(strobeNum),
          strobeLen(strobeLen),
          windowStart(windowStart),
          windowEnd(windowEnd),
          maxSpan(strobeLen + (strobeNum - 1) * (strobeLen - 1 + windowEnd)),
          syncmerSmerLen(syncmerSmerLen),
          syncmerSmerMask(0)
    {
        validateParameters();
        if (this->syncmerSmerLen < 0 || this->syncmerSmerLen > kmerSize) {
            std::cerr << "Error: Invalid s-mer length for strobemer syncmer selection." << std::endl;
            exit(EXIT_FAILURE);
        }
        if (this->syncmerSmerLen > 0) {
            syncmerSmerMask = (1ULL << (bitsPerAA * this->syncmerSmerLen)) - 1;
        }
        strobeStarts.resize(strobeNum);
        strobes.resize(strobeNum);
    }

    StrobemerScanner(
        const GeneticCode &geneticCode,
        int strobeNum,
        int strobeLen,
        int extractionWindow,
        int syncmerSmerLen = 0)
        : StrobemerScanner(
              geneticCode,
              strobeNum,
              strobeLen,
              strobeLen,
              extractionWindow,
              syncmerSmerLen)
    {}

    int getStrobeNum() const {
        return strobeNum;
    }

    int getStrobeLen() const {
        return strobeLen;
    }

    int getWindowStart() const {
        return windowStart;
    }

    int getWindowEnd() const {
        return windowEnd;
    }

    void initScanner(const char *seq, size_t seqStart, size_t seqEnd, bool isForward = true) override {
        MetamerScanner::initScanner(seq, seqStart, seqEnd, isForward);
        std::fill(strobeStarts.begin(), strobeStarts.end(), 0);
        std::fill(strobes.begin(), strobes.end(), StrobeValue());
    }

    Kmer next() override {
        while (posStart <= aaLen - maxSpan) {
            if (!buildRandstrobe(posStart)) {
                ++posStart;
                continue;
            }

            const uint64_t strobemer = packRandstrobe();
            if (!passesClosedSyncmer(strobemer >> dnaBits)) {
                ++posStart;
                continue;
            }
            const uint8_t strobeOffsets = packDefaultOffsets();
            const int firstStart = posStart++;
            if (isForward) {
                return {strobemer, QueryKmerInfo(0, seqStart + firstStart * 3, 0, strobeOffsets)};
            }

            const int lastStart = strobeStarts.back();
            return {strobemer, QueryKmerInfo(0, seqEnd - (lastStart + strobeLen) * 3 + 1, 0, strobeOffsets)};
        }

        return {UINT64_MAX, 0};
    }
};

#endif //METABULI_STROBEMER_SCANNER_H
