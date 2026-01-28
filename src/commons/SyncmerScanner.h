#ifndef METABULI_SYNCMER_SCANNER_H
#define METABULI_SYNCMER_SCANNER_H

#include <iostream>
#include <deque>

#include "KmerScanner.h"
#include "printBinary.h"

class SyncmerScanner : public MetamerScanner {
protected:
    // Internal values
    int smerLen;
    uint64_t smerMask;

    // Variables for syncmer scanning
    std::deque<Kmer> dq;
    int smerCnt;
    uint64_t smer;
    int prevPos;

public:
    SyncmerScanner(int kmerLen, int smerLen, const GeneticCode &geneticCode) 
        : MetamerScanner(geneticCode, kmerLen) 
    {
        this->smerLen = smerLen;
        this->smerMask = (1ULL << (5 * smerLen)) - 1;
    }

    void initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward) override {
        MetamerScanner::initScanner(seq, seqStart, seqEnd, isForward);
        this->dq.clear();
        this->smerCnt = 0;
        this->smer = 0;
        this->prevPos = -kmerSize;
    }

    Kmer next() override {
        bool syncmerFound = false;
        int aa = 0;
        while (posStart <= aaLen - kmerSize && !syncmerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            while (smerCnt < kmerSize - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    if (isForward) {
                        int ci = seqStart + (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    } else {
                        int ci = seqEnd - (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    }
                    if (aa < 0) { sawN = true; break; }
                    smer = (smer << bitsPerAA) | (uint64_t)aa;                    
                    loadedCharCnt++;
                }
                if (sawN) break;
                smer &= smerMask;
                while (!dq.empty() && dq.back().value > smer) dq.pop_back();
                dq.emplace_back(smer, posStart + smerCnt);
                smerCnt++;
            }
            if (sawN) {
                posStart += smerCnt + loadedCharCnt + 1;
                prevPos = posStart - kmerSize; // Reset previous position ??
                dq.clear(); 
                smerCnt = loadedCharCnt = 0; 
                smer = 0;
                continue;
            }
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            uint32_t anchor1 = static_cast<uint32_t>(posStart);
            uint32_t anchor2 = static_cast<uint32_t>(posStart + (kmerSize - smerLen));
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                if (isForward) {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqStart + (prevPos + kmerSize + i) * 3;
                        aaPart = (aaPart << bitsPerAA) | (uint64_t)geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                        dnaPart = (dnaPart << bitsPerCodon) | (uint64_t)geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    }
                } else {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqEnd - (prevPos + kmerSize + i) * 3;
                        aaPart = (aaPart << bitsPerAA) | (uint64_t)geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                        dnaPart = (dnaPart << bitsPerCodon) | (uint64_t)geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    }
                }
                prevPos = posStart;
                syncmerFound = true;
            }
            ++posStart;
        }
        if (syncmerFound) {
            if (isForward) {
                return {(aaPart << dnaBits) | (dnaPart & dnaMask), seqStart + prevPos * 3};
            } else {
                return {(aaPart << dnaBits) | (dnaPart & dnaMask), seqEnd - (prevPos + kmerSize) * 3 + 1};
            }
        } else {
            return {UINT64_MAX, 0}; // No more syncmers found
        }
    }
};

class SyncmerScanner_aa2aa : public KmerScanner_aa2aa {
protected:
    // Internal values
    int smerLen;
    uint64_t smerMask;

    // Variables for syncmer scanning
    std::deque<Kmer> dq;
    int smerCnt;
    uint64_t smer;
    int prevPos;

public:
    SyncmerScanner_aa2aa(int k, int s) : KmerScanner_aa2aa(k), smerLen(s) {
        smerMask = (1ULL << (5 * smerLen)) - 1;
    }

    ~SyncmerScanner_aa2aa() {}

    void initScanner(
        const char * seq, 
        size_t seqStart, 
        size_t seqEnd, 
        bool isForward = true) override 
    {
        KmerScanner_aa2aa::initScanner(seq, seqStart, seqEnd, isForward);
        this->dq.clear();
        this->smerCnt = 0;
        this->smer = 0;
        this->prevPos = -kmerSize;
    }

    Kmer next() override {
        bool syncymerFound = false;
        int aa = 0;
        while ((posStart <= seqLen - kmerSize) && !syncymerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            // Fill the deque with s-mers
            while (smerCnt < kmerSize - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    aa = aacids[seq[seqStart + posStart + smerCnt + loadedCharCnt]];
                    if (aa > 23) { sawN = true; break; }
                    smer = (smer << 5) | (uint64_t)aa;
                    loadedCharCnt++;
                }
                if (sawN) break;
                smer &= smerMask;
                while (!dq.empty() && dq.back().value > smer) dq.pop_back();
                dq.emplace_back(smer, posStart + smerCnt);
                smerCnt++;
            }
            if (sawN) {
                posStart += smerCnt + loadedCharCnt + 1;
                prevPos = posStart - kmerSize;
                loadedCharCnt = 0;
                smerCnt = 0;
                smer = 0;
                dq.clear();
                continue;
            }

            // Remove s-mers that are out of the current k-mer window
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            uint32_t anchor1 = static_cast<uint32_t>(posStart);
            uint32_t anchor2 = static_cast<uint32_t>(posStart + (kmerSize - smerLen));

            // Check if the minimum s-mer is at one of the anchor positions
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shift = posStart - prevPos;
                for (int i = 0; i < shift; ++i) {
                    int aa = aacids[seq[seqStart + prevPos + kmerSize + i]];
                    aaPart = (aaPart << 5) | (uint64_t)aa;
                }
                prevPos = posStart;
                syncymerFound = true;
            }
            ++posStart;
        }
        if (syncymerFound) {
            return { aaPart & mask, seqStart + (posStart - 1) };
        }
        return { UINT64_MAX, 0 }; // No more kmers found
    }
};

class SyncmerScanner_dna2aa : public KmerScanner_dna2aa {
protected:
    // Internal values
    int smerLen;
    uint64_t smerMask;

    // Variables for syncmer scanning
    std::deque<Kmer> dq;
    int smerCnt;
    uint64_t smer;
    int prevPos;

public:
    SyncmerScanner_dna2aa(const GeneticCode &geneticCode, int k, int s) 
        : KmerScanner_dna2aa(geneticCode, k), smerLen(s) 
    {
        this->smerMask = (1ULL << (bitsPerAA * smerLen)) - 1;
    }

    void initScanner(
        const char * seq, 
        size_t seqStart, 
        size_t seqEnd, 
        bool isForward = true) override 
    {
        KmerScanner_dna2aa::initScanner(seq, seqStart, seqEnd, isForward);
        this->dq.clear();
        this->smerCnt = 0;
        this->smer = 0;
        this->prevPos = -kmerSize;
    }

    Kmer next() override {
        bool syncymerFound = false;
        int aa = 0;
        while (posStart <= aaLen - kmerSize && !syncymerFound) {
            bool sawN = false;
            smerCnt -= (smerCnt > 0);
            // Fill the deque with s-mers
            while (smerCnt < kmerSize - smerLen + 1) {
                loadedCharCnt -= (loadedCharCnt == smerLen);
                while (loadedCharCnt < smerLen) {
                    if (isForward) {
                        int ci = seqStart + (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    } else {
                        int ci = seqEnd - (posStart + smerCnt + loadedCharCnt) * 3;
                        aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);            
                    }
                    if (aa < 0) { sawN = true; break; }
                    smer = (smer << bitsPerAA) | (uint64_t)aa;
                    loadedCharCnt++;
                }
                if (sawN) break;
                smer &= smerMask;
                while (!dq.empty() && dq.back().value > smer) dq.pop_back();
                dq.emplace_back(smer, posStart + smerCnt);
                smerCnt++;
            }
            if (sawN) {
                posStart += smerCnt + loadedCharCnt + 1;
                prevPos = posStart - kmerSize;
                loadedCharCnt = 0;
                smerCnt = 0;
                smer = 0;
                dq.clear();
                continue;
            }

            // Remove s-mers that are out of the current k-mer window
            if (!dq.empty() && dq.front().pos < posStart) dq.pop_front();
            uint32_t anchor1 = static_cast<uint32_t>(posStart);
            uint32_t anchor2 = static_cast<uint32_t>(posStart + (kmerSize - smerLen));

            // Check if the minimum s-mer is at one of the anchor positions
            if (!dq.empty() && (dq.front().pos == anchor1 || dq.front().pos == anchor2)) {
                int shifts = posStart - prevPos;
                if (isForward) {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqStart + (prevPos + kmerSize + i) * 3;
                        aaPart = (aaPart << bitsPerAA) | (uint64_t)geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    }
                } else {
                    for (int i = 0; i < shifts; ++i) {
                        int ci = seqEnd - (prevPos + kmerSize + i) * 3;
                        aaPart = (aaPart << bitsPerAA) | (uint64_t)geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    }
                }
                prevPos = posStart;
                syncymerFound = true;
            }
            ++posStart;
        }
        if (syncymerFound) {
            if (isForward) {
                return { aaPart & mask, seqStart + (prevPos * 3)};
            } else {
                return { aaPart & mask, seqEnd - (prevPos + kmerSize) * 3 + 1};
            }
        }
        return {UINT64_MAX, 0}; // No more kmers found
    }

};


class SpacedSyncmerScanner : public SpacedMetamerScanner {
    // Syncmer settings
    int smerLen;
    uint64_t smerMask;  // Mask for s-mer (s * bitsPerAA)
    
public:
    SpacedSyncmerScanner(
        const GeneticCode& gc, 
        int kmerLen,
        uint32_t spaceMask,
        int sLen) 
        : SpacedMetamerScanner(gc, kmerLen, spaceMask), smerLen(sLen) 
    {
        smerMask = (~0ULL) >> (64 - (smerLen * bitsPerAA));
    }

    Kmer next() override {
        int aa = 0;
        int codon = 0;

        while (posStart <= aaLen - windowSize) {
            bool sawN = false;
            
            // --- Phase 1: Load Contiguous Window ---
            loadedCharCnt -= (loadedCharCnt == windowSize);
            while (loadedCharCnt < windowSize) {
                int ci;
                if (isForward) {
                    ci    = seqStart + (posStart + loadedCharCnt) * 3;
                    aa    = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    codon = geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                } else {
                    ci    = seqEnd - (posStart + loadedCharCnt) * 3;
                    aa    = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                    codon = geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                }

                if (aa < 0) { sawN = true; break; }
                
                // Rolling buffer updates
                aaPart  = (aaPart  << bitsPerAA)    | (uint64_t)aa;
                dnaPart = (dnaPart << bitsPerCodon) | (uint64_t)codon;
                
                loadedCharCnt++;
            }

            if (sawN) {
                posStart += loadedCharCnt + 1;
                loadedCharCnt = 0;
                aaPart = dnaPart = 0;
                continue;
            }

            // --- Phase 2: Construct Packed AA K-mer ---
            // We need the packed AA version to check the syncmer condition
            uint64_t packedAA = 0;
            
            // OPTIMIZATION: If AVX2 is available, use _pext_u64 here!
            // packedAA = _pext_u64(aaPart, precomputedAaPextMask);
            
            // Portable fallback:
            uint32_t validPosMask = spaceMask;
            while (validPosMask) {
                int i = __builtin_ctz(validPosMask);
                validPosMask &= (validPosMask - 1);
                const uint64_t curAA = (aaPart >> (i * bitsPerAA)) & aaMaskPerCodon;
                packedAA |= (curAA << aaShifts[i]);
            }

            // --- Phase 3: Check "True Spaced Syncmer" ---
            // Find the minimum s-mer within the PACKED k-mer
            uint64_t minSmer = UINT64_MAX;
            int minPos = -1;
            
            // We iterate over the ACTIVE positions (0 to kmerLen - s)
            for (int i = 0; i <= kmerSize - smerLen; ++i) {
                // Extract s-mer from packed sequence
                const uint64_t sVal = (packedAA >> (i * bitsPerAA)) & smerMask;
                
                // Simple hash mixing to randomize AA order (optional but recommended)
                // sVal = hash64(sVal); 

                if (sVal < minSmer) {
                    minSmer = sVal;
                    minPos = i;
                }
            }

            // Check if min s-mer is at the target position (e.g., Start or End)
            // Here we check if it is at the "start" (LSB or MSB depending on your definition)
            // Assuming "Open Syncmer" at index 0 (relative to packed form)
            // bool isSelected = (minPos == 0); 
            
            // For "Closed Syncmer" (more sensitive):
            const bool isSelected = (minPos == 0 || minPos == (kmerSize - smerLen));

            if (isSelected) {
                // --- Phase 4: Construct Packed DNA K-mer ---
                // Only now do we expend the effort to pack the DNA part
                uint64_t packedDNA = 0;
                validPosMask = spaceMask;
                while (validPosMask) {
                    int i = __builtin_ctz(validPosMask);
                    validPosMask &= (validPosMask - 1);
                    const uint64_t curCodon = (dnaPart >> (i * bitsPerCodon)) & codonMask;
                    packedDNA |= (curCodon << dnaShifts[i]);
                }
                // int packedIdx = 0;  
                // for (int i = 0; i < windowSize; ++i) {
                //     if (spaceMask & (1u << i)) {
                //         uint64_t curCodon = (dnaPart >> (i * bitsPerCodon)) & codonMask;
                //         packedDNA |= (curCodon << (packedIdx * bitsPerCodon));
                //         packedIdx++;
                //     }
                // }

                // Prepare result
                const uint64_t finalKmer = (packedAA << dnaBits) | packedDNA;
                uint32_t finalPos;
                
                if (isForward) finalPos = seqStart + posStart * 3;
                else           finalPos = seqEnd - (posStart + windowSize) * 3 + 1;

                posStart++; 
                return { finalKmer, finalPos };
            }

            // Not selected, move to next
            posStart++;
        }
        
        return { UINT64_MAX, 0 };
    }
};

#endif //METABULI_SYNCMER_SCANNER_H