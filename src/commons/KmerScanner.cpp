#include "KmerScanner.h"
#include "MetamerPattern.h"
#include "GeneticCode.h"


MultiCodeScanner::MultiCodeScanner(const MultiCodePattern * pattern) 
    : KmerScanner(pattern->codePattern.size()), pattern(pattern) 
{
    for (size_t i = 0; i < pattern->geneticCodes.size(); ++i) {
        codonBitList.push_back(pattern->geneticCodes[i]->bitPerCodon);
        aaBitList.push_back(pattern->geneticCodes[i]->bitPerAA);
    }
    dnaPartList.resize(pattern->geneticCodes.size(), 0);
    aaPartList.resize(pattern->geneticCodes.size(), 0);
    for (size_t i = 0; i < pattern->codePattern.size(); ++i) {
        totalDNABits += codonBitList[pattern->codePattern[i]];
        totalAABits += aaBitList[pattern->codePattern[i]];
    }
    dnaMask = (1ULL << totalDNABits) - 1;
    aaMask = (1ULL << totalAABits) - 1;
}

void MultiCodeScanner::initScanner(const char * seq, size_t seqStart, size_t seqEnd, bool isForward) {
    KmerScanner::initScanner(seq, seqStart, seqEnd, isForward);
    this->aaLen = seqLen / 3;
    for (size_t i = 0; i < pattern->geneticCodes.size(); ++i) {
        dnaPartList[i] = 0;
        aaPartList[i] = 0;
    }
}

Kmer MultiCodeScanner::next() {
    while (posStart <= aaLen - kmerSize) {
        bool sawN = false;
        loadedCharCnt -= (loadedCharCnt == kmerSize);
        while (loadedCharCnt < kmerSize) {
            for (size_t c = 0; c < pattern->geneticCodes.size(); ++c) {
                int ci;
                int aa = 0;
                int codon = 0;
                if (isForward) {
                    ci = seqStart + (posStart + loadedCharCnt) * 3;
                    aa = pattern->geneticCodes[c]->getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                    codon = pattern->geneticCodes[c]->getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                } else {
                    ci = seqEnd - (posStart + loadedCharCnt) * 3;
                    aa = pattern->geneticCodes[c]->getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);          
                    codon = pattern->geneticCodes[c]->getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
                }
                if (aa < 0) { 
                    // std::cout << c << " " << posStart << " " << loadedCharCnt << " " << ci << std::endl;
                    sawN = true; 
                    break; 
                }
                dnaPartList[c] = (dnaPartList[c] << codonBitList[c]) | (uint64_t)codon;
                aaPartList[c] = (aaPartList[c] << aaBitList[c]) | (uint64_t)aa;
            }
            if (sawN) { break; }
            loadedCharCnt++;
        }
        if (sawN) {
            posStart += loadedCharCnt + 1;
            for (size_t i = 0; i < pattern->geneticCodes.size(); ++i) {
                dnaPartList[i] = 0;
                aaPartList[i] = 0;
            }
            loadedCharCnt = 0;
            continue;
        }
        uint64_t combinedAA = 0;
        uint64_t combinedDNA = 0;
        for (size_t i = 0; i < pattern->codePattern.size(); ++i) {
            const int codeIdx = pattern->codePattern[i];
            uint64_t aa = extract_bits(aaPartList[codeIdx], (kmerSize - 1 - i) * aaBitList[codeIdx], aaBitList[codeIdx]);
            uint64_t dna = extract_bits(dnaPartList[codeIdx], (kmerSize - 1 - i) * codonBitList[codeIdx], codonBitList[codeIdx]);
            combinedAA = (combinedAA << aaBitList[codeIdx]) | aa;
            combinedDNA = (combinedDNA << codonBitList[codeIdx]) | dna;
        }
        if (isForward) {
            return { ((combinedAA & aaMask) << totalDNABits) | (combinedDNA & dnaMask), seqStart + (posStart++) * 3 };
        } else {
            return { ((combinedAA & aaMask) << totalDNABits) | (combinedDNA & dnaMask), seqEnd - ((posStart++) + kmerSize) * 3 + 1 };
        }

    }
    return { UINT64_MAX, 0 }; // No more kmers found
}


Kmer MetamerScanner::next() 
{
    int aa = 0;
    int codon = 0;
    while (posStart <= aaLen - kmerSize) {
        bool sawN = false;
        loadedCharCnt -= (loadedCharCnt == kmerSize);
        while (loadedCharCnt < kmerSize) {
            int ci;
            if (isForward) {
                ci = seqStart + (posStart + loadedCharCnt) * 3;
                aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                codon = geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
            } else {
                ci = seqEnd - (posStart + loadedCharCnt) * 3;
                aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);          
                codon = geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
            }
            if (aa < 0) { sawN = true; break; }
            dnaPart = (dnaPart << bitsPerCodon) | (uint64_t)codon;
            aaPart = (aaPart << bitsPerAA) | (uint64_t)aa;
            loadedCharCnt++;
        }
        if (sawN) {
            posStart += loadedCharCnt + 1;
            dnaPart = aaPart = 0;
            loadedCharCnt = 0;
            continue;
        }
        if (isForward) {
            return { ((aaPart & aaMask) << dnaBits) | (dnaPart & dnaMask), seqStart + (posStart++) * 3 };
        } else {
            return { ((aaPart & aaMask) << dnaBits) | (dnaPart & dnaMask), seqEnd - ((posStart++) + kmerSize) * 3 + 1 };
        }
    }
    return { UINT64_MAX, 0 }; // No more kmers found
}


Kmer SpacedMetamerScanner::next() 
{
    int aa = 0;
    int codon = 0;
    while (posStart <= aaLen - windowSize) {
        bool sawN = false;
        loadedCharCnt -= (loadedCharCnt == windowSize);
        while (loadedCharCnt < windowSize) {
            int ci;
            if (isForward) {
                ci = seqStart + (posStart + loadedCharCnt) * 3;
                aa = geneticCode.getAA(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
                codon = geneticCode.getCodon(atcg[seq[ci]], atcg[seq[ci + 1]], atcg[seq[ci + 2]]);
            } else {
                ci = seqEnd - (posStart + loadedCharCnt) * 3;
                aa = geneticCode.getAA(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);          
                codon = geneticCode.getCodon(iRCT[atcg[seq[ci]]], iRCT[atcg[seq[ci - 1]]], iRCT[atcg[seq[ci - 2]]]);
            }
            if (aa < 0) { sawN = true; break; }
            dnaPart = (dnaPart << bitsPerCodon) | codon;
            aaPart = (aaPart << bitsPerAA) | aa;
            loadedCharCnt++;
        }
        if (sawN) {
            posStart += loadedCharCnt + 1;
            dnaPart = aaPart = 0;
            loadedCharCnt = 0;
            continue;
        }

        // const uint64_t currDnaPart = dnaPart & dnaMask;
        // const uint64_t currAaPart = aaPart & aaMask;

        uint64_t aaPartResult = 0;
        uint64_t dnaPartResult = 0;
        // int out = 0;
        // for (int i = 0; i < windowSize; ++i) {
        //     if (spaceMask & (1u << i)) {
        //         const uint64_t curAA    = (aaPart  >> (i * bitsPerAA))    & aaMaskPerCodon;
        //         const uint64_t curCodon = (dnaPart >> (i * bitsPerCodon)) & codonMask;
        //         aaPartResult  |= (curAA    <<  aaShifts[i]); //(out * bitsPerAA)
        //         dnaPartResult |= (curCodon << dnaShifts[i]);   // (out * bitsPerCodon)
        //         out++;
        //     }
        // }

        uint32_t validPosMask = spaceMask;
        while (validPosMask) {
            const int i = __builtin_ctz(validPosMask);
            validPosMask &= validPosMask - 1;
            const uint64_t myAA  = (aaPart  >> (i * bitsPerAA))    & aaMaskPerCodon;
            const uint64_t curCodon = (dnaPart >> (i * bitsPerCodon)) & codonMask;
            aaPartResult  |= (myAA  << aaShifts[i]);
            dnaPartResult |= (curCodon << dnaShifts[i]);   
        }

        if (isForward) {
            return { (aaPartResult << dnaBits) | (dnaPartResult), seqStart + (posStart++) * 3 };
        } else {
            return { (aaPartResult << dnaBits) | (dnaPartResult), seqEnd - ((posStart++) + windowSize) * 3 + 1 };
        }
    }
    return { UINT64_MAX, 0 }; // No more kmers found
}