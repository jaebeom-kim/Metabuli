#ifndef METAMER_PATTERN_H
#define METAMER_PATTERN_H

#include <iostream>
#include <vector>
#include <memory>
#include "KmerScanner.h"
#include "GeneticCode.h"

class MetamerPattern {
public:
    uint64_t dnaMask; // Value & dnaMask -> DNA part
    int totalDNABits;
    int totalAABits;
    int kmerLen;

    MetamerPattern() : dnaMask(0), totalDNABits(0), totalAABits(0), kmerLen(0) {}
    MetamerPattern(uint64_t dnaMask, int totalDNABits, int totalAABits, int kmerLen) 
        : dnaMask(dnaMask), totalDNABits(totalDNABits), totalAABits(totalAABits), kmerLen(kmerLen) {}
    
    virtual ~MetamerPattern() {}

    virtual std::vector<std::unique_ptr<KmerScanner>> createScanners(int num) const = 0;
    virtual bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const = 0;
    virtual void printAA(uint64_t value) const = 0;
    virtual void printDNA(uint64_t value) const = 0;
    virtual std::pair<uint8_t, uint32_t> getHammingDists(uint64_t kmer1, uint64_t kmer2) const = 0;
    virtual std::pair<uint8_t, uint32_t> getHammingDists_reverse(uint64_t kmer1, uint64_t kmer2) const = 0;
    
};

class SingleCodePattern : public MetamerPattern {
public:
    std::unique_ptr<GeneticCode> geneticCode;
    int bitPerCodon;
    int bitPerAA;
    const uint64_t codonMask;
    const uint64_t aaMask;
    

    SingleCodePattern(std::unique_ptr<GeneticCode> code, int kmerLen) 
        : MetamerPattern(0, 0, 0, kmerLen),
         geneticCode(std::move(code)),
         bitPerCodon(geneticCode->bitPerCodon), 
         bitPerAA(geneticCode->bitPerAA),
         codonMask((1ULL << bitPerCodon) - 1),
         aaMask((1ULL << bitPerAA) - 1)
    {
        totalDNABits = kmerLen * geneticCode->bitPerCodon;
        totalAABits = kmerLen * geneticCode->bitPerAA;
        dnaMask = (1ULL << totalDNABits) - 1;
    }

    std::vector<std::unique_ptr<KmerScanner>> createScanners(int num) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < num; ++i) {
            scanners.push_back(std::make_unique<MetamerScanner>(*geneticCode, kmerLen));
        }
        return scanners;
    }

    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;
        return (dnaPart1 & ((1U << (totalDNABits - bitPerCodon * shift)) - 1)) 
                == (dnaPart2 >> (bitPerCodon * shift));
    }

    void printAA(uint64_t value) const override {
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            std::cout << geneticCode->aminoacids[aa];
        }
    }
    
    void printDNA(uint64_t value) const override {
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            int codon = (dnaPart >> (i * bitPerCodon)) & codonMask;
            std::cout << geneticCode->aa2codon[aa][codon];
        }
    }

    std::pair<uint8_t, uint32_t> getHammingDists(uint64_t kmer1, uint64_t kmer2) const override {
        uint8_t hammingSum = 0;
        uint32_t hammings = 0;

        const uint64_t aaPart = kmer1 >> totalDNABits;  
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;

        for (int i = 0; i < kmerLen ; i++) {
            const int aa = (aaPart >> (totalAABits - (i + 1) * bitPerAA)) & aaMask;
            const int codon1 = (dnaPart1 >> (totalDNABits - (i + 1) * bitPerCodon)) & codonMask;
            const int codon2 = (dnaPart2 >> (totalDNABits - (i + 1) * bitPerCodon)) & codonMask;
            const int hammingDist = geneticCode->getHammingDist(aa, codon1, codon2);
            hammingSum += hammingDist;
            hammings |= hammingDist << (2 * (kmerLen - 1 -i));
        }  
        return {hammingSum, hammings};
    }

    std::pair<uint8_t, uint32_t> getHammingDists_reverse(uint64_t kmer1, uint64_t kmer2) const override {
        uint8_t hammingSum = 0;
        uint32_t hammings = 0;

        const uint64_t aaPart = kmer1 >> totalDNABits;  
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;

        for (int i = 0; i < kmerLen ; i++) {
            const int aa = (aaPart >> (totalAABits - (i + 1) * bitPerAA)) & aaMask;
            const int codon1 = (dnaPart1 >> (totalDNABits - (i + 1) * bitPerCodon)) & codonMask;
            const int codon2 = (dnaPart2 >> (totalDNABits - (i + 1) * bitPerCodon)) & codonMask;
            const int hammingDist = geneticCode->getHammingDist(aa, codon1, codon2);
            hammingSum += hammingDist;
            hammings |= hammingDist << (2 * i);
        }  
        return {hammingSum, hammings};
    }
};

class LegacyPattern : public SingleCodePattern {
public:
    LegacyPattern(std::unique_ptr<GeneticCode> geneticCode, int kmerLen) 
        : SingleCodePattern(std::move(geneticCode), kmerLen) {}

    std::vector<std::unique_ptr<KmerScanner>> createScanners(int num) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < num; ++i) {
            scanners.push_back(std::make_unique<OldMetamerScanner>(*geneticCode));
        }
        return scanners;
    }

    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;
        return (dnaPart1 >> (bitPerCodon * shift)) == (dnaPart2 & ((1U << (totalDNABits - bitPerCodon * shift)) - 1));
    }
};

class MultiCodePattern : public MetamerPattern {
public:
    std::vector<std::unique_ptr<GeneticCode>> geneticCodes;
    std::vector<int> codePattern;
    std::vector<int> codonBitList;
    std::vector<int> aaBitList;

    MultiCodePattern() {}
    ~MultiCodePattern() {}

    MultiCodePattern(const std::string & customFile);

    std::vector<std::unique_ptr<KmerScanner>> createScanners(int num) const override {
        std::vector<std::unique_ptr<KmerScanner>> scanners;
        for (int i = 0; i < num; ++i) {
            scanners.push_back(std::make_unique<MultiCodeScanner>(this));
        }
        return scanners;
    }

    MultiCodePattern(std::vector<std::unique_ptr<GeneticCode>> &&geneticCodes, const std::vector<int> & codePattern) 
        : MetamerPattern(0, 0, 0, codePattern.size()), 
          geneticCodes(std::move(geneticCodes)), 
          codePattern(codePattern) 
    {
        for (size_t i = 0; i < geneticCodes.size(); ++i) {
            codonBitList.push_back(geneticCodes[i]->bitPerCodon);
            aaBitList.push_back(geneticCodes[i]->bitPerAA);
        }
        for (size_t i = 0; i < codePattern.size(); ++i) {
            totalDNABits += codonBitList[codePattern[i]];
            totalAABits += aaBitList[codePattern[i]];
        }
        dnaMask = (1ULL << totalDNABits) - 1;
    }

    void init() {
        for (size_t i = 0; i < geneticCodes.size(); ++i) {
            codonBitList.push_back(geneticCodes[i]->bitPerCodon);
            aaBitList.push_back(geneticCodes[i]->bitPerAA);
        }
        for (size_t i = 0; i < codePattern.size(); ++i) {
            totalDNABits += codonBitList[codePattern[i]];
            totalAABits += aaBitList[codePattern[i]];
        }
        std::cout << "Total DNA bits: " << totalDNABits << ", Total AA bits: " << totalAABits << std::endl;
        dnaMask = (1ULL << totalDNABits) - 1;
    }

    void printAA(uint64_t value) const override {
        uint64_t aaPart = value >> totalDNABits;
        int k = codePattern.size();
        int aaBitSum = 0;
        for (int i = 0; i < k; ++i) {
            const int * currAABit = &aaBitList[codePattern[i]];
            aaBitSum += *currAABit;
            int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
            std::cout << geneticCodes[codePattern[i]]->aminoacids[aa];
        }
    }

    void printDNA(uint64_t value) const override {
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        int k = codePattern.size();
        int aaBitSum = 0;
        int dnaBitSum = 0;
        for (int i = 0; i < k; ++i) {
            const int * currAABit = &aaBitList[codePattern[i]];
            aaBitSum += *currAABit;
            int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
            
            const int * currCodonBit = &codonBitList[codePattern[i]];
            dnaBitSum += *currCodonBit;
            int codon = (dnaPart >> (totalDNABits - dnaBitSum)) & ((1 << *currCodonBit) - 1);

            std::cout << geneticCodes[codePattern[i]]->aa2codon[aa][codon];
        }
    }



    // It checks if rear part of kmer1 and front part of kmer2 overlap given a shift
    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        uint64_t dnaPart1 = kmer1 & dnaMask;
        uint64_t dnaPart2 = kmer2 & dnaMask;
        uint64_t aaPart1 = kmer1 >> totalDNABits;
        uint64_t aaPart2 = kmer2 >> totalDNABits;
        int k = codePattern.size();
        int range = k - shift;
        int aaBitSum1 = 0;
        int aaBitSum2 = 0;
        int dnaBitSum1 = 0;
        int dnaBitSum2 = 0;
        for (int i = 0; i < shift; ++i) {
            aaBitSum1 += aaBitList[codePattern[i]];
            dnaBitSum1 += codonBitList[codePattern[i]];
        }
        for (int i = 0; i < range; ++i) {
            // A.A. of k-mer 1
            int currAABit = aaBitList[codePattern[i + shift]];
            aaBitSum1 += currAABit;
            int aa1 = (aaPart1 >> (totalAABits - aaBitSum1)) & ((1 << currAABit) - 1);

            // A.A. of k-mer 2
            currAABit = aaBitList[codePattern[i]];
            aaBitSum2 += currAABit;
            int aa2 = (aaPart2 >> (totalAABits - aaBitSum2)) & ((1 << currAABit) - 1);

            // Codon of k-mer 1
            int currCodonBit = codonBitList[codePattern[i + shift]];
            dnaBitSum1 += currCodonBit;
            int codon1 = (dnaPart1 >> (totalDNABits - dnaBitSum1)) & ((1 << currCodonBit) - 1);
            std::string dna = geneticCodes[codePattern[i + shift]]->aa2codon[aa1][codon1];

            // Codon of k-mer 2
            currCodonBit = codonBitList[codePattern[i]];
            dnaBitSum2 += currCodonBit;
            int codon2 = (dnaPart2 >> (totalDNABits - dnaBitSum2)) & ((1 << currCodonBit) - 1);
            std::string dna2 = geneticCodes[codePattern[i]]->aa2codon[aa2][codon2];

            if (dna != dna2) {
                return false;
            }
        }
        return true;
    }

    inline std::pair<uint8_t, uint32_t> getHammingDists(uint64_t kmer1, uint64_t kmer2) const override {
        uint8_t hammingSum = 0;
        uint32_t hammings = 0;

        const uint64_t aaPart = kmer1 >> totalDNABits;  
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;

        int aaBitSum = 0;
        int dnaBitSum = 0;
        for (int i = 0; i < kmerLen ; i++) {
            const uint8_t patternIdx = codePattern[i];
            const int currAABit = aaBitList[patternIdx];
            aaBitSum += currAABit;
            const int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << currAABit) - 1);
            const int currCodonBit = codonBitList[patternIdx];
            dnaBitSum += currCodonBit;
            const uint64_t codonMask = (1ULL << currCodonBit) - 1;
            const int codon1 = (dnaPart1 >> (totalDNABits - dnaBitSum)) & codonMask;
            const int codon2 = (dnaPart2 >> (totalDNABits - dnaBitSum)) & codonMask;
            const int hammingDist = geneticCodes[patternIdx]->getHammingDist(aa, codon1, codon2);
            hammingSum += hammingDist;
            hammings |= hammingDist << (2 * (kmerLen - 1 -i));
        } 

      return {hammingSum, hammings};
    }

    inline std::pair<uint8_t, uint32_t> getHammingDists_reverse(uint64_t kmer1, uint64_t kmer2) const override {
        uint8_t hammingSum = 0;
        uint32_t hammings = 0;
        const uint64_t aaPart = kmer1 >> totalDNABits;  
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;

        int aaBitSum = 0;
        int dnaBitSum = 0;
        for (int i = 0; i < kmerLen ; i++) {
          const uint8_t patternIdx = codePattern[i];
          const int currAABit = aaBitList[patternIdx];
          aaBitSum += currAABit;
          const int aa = (aaPart >> (totalAABits - aaBitSum)) & ((1 << currAABit) - 1);
        
          const int currCodonBit = codonBitList[patternIdx];
          dnaBitSum += currCodonBit;
          const uint64_t codonMask = (1ULL << currCodonBit) - 1;
          const int codon1 = (dnaPart1 >> (totalDNABits - dnaBitSum)) & codonMask;
          const int codon2 = (dnaPart2 >> (totalDNABits - dnaBitSum)) & codonMask;
          const int hammingDist = geneticCodes[patternIdx]->getHammingDist(aa, codon1, codon2);
          hammingSum += hammingDist;
          hammings |= hammingDist << (2 * i);
        } 
    
        return {hammingSum, hammings};
    }

};

#endif // METAMER_PATTERN_H