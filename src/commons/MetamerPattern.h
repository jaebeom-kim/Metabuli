#ifndef METAMER_PATTERN_H
#define METAMER_PATTERN_H

#include <iostream>
#include <vector>
#include <memory>
#include "KmerScanner.h"
#include "GeneticCode.h"
#include "SubstitutionMatrix.h"
#include "TranslateNucl.h"

int getCodeNum(const std::string & customFile);

class MetamerPattern {
public:
    uint64_t dnaMask; // Value & dnaMask -> DNA part
    int totalDNABits;
    int totalAABits;
    int kmerLen;
    TranslateNucl * translateNucl = nullptr;

    MetamerPattern() : dnaMask(0), totalDNABits(0), totalAABits(0), kmerLen(0) {
        translateNucl = new TranslateNucl(TranslateNucl::CANONICAL);
    }
    MetamerPattern(uint64_t dnaMask, int totalDNABits, int totalAABits, int kmerLen) 
        : dnaMask(dnaMask), totalDNABits(totalDNABits), totalAABits(totalAABits), kmerLen(kmerLen) {
        translateNucl = new TranslateNucl(TranslateNucl::CANONICAL);
    }
    
    virtual ~MetamerPattern() {
        delete translateNucl;
    }

    virtual std::vector<std::unique_ptr<KmerScanner>> createScanners(int num) const = 0;    


    virtual float substitutionScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const = 0;
    virtual float hammingDistScore(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const = 0;
    virtual uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const = 0;
    virtual uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const  = 0;
    virtual MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const = 0;
    
    virtual bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const = 0;
    
    virtual void printAA(uint64_t value) const = 0;
    virtual void printDNA(uint64_t value) const = 0;
    virtual std::string toDnaString(uint64_t value) const = 0;
    
};

class SingleCodePattern : public MetamerPattern {
public:
    std::unique_ptr<GeneticCode> geneticCode;
    int bitPerCodon;
    int bitPerAA;
    uint64_t codonMask;
    uint64_t aaMask;
    
    SingleCodePattern(const std::string & customFile);
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

    float substitutionScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const override;
    float hammingDistScore(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const override;
    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const override;

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

    std::string toDnaString(uint64_t value) const override {
        std::string dnaStr;
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & aaMask;
            int codon = (dnaPart >> (i * bitPerCodon)) & codonMask;
            dnaStr += geneticCode->aa2codon[aa][codon];
        }
        return dnaStr;
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

    MatchScore calMatchScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const override;
    float substitutionScore(uint64_t kmer1, uint64_t kmer2, int count, const SubstitutionMatrix& matrix, bool fromR) const override;
    float hammingDistScore(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const override;
    uint8_t hammingDistSum(uint64_t kmer1, uint64_t kmer2) const override;

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

    std::string toDnaString(uint64_t value) const override {
        std::string dnaStr;
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

            dnaStr += geneticCodes[codePattern[i]]->aa2codon[aa][codon];
        }
        return dnaStr;
    }

    // It checks if rear part of kmer1 and front part of kmer2 overlap given a shift
    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override;

};

#endif // METAMER_PATTERN_H