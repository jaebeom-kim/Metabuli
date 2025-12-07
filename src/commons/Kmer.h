#ifndef ADKMER3_KMER_H
#define ADKMER3_KMER_H
#include <iostream>
#include "NcbiTaxonomy.h"
#include "GeneticCode.h"
#include <cstdint>
#include <bitset>

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

    virtual bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const = 0;
    virtual void printAA(uint64_t value) const = 0;
    virtual void printDNA(uint64_t value) const = 0;
    virtual std::pair<uint8_t, uint32_t> getHammingDists(uint64_t kmer1, uint64_t kmer2) const = 0;
    virtual std::pair<uint8_t, uint32_t> getHammingDists_reverse(uint64_t kmer1, uint64_t kmer2) const = 0;
    virtual std::vector<std::unique_ptr<KmerScanner>> createScanners(int num) const = 0;

};

class SingleCodePattern : public MetamerPattern {
public:
    const GeneticCode * geneticCode;
    int bitPerCodon;
    int bitPerAA;

    SingleCodePattern(const GeneticCode * geneticCode, int kmerLen) 
        : MetamerPattern(0, 0, 0, kmerLen), geneticCode(geneticCode), bitPerCodon(geneticCode->bitPerCodon), bitPerAA(geneticCode->bitPerAA)
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
        return (kmer1 & ((1U << (totalDNABits - bitPerCodon * shift)) - 1)) 
                == (kmer2 >> (bitPerCodon * shift));
    }

    void printAA(uint64_t value) const override {
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & ((1 << bitPerAA) - 1);
            std::cout << geneticCode->aminoacids[aa];
        }
    }
    
    void printDNA(uint64_t value) const override {
        uint64_t dnaPart = value & dnaMask;
        uint64_t aaPart = value >> totalDNABits;
        for (int i = kmerLen - 1; i >= 0; --i) {
            int aa = (aaPart >> (i * bitPerAA)) & ((1 << bitPerAA) - 1);
            int codon = (dnaPart >> (i * bitPerCodon)) & ((1 << bitPerCodon) - 1);
            std::cout << geneticCode->aa2codon[aa][codon];
        }
    }

    std::pair<uint8_t, uint32_t> getHammingDists(uint64_t kmer1, uint64_t kmer2) const override {
        uint8_t hammingSum = 0;
        uint32_t hammings = 0;

        const uint64_t aaPart = kmer1 >> totalDNABits;  
        const uint64_t dnaPart1 = kmer1 & dnaMask;
        const uint64_t dnaPart2 = kmer2 & dnaMask;

        for (size_t i = 0; i < kmerLen ; i++) {
            const int aa = (aaPart >> (totalAABits - (i + 1) * bitPerAA)) & ((1 << bitPerAA) - 1);
            const int codon1 = (dnaPart1 >> (totalDNABits - (i + 1) * bitPerCodon)) & ((1 << bitPerCodon) - 1);
            const int codon2 = (dnaPart2 >> (totalDNABits - (i + 1) * bitPerCodon)) & ((1 << bitPerCodon) - 1);
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

        for (size_t i = 0; i < kmerLen ; i++) {
            const int aa = (aaPart >> (i * bitPerAA)) & ((1 << bitPerAA) - 1);
            const int codon1 = (dnaPart1 >> (i * bitPerCodon)) & ((1 << bitPerCodon) - 1);
            const int codon2 = (dnaPart2 >> (i * bitPerCodon)) & ((1 << bitPerCodon) - 1);
            const int hammingDist = geneticCode->getHammingDist(aa, codon1, codon2);
            hammingSum += hammingDist;
            hammings |= hammingDist << (2 * i);
        }  
        return {hammingSum, hammings};
    }
};

class LegacyPattern : public SingleCodePattern {
public:
    LegacyPattern(const GeneticCode * geneticCode, int kmerLen) 
        : SingleCodePattern(geneticCode, kmerLen) {}
    
    bool checkOverlap(uint64_t kmer1, uint64_t kmer2, int shift) const override {
        return (kmer1 >> (bitPerCodon * shift)) == (kmer2 & ((1U << (totalDNABits - bitPerCodon * shift)) - 1));
    }
};

class MultiCodePattern : public MetamerPattern {
public:
    std::vector<const GeneticCode *> geneticCodes;
    std::vector<int> codePattern;
    std::vector<int> codonBitList;
    std::vector<int> aaBitList;

    MultiCodePattern() {}
    ~MultiCodePattern() {}
    MultiCodePattern(const std::vector<const GeneticCode *> &geneticCodes, const std::vector<int> & codePattern) 
        : MetamerPattern(0, 0, 0, codePattern.size()), geneticCodes(geneticCodes), codePattern(codePattern) 
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
        for (size_t i = 0; i < kmerLen ; i++) {
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
        for (size_t i = 0; i < kmerLen ; i++) {
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

struct QueryKmerInfo {
    QueryKmerInfo(uint32_t sequenceID, uint32_t pos, uint8_t frame)
        : sequenceID(sequenceID), pos(pos), frame(frame) {}
    QueryKmerInfo() = default;
    struct {
        uint64_t pos        : 32; // 32 bits
        uint64_t sequenceID : 29; // 29 bits
        uint64_t frame      : 3;  // 3 bits
    };
}; // 8 byte

struct TargetKmerInfo {
    explicit TargetKmerInfo(TaxID taxId = 0, TaxID speciesId = 0) : taxId(taxId), speciesId(speciesId) {}
    TaxID taxId;     // 4 byte
    TaxID speciesId; // 4 byte
};

struct Kmer {
    uint64_t value;
    union {
        uint32_t pos;
        uint32_t id;
        QueryKmerInfo qInfo;
        TargetKmerInfo tInfo;
    };

    Kmer() : value(0), id(0) {}

    Kmer(uint64_t value, TaxID taxid) : value(value), id(uint32_t(taxid)){}

    Kmer(uint64_t value, uint32_t id) : value(value), id(id) {}

    Kmer(uint64_t value, const QueryKmerInfo & qInfo) : value(value), qInfo(qInfo) {}

    Kmer(uint64_t value, const TargetKmerInfo & tInfo) : value(value), tInfo(tInfo) {}

    Kmer(uint64_t value, TaxID taxId, TaxID speciesId) : value(value), tInfo(taxId, speciesId) {}

    Kmer(uint64_t value, uint32_t seqId, uint32_t pos, uint8_t frame) 
        : value(value), qInfo(seqId, pos, frame) {}

    bool isEmpty() const {
        return value == 0 && id == 0;
    }

    void printAA(const GeneticCode * code) const {
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            std::cout << code->aminoacids[aa];
        }
    }

    void printAA(const GeneticCode * code, int k) const {
        for (int i = 0; i < k; ++i) {
            int aa = (value >> (((k - 1) * 5) - 5 * i)) & 0x1F;
            std::cout << code->aminoacids[aa];
        }
    }
        
    // void printAA(const MetamerPattern * pattern) const {
    //     uint64_t aaPart = value >> (pattern->totalDNABits);
    //     int k = pattern->kmerLen;
    //     int aaBitSum = 0;
    //     for (int i = 0; i < k; ++i) {
    //         const int * currAABit = &pattern->aaBitList[pattern->codePattern[i]];
    //         aaBitSum += *currAABit;
    //         int aa = (aaPart >> (pattern->totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
    //         std::cout << pattern->geneticCodes[pattern->codePattern[i]]->aminoacids[aa];
    //     }
    // }

    // void printDNA(const MetamerPattern * pattern) const {
    //     uint64_t dnaPart = value & pattern->dnaMask;
    //     uint64_t aaPart = value >> (pattern->totalDNABits);
    //     int k = pattern->codePattern.size();
    //     int aaBitSum = 0;
    //     int dnaBitSum = 0;
    //     for (int i = 0; i < k; ++i) {
    //         const int * currAABit = &pattern->aaBitList[pattern->codePattern[i]];
    //         aaBitSum += *currAABit;
    //         int aa = (aaPart >> (pattern->totalAABits - aaBitSum)) & ((1 << *currAABit) - 1);
            
    //         const int * currCodonBit = &pattern->codonBitList[pattern->codePattern[i]];
    //         dnaBitSum += *currCodonBit;
    //         int codon = (dnaPart >> (pattern->totalDNABits - dnaBitSum)) & ((1 << *currCodonBit) - 1);

    //         std::cout << pattern->geneticCodes[pattern->codePattern[i]]->aa2codon[aa][codon];
    //     }
    // }

    void printDNA(const GeneticCode * code) const {
        uint64_t dnaPart = value & 0xFFFFFF;
        uint64_t aaPart = value >> 24;
        for (int i = 0; i < 8; ++i) {
            int aa = (aaPart >> (35 - 5 * i)) & 0x1F;
            int codon = (dnaPart >> (21 - 3 * i)) & 0x7;
            std::cout << code->aa2codon[aa][codon];
        }
    }

    static bool compareTargetKmer(const Kmer & a, const Kmer & b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }

        if (a.tInfo.speciesId != b.tInfo.speciesId) {
            return a.tInfo.speciesId < b.tInfo.speciesId;
        }

        return a.tInfo.taxId < b.tInfo.taxId;
    }

    static bool compareQueryKmer(const Kmer &a, const Kmer &b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }

        if (a.qInfo.sequenceID != b.qInfo.sequenceID) {
            return a.qInfo.sequenceID < b.qInfo.sequenceID;
        }

        if (a.qInfo.frame != b.qInfo.frame) {
            return a.qInfo.frame < b.qInfo.frame;
        }
        
        return a.qInfo.pos < b.qInfo.pos;
    }

    static bool compareQKmerByIdAndPos(const Kmer &a, const Kmer &b) {
        if (a.qInfo.sequenceID != b.qInfo.sequenceID) {
            return a.qInfo.sequenceID < b.qInfo.sequenceID;
        }
        return a.qInfo.pos < b.qInfo.pos;
    }

    static bool compareKmer(const Kmer &a, const Kmer &b) {
        if (a.value != b.value) {
            return a.value < b.value;
        }
        return a.id < b.id;
    }
};


struct DiffIdxSplit{
    DiffIdxSplit(uint64_t ADkmer, size_t diffIdxOffset, size_t infoIdxOffset) : ADkmer(ADkmer), diffIdxOffset(diffIdxOffset), infoIdxOffset(infoIdxOffset) { }
    DiffIdxSplit(const DiffIdxSplit & copy) {ADkmer = copy.ADkmer; diffIdxOffset = copy.diffIdxOffset; infoIdxOffset=copy.infoIdxOffset;}
    DiffIdxSplit() {};
    DiffIdxSplit& operator=(const DiffIdxSplit&) = default;
    uint64_t ADkmer;
    size_t diffIdxOffset;
    size_t infoIdxOffset;
};

struct Metamer {
    Metamer() : metamer(0), id(0) {}
    Metamer(uint64_t metamer, uint32_t id) : metamer(metamer), id(id) {}
    uint64_t metamer;
    uint32_t id; // it is mapped to taxonomy ID and protein ID

    static std::bitset<96> substract(const Metamer & metamer1, const Metamer & metamer2) {
        // metamer 1 is the same or greater than metamer2
        if (metamer1.metamer == metamer2.metamer) {
            return std::bitset<96>(metamer1.id - metamer2.id);
        }
        if (metamer1.id >= metamer2.id) {
            std::bitset<96> result;
            result = metamer1.metamer - metamer2.metamer;
            result <<= 30;
            result |= (metamer1.id - metamer2.id);
            return result;
        }
        std::bitset<96> result;
        uint64_t diff = metamer1.metamer - metamer2.metamer - 1;
        result = diff;
        result <<= 30;
        result |= (((1U << 30) - 1) - metamer2.id + metamer1.id + 1);
        return result;    
    }


    Metamer add(const std::bitset<96> & diff) const {
        uint64_t idSum = this->id + (diff & std::bitset<96>(0x3FFFFFFF)).to_ullong();
        uint64_t metamerSum = this->metamer + (diff >> 30).to_ullong() + (idSum >> 30);
        idSum &= 0x3FFFFFFF;
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


#endif //ADKMER3_KMER_H
