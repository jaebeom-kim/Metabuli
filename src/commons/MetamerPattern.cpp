#include "MetamerPattern.h"
#include "json.h"

using json = nlohmann::json;

float SingleCodePattern::substitutionScore_R(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix & matrix) const 
{
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = 0; i < count; ++i) {
        const int myAA = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        const int idx1 = static_cast<int>(matrix.aa2num[translateNucl->translateSingleCodon(geneticCode->aa2codon[myAA][codon1].c_str())]); // ToDo: precompute it
        const int idx2 = static_cast<int>(matrix.aa2num[translateNucl->translateSingleCodon(geneticCode->aa2codon[myAA][codon2].c_str())]);
        score += matrix.subMatrix[idx1][idx2];
        aaBitSum += bitPerAA;
        dnaBitSum += bitPerCodon;
    }  
    return score;
}

float SingleCodePattern::substitutionScore_L(
    uint64_t kmer1,
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix & matrix) const
{
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA;
        dnaBitSum -= bitPerCodon;
        const int myAA   = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        const int idx1 = matrix.aa2num[
            translateNucl->translateSingleCodon(geneticCode->aa2codon[myAA][codon1].c_str())
        ];
        const int idx2 = matrix.aa2num[
            translateNucl->translateSingleCodon(geneticCode->aa2codon[myAA][codon2].c_str())
        ]; // ToDo: precompute it

        score += matrix.subMatrix[idx1][idx2];
    }

    return score;
}

float SingleCodePattern::hammingDistScore_R(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count) const 
{
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = 0; i < count; ++i) {
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        const int hammingDist = geneticCode->getHammingDist(aa, codon1, codon2);
        if (hammingDist == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * hammingDist;
        }
        aaBitSum += bitPerAA;
        dnaBitSum += bitPerCodon; 
    }  
    return score;
}

float SingleCodePattern::hammingDistScore_L(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count) const 
{
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA;
        dnaBitSum -= bitPerCodon;
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        const int hammingDist = geneticCode->getHammingDist(aa, codon1, codon2);
        if (hammingDist == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * hammingDist;
        }
    }  
    return score;
}

uint8_t SingleCodePattern::hammingDistSum_R(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int count) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = 0; i < count; ++i) {
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        hammingSum += geneticCode->getHammingDist(aa, codon1, codon2);
        aaBitSum += bitPerAA;
        dnaBitSum += bitPerCodon; 
    }  
    return hammingSum;
}

uint8_t SingleCodePattern::hammingDistSum_L(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int count) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA;
        dnaBitSum -= bitPerCodon;
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        hammingSum += geneticCode->getHammingDist(aa, codon1, codon2);
    }  
    return hammingSum;
}

MultiCodePattern::MultiCodePattern(const std::string & customFile) {
    std::ifstream file(customFile);
    if (!file.is_open()) {
        std::cout << "Error: Could not open custom metamer pattern file: " << customFile << std::endl;
    }

    std::string line;
    std::stringstream jsonbuf;
    bool inJson = false;
    while (std::getline(file, line)) {
        if (line == "===BEGIN_CUSTOM_METAMER===") {
            inJson = true;
            continue;   
        }
        if (line == "===END_CUSTOM_METAMER===") {
            inJson = false;
            break;  
        }
        if (inJson) {
            jsonbuf << line << "\n";
        }
    }

    json j = json::parse(jsonbuf.str());

    std::string name = j["name"];
    kmerLen = j["length"];
    codePattern = j["position_codes"].get<std::vector<int>>();

    if (kmerLen != static_cast<int>(codePattern.size())) {
        std::cout << "Error: Length of code pattern does not match the specified length." << std::endl;
    }


    for (const auto& codeEntry : j["codes"]) {

        const auto& codons_json = codeEntry["codons"];

        std::cout << "Loading custom genetic code: " << codeEntry["id"] << std::endl;

        std::vector<std::pair<std::string, std::vector<std::string>>> translationTable;

        for (const auto& entry : codons_json) {
            std::string aaStr = entry[0].get<std::string>();
            std::vector<std::string> codons = entry[1].get<std::vector<std::string>>();
            translationTable.emplace_back(std::string{aaStr[0]}, codons);
        }

        geneticCodes.push_back(std::make_unique<CustomGeneticCode>(translationTable));
    }

    // Validate total bits do not exceed 64 bits
    int totalBits = 0;
    for (size_t i = 0; i < codePattern.size(); ++i) {
        totalBits += geneticCodes[codePattern[i]]->bitPerAA;
        totalBits += geneticCodes[codePattern[i]]->bitPerCodon;
        if (totalBits > 64) {
            std::cout << "Error: Total bits exceed 64 bits." << std::endl;
        }
    }

    init();
}


bool MultiCodePattern::checkOverlap(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int shift) const 
{
    uint64_t aaPart1 = kmer1 >> totalDNABits;
    uint64_t aaPart2 = kmer2 >> totalDNABits;
    int range = kmerLen - shift;
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
        int codon1 = (kmer1 >> (totalDNABits - dnaBitSum1)) & ((1 << currCodonBit) - 1);
        std::string dna = geneticCodes[codePattern[i + shift]]->aa2codon[aa1][codon1];

        // Codon of k-mer 2
        currCodonBit = codonBitList[codePattern[i]];
        dnaBitSum2 += currCodonBit;
        int codon2 = (kmer2 >> (totalDNABits - dnaBitSum2)) & ((1 << currCodonBit) - 1);
        std::string dna2 = geneticCodes[codePattern[i]]->aa2codon[aa2][codon2];

        if (dna != dna2) {
            return false;
        }
    }
    return true;
}

float MultiCodePattern::substitutionScore_R(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix & matrix) const 
{
    float totalScore = 0.0f;
    uint64_t aaPart = kmer1 >> totalDNABits; 

    // Get substitution score for count times from the right end.
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = kmerLen - 1; i >= kmerLen - count; --i) {
        const int currAABit = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];
        const int myAA = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int idx1 = static_cast<int>(matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCodes[codePattern[i]]->aa2codon[myAA][codon1].c_str()
            )
        ]);
        const int idx2 = static_cast<int>(matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCodes[codePattern[i]]->aa2codon[myAA][codon2].c_str()
            )
        ]); // To Do: Get A.A. using a canonical translation table before doing substitution

        // std::cout << "myAA: " << geneticCodes[codePattern[i]]->aminoacids[myAA] << "(" << myAA << ")"  
        //     << " codon1: " << geneticCodes[codePattern[i]]->aa2codon[myAA][codon1] << "(" << codon1 << ")" 
        //     << " codon2: " << geneticCodes[codePattern[i]]->aa2codon[myAA][codon2] << "(" << codon2 << ")" << std::endl;
        // std::cout << "aa1: " << translateNucl->translateSingleCodon(
        //         geneticCodes[codePattern[i]]->aa2codon[myAA][codon1].c_str()
        //     ) << " idx1: " << idx1 << std::endl;
        // std::cout << "aa2: " << translateNucl->translateSingleCodon(
        //         geneticCodes[codePattern[i]]->aa2codon[myAA][codon2].c_str()
        //     ) << " idx2: " << idx2 << std::endl;
        // std::cout << "Score: " << matrix.subMatrix[idx1][idx2] << std::endl;

        totalScore += matrix.subMatrix[idx1][idx2]; // To Do: Get A.A. using a canonical translation table before doing substitution score lookup.
        aaBitSum += currAABit;
        dnaBitSum += currDnaBit;
    }
    return totalScore;
}

float MultiCodePattern::substitutionScore_L(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix & matrix) const 
{
    float totalScore = 0.0f;
    uint64_t aaPart = kmer1 >> totalDNABits; 

    // Get substitution score for count times from the right end.
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < count; ++i) {
        const int currAABit = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];
        aaBitSum -= currAABit;
        dnaBitSum -= currDnaBit;
        const int myAA   = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int idx1 = static_cast<int>(matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCodes[codePattern[i]]->aa2codon[myAA][codon1].c_str()
            )
        ]);
        const int idx2 = static_cast<int>(matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCodes[codePattern[i]]->aa2codon[myAA][codon2].c_str()
            )
        ]); // To Do: Get A.A. using a canonical translation table before doing substitution

        // std::cout << "myAA: " << geneticCodes[codePattern[i]]->aminoacids[myAA] << "(" << myAA << ")"  
        //     << " codon1: " << geneticCodes[codePattern[i]]->aa2codon[myAA][codon1] << "(" << codon1 << ")" 
        //     << " codon2: " << geneticCodes[codePattern[i]]->aa2codon[myAA][codon2] << "(" << codon2 << ")" << std::endl;
        // std::cout << "aa1: " << translateNucl->translateSingleCodon(
        //         geneticCodes[codePattern[i]]->aa2codon[myAA][codon1].c_str()
        //     ) << " idx1: " << idx1 << std::endl;
        // std::cout << "aa2: " << translateNucl->translateSingleCodon(
        //         geneticCodes[codePattern[i]]->aa2codon[myAA][codon2].c_str()
        //     ) << " idx2: " << idx2 << std::endl;
        // std::cout << "Score: " << matrix.subMatrix[idx1][idx2] << std::endl;
        
        totalScore += matrix.subMatrix[idx1][idx2]; // To Do: Get A.A. using a canonical translation table before doing substitution score lookup.
    }
    return totalScore;
}



float MultiCodePattern::hammingDistScore_R(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count) const 
{
    float totalScore = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = kmerLen - 1; i >= kmerLen - count; --i) {
        const int currAABit = aaBitList[codePattern[i]];
        const int currCodonBit = codonBitList[codePattern[i]];

        const int aa = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currCodonBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currCodonBit) - 1);

        const int hammingDist = geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
        if (hammingDist == 0) {
            totalScore += 3.0f;
        } else {
            totalScore += 2.0f - 0.5f * hammingDist;
        }
        aaBitSum += currAABit;
        dnaBitSum += currCodonBit; 
    }
   
    return totalScore;
}

float MultiCodePattern::hammingDistScore_L(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count) const 
{
    float totalScore = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;

    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < count; ++i) {
        const int currAABit = aaBitList[codePattern[i]];
        const int currCodonBit = codonBitList[codePattern[i]];
        aaBitSum -= currAABit;
        dnaBitSum -= currCodonBit;
        
        const int aa = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currCodonBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currCodonBit) - 1);
        
        const int hammingDist = geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
        if (hammingDist == 0) {
            totalScore += 3.0f;
        } else {
            totalScore += 2.0f - 0.5f * hammingDist;
        }
    }
   
    return totalScore;
}

uint8_t MultiCodePattern::hammingDistSum_R(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int count) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = kmerLen - 1; i >= kmerLen - count; --i) {
        const int aa = (aaPart >> aaBitSum) & ((1 << aaBitList[codePattern[i]]) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << codonBitList[codePattern[i]]) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << codonBitList[codePattern[i]]) - 1);
        hammingSum += geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
        aaBitSum += aaBitList[codePattern[i]];
        dnaBitSum += codonBitList[codePattern[i]]   ; 
    }  
    return hammingSum;
}

uint8_t MultiCodePattern::hammingDistSum_L(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int count) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < count; ++i) {
        aaBitSum -= aaBitList[codePattern[i]];
        dnaBitSum -= codonBitList[codePattern[i]];
        const int aa = (aaPart >> aaBitSum) & ((1 << aaBitList[codePattern[i]]) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << codonBitList[codePattern[i]]) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << codonBitList[codePattern[i]]) - 1);
        hammingSum += geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
    }  
    return hammingSum;
}