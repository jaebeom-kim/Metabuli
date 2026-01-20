#include "MetamerPattern.h"
#include "json.h"

using json = nlohmann::json;

int getCodeNum(const std::string & customFile) {
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
    int code_count = j["code_count"];
    file.close();
    return code_count;
}

SingleCodePattern::SingleCodePattern(const std::string & customFile) {
    std::cout << "Loading single-code metamer pattern from file: " << customFile << std::endl;
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

    for (const auto& codeEntry : j["codes"]) {

        const auto& codons_json = codeEntry["codons"];

        std::cout << "Loading custom genetic code: " << codeEntry["id"] << std::endl;

        std::vector<std::pair<std::string, std::vector<std::string>>> translationTable;

        for (const auto& entry : codons_json) {
            std::string aaStr = entry[0].get<std::string>();
            std::vector<std::string> codons = entry[1].get<std::vector<std::string>>();
            translationTable.emplace_back(std::string{aaStr[0]}, codons);
        }
        geneticCode = std::make_unique<CustomGeneticCode>(translationTable);
    }

    // Validate total bits do not exceed 64 bits
    int totalBits = 0;
    totalBits += geneticCode->bitPerAA;
    totalBits += geneticCode->bitPerCodon;
    if (totalBits > 64) {
        std::cout << "Error: Total bits exceed 64 bits." << std::endl;
    }

    bitPerCodon = geneticCode->bitPerCodon;
    bitPerAA = geneticCode->bitPerAA;
    totalDNABits = kmerLen * bitPerCodon;
    totalAABits = kmerLen * bitPerAA;
    dnaMask = (1ULL << totalDNABits) - 1;
    codonMask = (1ULL << bitPerCodon) - 1;
    aaMask = (1ULL << bitPerAA) - 1;
    file.close();
}

MatchScore SingleCodePattern::calMatchScore(
    uint64_t kmer1,
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix& matrix,
    bool fromR) const 
{
    float idScore = 0.0f;
    float subScore = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;

    int aaBitSum  = (fromR) ? 0 : totalAABits;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;

    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA * !fromR;
        dnaBitSum -= bitPerCodon * !fromR;

        const int myAA   = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;

        const int idx1 = matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCode->aa2codon[myAA][codon1].c_str())
        ];
        const int idx2 = matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCode->aa2codon[myAA][codon2].c_str())
        ];

        const int hammingDist = geneticCode->getHammingDist(myAA, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        subScore += matrix.subMatrix[idx1][idx2];

        aaBitSum += bitPerAA * fromR;
        dnaBitSum += bitPerCodon * fromR;
    }
    return {idScore, subScore};
}
    


float SingleCodePattern::substitutionScore(
    uint64_t kmer1,
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix& matrix,
    bool fromR
) const {
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;

    int aaBitSum  = (fromR) ? 0 : totalAABits;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;

    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA * !fromR;
        dnaBitSum -= bitPerCodon * !fromR;

        const int myAA   = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;

        const int idx1 = matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCode->aa2codon[myAA][codon1].c_str())
        ];
        const int idx2 = matrix.aa2num[
            translateNucl->translateSingleCodon(
                geneticCode->aa2codon[myAA][codon2].c_str())
        ];

        score += matrix.subMatrix[idx1][idx2];

        aaBitSum += bitPerAA * fromR;
        dnaBitSum += bitPerCodon * fromR;
    }

    return score;
}

float SingleCodePattern::hammingDistScore(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count,
    bool fromR) const 
{
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum  = (fromR) ? 0 : totalAABits;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;
    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA * !fromR;
        dnaBitSum -= bitPerCodon * !fromR;
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        const int hammingDist = geneticCode->getHammingDist(aa, codon1, codon2);
        if (hammingDist == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * hammingDist;
        }
        aaBitSum += bitPerAA * fromR;
        dnaBitSum += bitPerCodon * fromR;
    }  
    return score;
}

uint8_t SingleCodePattern::hammingDistSum(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int count,
    bool fromR) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = (fromR) ? 0 : totalAABits;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;
    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA * !fromR;
        dnaBitSum -= bitPerCodon * !fromR;
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        hammingSum += geneticCode->getHammingDist(aa, codon1, codon2);
        aaBitSum += bitPerAA * fromR;
        dnaBitSum += bitPerCodon * fromR;
    }  
    return hammingSum;
}

uint8_t SingleCodePattern::hammingDistSum(
    uint64_t kmer1, 
    uint64_t kmer2) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = 0; i < kmerLen; ++i) {
        const int aa = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        hammingSum += geneticCode->getHammingDist(aa, codon1, codon2);
        aaBitSum += bitPerAA;
        dnaBitSum += bitPerCodon;
    }  
    return hammingSum;
}

MultiCodePattern::MultiCodePattern(const std::string & customFile) {
    std::cout << "Loading multi-code metamer pattern from file: " << customFile << std::endl;
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

MatchScore MultiCodePattern::calMatchScore(
    uint64_t kmer1,
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix& matrix,
    bool fromR) const
{
    float idScore = 0.0f;
    float subScore = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;

    int aaBitSum  = fromR ? 0 : totalAABits;
    int dnaBitSum = fromR ? 0 : totalDNABits;

    for (int step = 0; step < count; ++step) {
        const int i = fromR
            ? (kmerLen - 1 - step)
            : step;

        const int currAABit  = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];

        aaBitSum -= currAABit * (!fromR);
        dnaBitSum -= currDnaBit * (!fromR);

        const int myAA   = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);

        const auto& gc = geneticCodes[codePattern[i]];

        const int idx1 = matrix.aa2num[
            translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon1].c_str())
        ];
        const int idx2 = matrix.aa2num[
            translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon2].c_str())
        ];

        const int hammingDist = gc->getHammingDist(myAA, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        subScore += matrix.subMatrix[idx1][idx2];

        aaBitSum += currAABit * fromR;
        dnaBitSum += currDnaBit * fromR;
    }

    return {idScore, subScore};
}

float MultiCodePattern::substitutionScore(
    uint64_t kmer1,
    uint64_t kmer2,
    int count,
    const SubstitutionMatrix& matrix,
    bool fromR
) const {
    float totalScore = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;

    int aaBitSum  = fromR ? 0 : totalAABits;
    int dnaBitSum = fromR ? 0 : totalDNABits;

    for (int step = 0; step < count; ++step) {
        const int i = fromR
            ? (kmerLen - 1 - step)
            : step;

        const int currAABit  = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];

        aaBitSum -= currAABit * (!fromR);
        dnaBitSum -= currDnaBit * (!fromR);

        const int myAA   = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);

        const auto& gc = geneticCodes[codePattern[i]];

        const int idx1 = matrix.aa2num[
            translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon1].c_str())
        ];
        const int idx2 = matrix.aa2num[
            translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon2].c_str())
        ];

        totalScore += matrix.subMatrix[idx1][idx2];
        aaBitSum += currAABit * fromR;
        dnaBitSum += currDnaBit * fromR;
    }

    return totalScore;
}

float MultiCodePattern::hammingDistScore(
    uint64_t kmer1, 
    uint64_t kmer2,
    int count,
    bool fromR) const 
{
    float score = 0.0f;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum  = fromR ? 0 : totalAABits;
    int dnaBitSum = fromR ? 0 : totalDNABits;
    for (int step = 0; step < count; ++step) {
        const int i = fromR
            ? (kmerLen - 1 - step)
            : step;
        const int currAABit  = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];
        aaBitSum -= currAABit * (!fromR);
        dnaBitSum -= currDnaBit * (!fromR);
        const int aa = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int hammingDist = geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
        if (hammingDist == 0) {
            score += 3.0f;
        } else {
            score += 2.0f - 0.5f * hammingDist;
        }
        aaBitSum += currAABit * fromR;
        dnaBitSum += currDnaBit * fromR;
    }
    return score;
}

uint8_t MultiCodePattern::hammingDistSum(
    uint64_t kmer1, 
    uint64_t kmer2, 
    int count,
    bool fromR) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = (fromR) ? 0 : totalAABits;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;
    for (int step = 0; step < count; ++step) {
        const int i = fromR
            ? (kmerLen - 1 - step)
            : step;
        const int currAABit  = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];
        aaBitSum -= currAABit * (!fromR);
        dnaBitSum -= currDnaBit * (!fromR);
        const int aa = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        hammingSum += geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
        aaBitSum += currAABit * fromR;
        dnaBitSum += currDnaBit * fromR;
    }  
    return hammingSum;
}


uint8_t MultiCodePattern::hammingDistSum(
    uint64_t kmer1, 
    uint64_t kmer2) const 
{
    uint8_t hammingSum = 0;
    const uint64_t aaPart = kmer1 >> totalDNABits;  
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;
    for (int i = 0; i < kmerLen; ++i) {
        const int currAABit  = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];
        aaBitSum -= currAABit;
        dnaBitSum -= currDnaBit;
        const int aa = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        hammingSum += geneticCodes[codePattern[i]]->getHammingDist(aa, codon1, codon2);
    }  
    return hammingSum;
}