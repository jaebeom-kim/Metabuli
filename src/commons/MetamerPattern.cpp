#include "MetamerPattern.h"
#include "json.h"
#include <array>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#endif
#if defined(__SSSE3__)
#include <tmmintrin.h>
#endif

using json = nlohmann::json;

namespace {
constexpr size_t MAX_BATCH_KMER_LEN = 64;

inline bool usesRegularEightCodonLookup(const GeneticCode *geneticCode) {
    return dynamic_cast<const RegularGeneticCode *>(geneticCode) != nullptr;
}

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
size_t hammingDistSumBatchRegularNeon(
    uint64_t kmer1,
    const Kmer *candidates,
    size_t candidateCnt,
    uint8_t *hammings,
    int kmerLen,
    const uint8_t *lookup) {
    if (kmerLen > static_cast<int>(MAX_BATCH_KMER_LEN)) {
        return 0;
    }

    uint8x8_t rows[MAX_BATCH_KMER_LEN];
    uint64_t queryCodons = kmer1;
    for (int pos = 0; pos < kmerLen; ++pos) {
        rows[pos] = vld1_u8(lookup + ((queryCodons & 0x7U) << 3U));
        queryCodons >>= 3U;
    }

    size_t i = 0;
    for (; i + 8 <= candidateCnt; i += 8) {
        uint64_t streams[8] = {
            candidates[i + 0].value,
            candidates[i + 1].value,
            candidates[i + 2].value,
            candidates[i + 3].value,
            candidates[i + 4].value,
            candidates[i + 5].value,
            candidates[i + 6].value,
            candidates[i + 7].value,
        };
        uint8x8_t sums = vdup_n_u8(0);

        for (int pos = 0; pos < kmerLen; ++pos) {
            alignas(8) uint8_t targetCodons[8] = {
                static_cast<uint8_t>(streams[0] & 0x7U),
                static_cast<uint8_t>(streams[1] & 0x7U),
                static_cast<uint8_t>(streams[2] & 0x7U),
                static_cast<uint8_t>(streams[3] & 0x7U),
                static_cast<uint8_t>(streams[4] & 0x7U),
                static_cast<uint8_t>(streams[5] & 0x7U),
                static_cast<uint8_t>(streams[6] & 0x7U),
                static_cast<uint8_t>(streams[7] & 0x7U),
            };

            streams[0] >>= 3U;
            streams[1] >>= 3U;
            streams[2] >>= 3U;
            streams[3] >>= 3U;
            streams[4] >>= 3U;
            streams[5] >>= 3U;
            streams[6] >>= 3U;
            streams[7] >>= 3U;

            sums = vadd_u8(sums, vtbl1_u8(rows[pos], vld1_u8(targetCodons)));
        }

        vst1_u8(hammings + i, sums);
    }

    return i;
}
#endif

#if defined(__SSSE3__)
size_t hammingDistSumBatchRegularSsse3(
    uint64_t kmer1,
    const Kmer *candidates,
    size_t candidateCnt,
    uint8_t *hammings,
    int kmerLen,
    const uint8_t *lookup) {
    if (kmerLen > static_cast<int>(MAX_BATCH_KMER_LEN)) {
        return 0;
    }

    __m128i rows[MAX_BATCH_KMER_LEN];
    uint64_t queryCodons = kmer1;
    for (int pos = 0; pos < kmerLen; ++pos) {
        alignas(16) uint8_t row[16] = {};
        const uint8_t *src = lookup + ((queryCodons & 0x7U) << 3U);
        for (int j = 0; j < 8; ++j) {
            row[j] = src[j];
        }
        rows[pos] = _mm_load_si128(reinterpret_cast<const __m128i *>(row));
        queryCodons >>= 3U;
    }

    size_t i = 0;
    for (; i + 16 <= candidateCnt; i += 16) {
        uint64_t streams[16] = {
            candidates[i + 0].value,
            candidates[i + 1].value,
            candidates[i + 2].value,
            candidates[i + 3].value,
            candidates[i + 4].value,
            candidates[i + 5].value,
            candidates[i + 6].value,
            candidates[i + 7].value,
            candidates[i + 8].value,
            candidates[i + 9].value,
            candidates[i + 10].value,
            candidates[i + 11].value,
            candidates[i + 12].value,
            candidates[i + 13].value,
            candidates[i + 14].value,
            candidates[i + 15].value,
        };
        __m128i sums = _mm_setzero_si128();

        for (int pos = 0; pos < kmerLen; ++pos) {
            alignas(16) uint8_t targetCodons[16] = {
                static_cast<uint8_t>(streams[0] & 0x7U),
                static_cast<uint8_t>(streams[1] & 0x7U),
                static_cast<uint8_t>(streams[2] & 0x7U),
                static_cast<uint8_t>(streams[3] & 0x7U),
                static_cast<uint8_t>(streams[4] & 0x7U),
                static_cast<uint8_t>(streams[5] & 0x7U),
                static_cast<uint8_t>(streams[6] & 0x7U),
                static_cast<uint8_t>(streams[7] & 0x7U),
                static_cast<uint8_t>(streams[8] & 0x7U),
                static_cast<uint8_t>(streams[9] & 0x7U),
                static_cast<uint8_t>(streams[10] & 0x7U),
                static_cast<uint8_t>(streams[11] & 0x7U),
                static_cast<uint8_t>(streams[12] & 0x7U),
                static_cast<uint8_t>(streams[13] & 0x7U),
                static_cast<uint8_t>(streams[14] & 0x7U),
                static_cast<uint8_t>(streams[15] & 0x7U),
            };

            streams[0] >>= 3U;
            streams[1] >>= 3U;
            streams[2] >>= 3U;
            streams[3] >>= 3U;
            streams[4] >>= 3U;
            streams[5] >>= 3U;
            streams[6] >>= 3U;
            streams[7] >>= 3U;
            streams[8] >>= 3U;
            streams[9] >>= 3U;
            streams[10] >>= 3U;
            streams[11] >>= 3U;
            streams[12] >>= 3U;
            streams[13] >>= 3U;
            streams[14] >>= 3U;
            streams[15] >>= 3U;

            const __m128i idx =
                _mm_load_si128(reinterpret_cast<const __m128i *>(targetCodons));
            sums = _mm_add_epi8(sums, _mm_shuffle_epi8(rows[pos], idx));
        }

        _mm_storeu_si128(reinterpret_cast<__m128i *>(hammings + i), sums);
    }

    return i;
}
#endif
}

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


int getWeightedPosNum(const std::string & customFile) {
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
    int weightedPosNum = j["length"];
    file.close();
    return weightedPosNum;
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
    windowSize = kmerLen;
    std::cout << "K-mer length: " << kmerLen << std::endl;

    for (const auto& codeEntry : j["codes"]) {

        const auto& codons_json = codeEntry["codons"];

        std::cout << "Loading custom genetic code: " << codeEntry["id"] << std::endl;

        std::vector<std::pair<std::string, std::vector<std::string>>> translationTable;

        for (const auto& entry : codons_json) {
            std::string aaStr = entry[0].get<std::string>();
            std::vector<std::string> codons = entry[1].get<std::vector<std::string>>();
            for (size_t i = 1; i < aaStr.size(); ++i) {
                aaFreq[aaStr[0] - 'A'] += aaFreq[aaStr[i] - 'A'];
                aaFreq[aaStr[i] - 'A'] = 0.0;
            }
            translationTable.emplace_back(std::string{aaStr[0]}, codons);
        }
        geneticCode = std::make_unique<CustomGeneticCode>(translationTable);
    }
    initializeLnFreq();

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
    std::cout << "Total DNA bits: " << totalDNABits << ", Total AA bits: " << totalAABits << std::endl;
    dnaMask = (1ULL << totalDNABits) - 1;
    codonMask = (1ULL << bitPerCodon) - 1;
    aaMask = (1ULL << bitPerAA) - 1;
    file.close();
}

MatchScore SingleCodePattern::calMatchScore(
    uint64_t kmer1, 
    uint64_t kmer2) const 
{
    float idScore = 0.0f;
    // float subScore = 0.0f;
    double prob = 0.0;
    const uint64_t aaPart = kmer1 >> totalDNABits;
    int aaBitSum = 0;
    int dnaBitSum = 0;
    for (int i = 0; i < kmerLen; ++i) { // kmerLen is correct. Not windowSize 
        const int myAA = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        // const int idx1 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(
        //         geneticCode->aa2codon[myAA][codon1].c_str())
        // ];
        // const int idx2 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(
        //         geneticCode->aa2codon[myAA][codon2].c_str())
        // ];

        const int hammingDist = geneticCode->getHammingDist(myAA, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        // subScore += matrix.subMatrix[idx1][idx2];
        prob += lnFreq[geneticCode->aminoacids[myAA] - 'A'];

        aaBitSum += bitPerAA;
        dnaBitSum += bitPerCodon;
        
    }  
    return {idScore, 0.0f, prob};
}

MatchScore SingleCodePattern::calMatchScore(
    uint64_t kmer1,
    uint64_t kmer2,
    uint32_t count,
    bool fromR) const 
{
    float idScore = 0.0f;
    // float subScore = 0.0f;
    double prob = 0.0;
    const uint64_t aaPart = kmer1 >> totalDNABits;

    int aaBitSum  = (fromR) ? 0 : totalAABits;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;

    for (int i = 0; i < count; ++i) {
        aaBitSum -= bitPerAA * !fromR;
        dnaBitSum -= bitPerCodon * !fromR;

        const int myAA   = (aaPart >> aaBitSum) & aaMask;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;

        // const int idx1 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(
        //         geneticCode->aa2codon[myAA][codon1].c_str())
        // ];
        // const int idx2 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(
        //         geneticCode->aa2codon[myAA][codon2].c_str())
        // ];

        const int hammingDist = geneticCode->getHammingDist(myAA, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        // subScore += matrix.subMatrix[idx1][idx2];
        prob += lnFreq[geneticCode->aminoacids[myAA] - 'A'];

        aaBitSum += bitPerAA * fromR;
        dnaBitSum += bitPerCodon * fromR;
    }
    return {idScore, 0.0f, prob};
}

MatchScore SingleCodePattern::calMatchScore(
    uint64_t aa, 
    uint64_t codon1, 
    uint64_t codon2, 
    uint32_t validPosMask) const  
{
    float idScore = 0.0f;
    double prob = 0.0;
    while (validPosMask) {
        const int i = __builtin_ctz(validPosMask);
        validPosMask &= validPosMask - 1;

        const int sAA = bitPerAA * i;
        const int sDNA = bitPerCodon * i;

        const int myAA   = (aa >> sAA)  & aaMask;
        const int myCodon1 = (codon1  >> sDNA) & codonMask;
        const int myCodon2 = (codon2  >> sDNA) & codonMask;
        
        const int hammingDist = geneticCode->getHammingDist(myAA, myCodon1, myCodon2);
      
        idScore  += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        prob     += lnFreq[geneticCode->aminoacids[myAA] - 'A'];
    }

    return {idScore, 0.0f, prob};

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
    
    // 1. Prepare Streams
    // Extract the AA part to a separate variable so we can shift it independently.
    uint64_t aaStream = kmer1 >> totalDNABits; 
    
    // We will shift these destructively. 
    // It is safe that kmer1/kmer2 have AA bits at the top because the loop 
    // terminates before those bits are shifted down into the codon position.

    // 2. Cache Constants 
    // Loading these into local consts helps the compiler keep them in registers.
    const int bAA = bitPerAA;
    const int bCodon = bitPerCodon;

    // 3. Loop (LSB to MSB)
    for (int i = 0; i < kmerLen; ++i) { 
        // Step A: Extract Bottom Bits (Fast AND)
        const int aa = aaStream & aaMask;
        const int c1 = kmer1 & codonMask;
        const int c2 = kmer2 & codonMask;

        // Step B: Calculate & Sum
        // Note: Ensure geneticCode->getHammingDist is INLINED for max speed.
        hammingSum += geneticCode->getHammingDist(aa, c1, c2);

        // Step C: Shift Streams Down (Fast Constant Shift)
        // This brings the next codon/AA into the LSB position for the next iteration.
        aaStream >>= bAA;
        kmer1 >>= bCodon;
        kmer2 >>= bCodon;
    }  
    return hammingSum;
}

void SingleCodePattern::hammingDistSumBatch(
    uint64_t kmer1,
    const Kmer *candidates,
    size_t candidateCnt,
    uint8_t *hammings) const {
    if (candidateCnt == 0) {
        return;
    }

    const bool regularLookup = usesRegularEightCodonLookup(geneticCode.get());
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    size_t i = 0;
    if (regularLookup && bitPerCodon == 3 && codonMask == 0x7U) {
        i = hammingDistSumBatchRegularNeon(
            kmer1, candidates, candidateCnt, hammings, kmerLen, geneticCode->hammingLookup);
    }
#elif defined(__SSSE3__)
    size_t i = 0;
    if (regularLookup && bitPerCodon == 3 && codonMask == 0x7U) {
        i = hammingDistSumBatchRegularSsse3(
            kmer1, candidates, candidateCnt, hammings, kmerLen, geneticCode->hammingLookup);
    }
#else
    size_t i = 0;
#endif

    if (kmerLen > static_cast<int>(MAX_BATCH_KMER_LEN)) {
        for (; i < candidateCnt; ++i) {
            hammings[i] = hammingDistSum(kmer1, candidates[i].value);
        }
        return;
    }

    std::array<const uint8_t *, MAX_BATCH_KMER_LEN> rows{};
    uint64_t aaStream = kmer1 >> totalDNABits;
    uint64_t queryCodons = kmer1;
    const uint8_t *lookup = geneticCode->hammingLookup;
    const int maxCodon = geneticCode->maxCodonPerAA;
    const int rowStride = maxCodon * maxCodon;

    for (int pos = 0; pos < kmerLen; ++pos) {
        const int aa = aaStream & aaMask;
        const int queryCodon = queryCodons & codonMask;
        rows[pos] = regularLookup
            ? lookup + queryCodon * 8
            : lookup + aa * rowStride + queryCodon * maxCodon;
        aaStream >>= bitPerAA;
        queryCodons >>= bitPerCodon;
    }

    for (; i < candidateCnt; ++i) {
        uint64_t targetCodons = candidates[i].value;
        uint8_t hammingSum = 0;
        for (int pos = 0; pos < kmerLen; ++pos) {
            hammingSum += rows[pos][targetCodons & codonMask];
            targetCodons >>= bitPerCodon;
        }
        hammings[i] = hammingSum;
    }
}

// uint8_t SingleCodePattern::hammingDistSum(
//     uint64_t kmer1, 
//     uint64_t kmer2) const 
// {
//     uint8_t hammingSum = 0;
//     const uint64_t aaPart = kmer1 >> totalDNABits;  
//     int aaBitSum = 0;
//     int dnaBitSum = 0;
//     for (int i = 0; i < kmerLen; ++i) { // kmerLen is correct. Not windowSize 
//         const int aa = (aaPart >> aaBitSum) & aaMask;
//         const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
//         const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
//         hammingSum += geneticCode->getHammingDist(aa, codon1, codon2);
//         aaBitSum += bitPerAA;
//         dnaBitSum += bitPerCodon;
//     }  
//     return hammingSum;
// }



uint8_t SpacedPattern::hammingDistSum(
    uint64_t kmer1,    // Query k-mer
    uint64_t kmer2,    // Target k-mer
    int count,         // Window shift count
    bool fromR) const 
{
    const uint64_t aaPart   = kmer1 >> totalDNABits;

    uint32_t validPosMask;
    if (fromR) {
        validPosMask = spaceMask & ~(spaceMask << count);
    } else {
        validPosMask = spaceMask & ~(spaceMask >> count);
    }
    uint8_t hammingSum = 0; 

    while (validPosMask) {
        int i = __builtin_ctz(validPosMask); // Get index of next set bit
        validPosMask &= validPosMask - 1; // Clear the processed bit

        const int sAA = aaShifts[i];
        const int sDNA = dnaShifts[i];

        const int myAA   = (aaPart  >> sAA)  & aaMask;
        const int codon1 = (kmer1   >> sDNA) & codonMask;
        const int codon2 = (kmer2   >> sDNA) & codonMask;

        hammingSum += geneticCode->getHammingDist(myAA, codon1, codon2);
    }

    return hammingSum;
}


MatchScore SpacedPattern::calMatchScore(
    uint64_t kmer1,
    uint64_t kmer2,
    uint32_t validPosMask,
    bool fromR) const 
{
    float idScore = 0.0f;
    // float subScore = 0.0f;
    double prob = 0.0;
    const uint64_t aaPart = kmer1 >> totalDNABits;
    while (validPosMask) {
        const int i = __builtin_ctz(validPosMask);
        validPosMask &= validPosMask - 1;

        const int sAA = aaShifts[i];
        const int sDNA = dnaShifts[i];

        const int myAA   = (aaPart >> sAA)  & aaMask;
        const int codon1 = (kmer1  >> sDNA) & codonMask;
        const int codon2 = (kmer2  >> sDNA) & codonMask;
        
        // const int idx1 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(
        //         geneticCode->aa2codon[myAA][codon1].c_str())
        // ];

        // const int idx2 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(
        //         geneticCode->aa2codon[myAA][codon2].c_str())
        // ];

        const int hammingDist = geneticCode->getHammingDist(myAA, codon1, codon2);
        idScore  += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        // subScore += matrix.subMatrix[idx1][idx2];
        prob     += lnFreq[geneticCode->aminoacids[myAA] - 'A'];
    }

    return {idScore, 0.0f, prob};
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
    initializeLnFreq();

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


/* MultiCodePattern member functions */

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
    uint32_t count,
    bool fromR) const
{
    float idScore = 0.0f;
    // float subScore = 0.0f;
    double prob = 0.0;
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

        // const int idx1 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon1].c_str())
        // ];
        // const int idx2 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon2].c_str())
        // ];

        const int hammingDist = gc->getHammingDist(myAA, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        // subScore += matrix.subMatrix[idx1][idx2];
        prob += lnFreq[gc->aminoacids[myAA] - 'A'];

        aaBitSum += currAABit * fromR;
        dnaBitSum += currDnaBit * fromR;
    }

    return {idScore, 0.0f, prob};
}

MatchScore MultiCodePattern::calMatchScore(
    uint64_t kmer1,
    uint64_t kmer2) const
{
    float idScore = 0.0f;
    // float subScore = 0.0f;
    double prob = 0.0;
    const uint64_t aaPart = kmer1 >> totalDNABits;
    int aaBitSum = totalAABits;
    int dnaBitSum = totalDNABits;

    for (int i = 0; i < kmerLen; ++i) {
        const int currAABit  = aaBitList[codePattern[i]];
        const int currDnaBit = codonBitList[codePattern[i]];
        aaBitSum -= currAABit;
        dnaBitSum -= currDnaBit;
        const int myAA   = (aaPart >> aaBitSum) & ((1 << currAABit) - 1);
        const int codon1 = (kmer1 >> dnaBitSum) & ((1 << currDnaBit) - 1);
        const int codon2 = (kmer2 >> dnaBitSum) & ((1 << currDnaBit) - 1);

        const auto& gc = geneticCodes[codePattern[i]];

        // const int idx1 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon1].c_str())
        // ];
        // const int idx2 = matrix.aa2num[
        //     translateNucl->translateSingleCodon(gc->aa2codon[myAA][codon2].c_str())
        // ];

        const int hammingDist = gc->getHammingDist(myAA, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        // subScore += matrix.subMatrix[idx1][idx2];
        prob += lnFreq[gc->aminoacids[myAA] - 'A'];
    }

    return {idScore, 0.0f, prob};
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


uint8_t LegacyPattern::hammingDistSum(uint64_t kmer1, uint64_t kmer2, int count, bool fromR) const {
    uint8_t hammingSum = 0;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;
    for (int i = 0; i < count; ++i) {
        dnaBitSum -= bitPerCodon * !fromR;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        hammingSum += geneticCode->getHammingDist(0, codon1, codon2);
        dnaBitSum += bitPerCodon * fromR;
    }  
    return hammingSum;
}

uint8_t LegacyPattern::hammingDistSum(uint64_t kmer1, uint64_t kmer2) const {
    uint8_t hammingSum = 0;
    const int bCodon = bitPerCodon;
    for (int i = 0; i < kmerLen; ++i) {
        hammingSum += geneticCode->getHammingDist(0, kmer1 & codonMask, kmer2 & codonMask);
        kmer1 >>= bCodon;
        kmer2 >>= bCodon;
    }  
    return hammingSum;

}

void LegacyPattern::hammingDistSumBatch(
    uint64_t kmer1,
    const Kmer *candidates,
    size_t candidateCnt,
    uint8_t *hammings) const {
    if (candidateCnt == 0) {
        return;
    }

    const bool regularLookup = usesRegularEightCodonLookup(geneticCode.get());
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    size_t i = 0;
    if (regularLookup && bitPerCodon == 3 && codonMask == 0x7U) {
        i = hammingDistSumBatchRegularNeon(
            kmer1, candidates, candidateCnt, hammings, kmerLen, geneticCode->hammingLookup);
    }
#elif defined(__SSSE3__)
    size_t i = 0;
    if (regularLookup && bitPerCodon == 3 && codonMask == 0x7U) {
        i = hammingDistSumBatchRegularSsse3(
            kmer1, candidates, candidateCnt, hammings, kmerLen, geneticCode->hammingLookup);
    }
#else
    size_t i = 0;
#endif

    if (kmerLen > static_cast<int>(MAX_BATCH_KMER_LEN)) {
        for (; i < candidateCnt; ++i) {
            hammings[i] = hammingDistSum(kmer1, candidates[i].value);
        }
        return;
    }

    std::array<const uint8_t *, MAX_BATCH_KMER_LEN> rows{};
    uint64_t queryCodons = kmer1;
    const uint8_t *lookup = geneticCode->hammingLookup;
    const int maxCodon = geneticCode->maxCodonPerAA;

    for (int pos = 0; pos < kmerLen; ++pos) {
        const int queryCodon = queryCodons & codonMask;
        rows[pos] = regularLookup
            ? lookup + queryCodon * 8
            : lookup + queryCodon * maxCodon;
        queryCodons >>= bitPerCodon;
    }

    for (; i < candidateCnt; ++i) {
        uint64_t targetCodons = candidates[i].value;
        uint8_t hammingSum = 0;
        for (int pos = 0; pos < kmerLen; ++pos) {
            hammingSum += rows[pos][targetCodons & codonMask];
            targetCodons >>= bitPerCodon;
        }
        hammings[i] = hammingSum;
    }
}

MatchScore LegacyPattern::calMatchScore(uint64_t kmer1, uint64_t kmer2, uint32_t count, bool fromR) const {
    float idScore = 0.0f;
    int dnaBitSum = (fromR) ? 0 : totalDNABits;
    for (int i = 0; i < count; ++i) {
        dnaBitSum -= bitPerCodon * !fromR;
        const int codon1 = (kmer1 >> dnaBitSum) & codonMask;
        const int codon2 = (kmer2 >> dnaBitSum) & codonMask;
        const int hammingDist = geneticCode->getHammingDist(0, codon1, codon2);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        dnaBitSum += bitPerCodon * fromR;
    }
    return {idScore, 0.0f, 0.0};
}
MatchScore LegacyPattern::calMatchScore(uint64_t kmer1, uint64_t kmer2) const {
    float idScore = 0.0f;
    const int bCodon = bitPerCodon;
    for (int i = 0; i < kmerLen; ++i) {
        const int hammingDist = geneticCode->getHammingDist(0, kmer1 & codonMask, kmer2 & codonMask);
        idScore += (hammingDist == 0) ? 3.0f : (2.0f - 0.5f * hammingDist);
        kmer1 >>= bCodon;
        kmer2 >>= bCodon;
    }
    return {idScore, 0.0f, 0.0};
}
