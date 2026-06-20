#include "Classifier.h"
#include "FileUtil.h"
#include "QueryIndexer.h"
#include "common.h"

#include <variant>

#if defined(__ARM_NEON) || defined(__ARM_NEON__)
#include <arm_neon.h>
#elif defined(__SSE2__)
#include <emmintrin.h>
#endif

namespace {
inline void saturatingAddBins(uint8_t *dst, const uint8_t *src, size_t count) {
    size_t i = 0;
#if defined(__ARM_NEON) || defined(__ARM_NEON__)
    for (; i + 16 <= count; i += 16) {
        uint8x16_t a = vld1q_u8(dst + i);
        uint8x16_t b = vld1q_u8(src + i);
        vst1q_u8(dst + i, vqaddq_u8(a, b));
    }
#elif defined(__SSE2__)
    for (; i + 16 <= count; i += 16) {
        __m128i a = _mm_loadu_si128(reinterpret_cast<const __m128i *>(dst + i));
        __m128i b = _mm_loadu_si128(reinterpret_cast<const __m128i *>(src + i));
        _mm_storeu_si128(reinterpret_cast<__m128i *>(dst + i), _mm_adds_epu8(a, b));
    }
#endif
    for (; i < count; ++i) {
        const int sum = dst[i] + src[i];
        dst[i] = (sum > 255) ? 255 : static_cast<uint8_t>(sum);
    }
}

inline void addScoreToClade(
    TaxonomyWrapper *taxonomy,
    TaxID taxId,
    double score,
    std::unordered_map<TaxID, double> &cladeScoreSums)
{
    if (taxId == 0 || !taxonomy->nodeExists(taxId)) {
        return;
    }

    const TaxonNode *taxon = taxonomy->taxonNode(taxId);
    cladeScoreSums[taxon->taxId] += score;
    while (taxon->parentTaxId != taxon->taxId && taxonomy->nodeExists(taxon->parentTaxId)) {
        taxon = taxonomy->taxonNode(taxon->parentTaxId);
        cladeScoreSums[taxon->taxId] += score;
    }
}

inline void addAverageScoresToMap(
    const std::unordered_map<TaxID, double> &cladeScoreSums,
    const std::unordered_map<TaxID, TaxonCounts> &cladeCounts,
    std::unordered_map<TaxID, double> &avgScores)
{
    avgScores.reserve(avgScores.size() + cladeScoreSums.size());
    for (const auto &taxScore : cladeScoreSums) {
        auto countIt = cladeCounts.find(taxScore.first);
        if (countIt == cladeCounts.end() || countIt->second.cladeCount == 0) {
            continue;
        }

        avgScores[taxScore.first] = taxScore.second / static_cast<double>(countIt->second.cladeCount);
    }
}
}

Classifier::Classifier(LocalParameters & par) : par(par) {
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    matchPerKmer = par.matchPerKmer;
    loadDbParameters(par, par.filenames[1 + (par.seqMode == 2)]);
    kmerFormat = par.kmerFormat;
    if (par.dbTotalLength == 0) {
        par.dbTotalLength = readDbSize(par.filenames[1 + (par.seqMode == 2)]);
    }

    cout << "Database name : " << par.dbName << endl;
    cout << "Creation date : " << par.dbDate << endl;
    
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    if (par.precisionMode != 0) {
        preciseModePreset(par);
    }

    if (kmerFormat == 1) {
        metamerPattern = new LegacyPattern(std::make_unique<RegularGeneticCode>(), 8);
        cerr << "Warning: the specified database uses an old k-mer format." << endl;
        cerr << "   E-value calculation is not supported." << endl;
        par.maxEValue = 0;
    } else if (!par.customMetamer.empty()) {
        int codeNum = getCodeNum(par.customMetamer);
        if (codeNum == 1) {
            if (par.spaceMask.empty()) {
                cout << "Using SingleCodePattern with custom metamer." << endl;
                metamerPattern = new SingleCodePattern(par.customMetamer);
            } else {
                cout << "Using SpacedPattern with custom metamer." << endl;
                uint32_t mask = parseMask(par.spaceMask.c_str());
                metamerPattern = new SpacedPattern(par.customMetamer, mask);
            }
        } else if (codeNum > 1) {
            metamerPattern = new MultiCodePattern(par.customMetamer);
        }
    } else if (kmerFormat == 2) {
        if (par.spaceMask.empty()) {
            cout << "Using SingleCodePattern with RegularGeneticCode" << endl;
            metamerPattern = new SingleCodePattern(std::make_unique<RegularGeneticCode>(), 8);
        } else {
            cout << "Using SpacedPattern with RegularGeneticCode" << endl;
            uint32_t mask = parseMask(par.spaceMask.c_str());
            metamerPattern = new SpacedPattern(std::make_unique<RegularGeneticCode>(), __builtin_popcount(mask), mask);
        }
    }

    if (par.storeKmerPos) {
        C_LOG2_C.resize(256);
        for (int i = 1; i < 256; ++i) {
            C_LOG2_C[i] = i * std::log2(static_cast<double>(i));
        }
    }

    kmerExtractor = new KmerExtractor(par, metamerPattern);
    queryIndexer = new QueryIndexer(par);
    kmerMatcher = new KmerMatcher(par, taxonomy, metamerPattern);
    reporter = new Reporter(par, taxonomy);
}

Classifier::~Classifier() {
    delete taxonomy;
    delete queryIndexer;
    delete kmerExtractor;
    delete kmerMatcher;
    delete reporter;
    delete geneticCode;
    delete metamerPattern;
    if (mappingResList) delete[] mappingResList;
}


uint64_t Classifier::calculateBufferSize(
    uint64_t queryListSize,
    uint64_t matchPerKmer,
    int mode) 
{
    int matchSize = mode == 1 ? sizeof(Match) : sizeof(MatchWithPos);

    size_t totalBytes = (size_t) par.ramUsage * 1024 * 1024 * 1024;
    size_t bytesPerThread = 
        1024 * 1024 * matchSize + // local match buffer 32 MB
        1024 * 1024 * 4 * sizeof(Kmer) +  // DeltaIdxReader::valueBuffer 64 MB
        1024 * 1024 * 4 * sizeof(uint16_t) +   // DeltaIdxReader::deltaBuffer 2 MB
        1024 * 1024 * 4 * sizeof(uint32_t);    // DeltaIdxReader::infoBuffer 4 MB
    
    if (mode != 1) {
        bytesPerThread += 1024 * 1024 * 4 * sizeof(uint16_t); // DeltaIdxReader::posBuffer 2 MB
    }

    size_t overhead = 128 * 1024 * 1024; // 128MB
    size_t queryListBytes = queryListSize * (sizeof(Query) + 150); //  104,857,600
    size_t availableBytes = totalBytes - (par.threads * bytesPerThread) - overhead - queryListBytes; 
    size_t bytePerKmer = sizeof(Kmer) + matchPerKmer * matchSize;
    size_t totalSize = availableBytes / bytePerKmer;
    return totalSize;
}
    
void Classifier::preciseModePreset(LocalParameters & par) {
    uint32_t mask = parseMask(par.spaceMask.c_str());
    size_t windowSize = par.spaceMask.length();
    size_t kmerLen = __builtin_popcount(mask);

    float minScoreCp = par.minScore;
    float minSpScoreCp = par.minSpScore;
    double maxEValueCp = par.maxEValue;

    if (par.precisionMode == 1) { // short-read preset
        if (windowSize == kmerLen || kmerLen == 0) {
            std::cout << "Using short-read presets for contiguous k-mer search: " << std::endl;
            std::cout << "   --min-score 0.15 --min-sp-score 0.5 -e 0.001" << std::endl;
            par.minScore = 0.15;
            par.minSpScore = 0.5;
            par.maxEValue = 0.001;
        } else {
            std::cout << "Using short-read presets for spaced k-mer search: " << std::endl;
            std::cout << "   --min-score 0.2 --min-sp-score 0.6 -e 0.001" << std::endl;
            par.minScore = 0.2;
            par.minSpScore = 0.6;
            par.maxEValue = 0.001;
        }
    } else if (par.precisionMode == 2) { // HiFi long-read preset
        std::cout << "Using HiFi long-read presets: " << std::endl;
        std::cout << "   --min-score 0.07 --min-sp-score 0.3 -e 0.001" << std::endl;
        par.minScore = 0.07;
        par.minSpScore = 0.3;
        par.maxEValue = 0.001;
    }

    if (minScoreCp != 0 && minScoreCp != par.minScore) {
        std::cout << "Overriding preset --min-score " << par.minScore << " with user specified value " << minScoreCp << std::endl;
        par.minScore = minScoreCp;
    }
    if (minSpScoreCp != 0 && minSpScoreCp != par.minSpScore) {
        std::cout << "Overriding preset --min-sp-score " << par.minSpScore << " with user specified value " << minSpScoreCp << std::endl;
        par.minSpScore = minSpScoreCp;
    }
    if (maxEValueCp != 1 && maxEValueCp != par.maxEValue) {
        std::cout << "Overriding preset --e-value " << par.maxEValue << " with user specified value " << maxEValueCp << std::endl;
        par.maxEValue = maxEValueCp;
    }

}



void Classifier::classifyReads() {
    Buffer<Kmer> queryKmerBuffer;
    Buffer<Match> matchBuffer;
    vector<Query> queryList;
    queryList.reserve(512 * 1024);

    reporter->openReadClassificationFile();

    bool complete = false;
    SeqEntry * savedSeq_1 = new SeqEntry();
    SeqEntry * savedSeq_2 = new SeqEntry(); 
    uint64_t processedReadCnt = 0;
    std::unordered_map<TaxID, double> cladeScoreSums;

    std::cout << "--------------------" << std::endl;
    while (!complete) {
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = par.seqMode == 2 ? KSeqFactory(par.filenames[1].c_str()) : nullptr;

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        size_t queryKmerBufferSize = calculateBufferSize(queryList.size(), matchPerKmer, 1);
        queryKmerBuffer.reallocateMemory(queryKmerBufferSize);
        matchBuffer.reallocateMemory(queryKmerBufferSize * matchPerKmer);

        bool moreReads = true;
        while (moreReads) {
            queryList.clear(); // TODO: use reserve and push_back instead
            queryKmerBuffer.init();
            matchBuffer.startIndexOfReserve = 0;
            // matchBuffer.init();

            // 1) Extract query k-mers
            time_t start = time(nullptr);
            std::cout << "Query k-mer extraction : " << std::flush;
            uint64_t seqCnt = 0;
            moreReads = kmerExtractor->extractQueryKmers(
                    queryKmerBuffer,
                    queryList,
                    seqCnt,
                    savedSeq_1,
                    savedSeq_2,
                    kseq1,
                    kseq2);
            std::cout << difftime(time(nullptr), start) << " s" << std::endl;

            // 2) Sort k-mers
            start = time(nullptr);
            std::cout << "Query k-mer sorting    : " << std::flush;
            SORT_PARALLEL(queryKmerBuffer.buffer, queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, Kmer::compareQueryKmer);
            std::cout << difftime(time(nullptr), start) << " s" << std::endl;

            // for (size_t i = 0; i < queryKmerBuffer.startIndexOfReserve; i++) {
            //     if (queryKmerBuffer.buffer[i].qInfo.sequenceID == 1003) {
            //         cout << "Sorted Kmer: ";
            //         metamerPattern->printDNA(queryKmerBuffer.buffer[i].value); cout << endl;
            //     }
            // }
            // 3) Match k-mers
            start = time(nullptr);
            bool searchComplete = kmerMatcher->matchKmers(&queryKmerBuffer, &matchBuffer, dbDir);
            if (searchComplete) {
                std::cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << std::endl;

                // 4) Sort matches
                kmerMatcher->sortMatches(&matchBuffer);

                // 5) Assign taxonomy
                assignTaxonomy(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
                for (const Query &query : queryList) {
                    if (query.classification != 0) {
                        addScoreToClade(taxonomy, query.classification, query.idScore, cladeScoreSums);
                    }
                }

                // 6) Write classification results
                start = time(nullptr);
                std::cout << "Writing results        : " << std::flush;
                reporter->writeReadClassification(queryList);
                processedReadCnt += seqCnt;                    
                std::cout << difftime(time(nullptr), start) << " s" << std::endl;
                std::cout << "Processed read count   : " << processedReadCnt << std::endl;
            } else {
                matchPerKmer *= 2;
                moreReads = true;
                std::cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << std::endl;
                break;
            }
            std::cout << "--------------------" << std::endl;

        }
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (!moreReads) {
            complete = true;
        }
        std::cout << "--------------------" << std::endl;
    }
    delete savedSeq_1;
    delete savedSeq_2;
    std::cout << "Total k-mer match count: " << kmerMatcher->getTotalMatchCnt() << std::endl;    
    reporter->closeReadClassificationFile();

    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCounts, parentToChildren);
    std::unordered_map<TaxID, double> avgScores;
    addAverageScoresToMap(cladeScoreSums, cladeCounts, avgScores);
    reporter->filterClassificationFile(
        reporter->getClassificationFileName(),
        reporter->getClassificationFileName() + ".filtered",
        avgScores,
        0.2
    );
    reporter->writeReportFile(processedReadCnt, cladeCounts, avgScores, ReportType::Default);
    std::cout << "Taxonomic classification completed." << std::endl;

}

void Classifier::classifyReadsWithPos() {
    Buffer<Kmer> queryKmerBuffer;
    Buffer<MatchWithPos> matchBuffer;
    vector<Query> queryList;
    queryList.reserve(512 * 1024);

    reporter->openReadClassificationFile();

    bool complete = false;
    SeqEntry * savedSeq_1 = new SeqEntry();
    SeqEntry * savedSeq_2 = new SeqEntry(); 
    uint64_t processedReadCnt = 0;
    std::unordered_map<TaxID, double> cladeScoreSums;

    std::cout << "--------------------" << std::endl;
    while (!complete) {
        KSeqWrapper* kseq1 = KSeqFactory(par.filenames[0].c_str());
        KSeqWrapper* kseq2 = par.seqMode == 2 ? KSeqFactory(par.filenames[1].c_str()) : nullptr;

        // Move kseq to unprocessed reads
        for (size_t i = 0; i < processedReadCnt; i++) {
            kseq1->ReadEntry();
            if (par.seqMode == 2) { kseq2->ReadEntry(); }
        }

        size_t queryKmerBufferSize = calculateBufferSize(queryList.size(), matchPerKmer, 2);
        queryKmerBuffer.reallocateMemory(queryKmerBufferSize);
        matchBuffer.reallocateMemory(queryKmerBufferSize * matchPerKmer);

        bool moreReads = true;
        while (moreReads) {
            queryList.clear(); 
            queryKmerBuffer.init();
            matchBuffer.startIndexOfReserve = 0;

            // 1) Extract query k-mers
            time_t start = time(nullptr);
            std::cout << "Query k-mer extraction : " << std::flush;
            uint64_t seqCnt = 0;
            moreReads = kmerExtractor->extractQueryKmers(
                    queryKmerBuffer,
                    queryList,
                    seqCnt,
                    savedSeq_1,
                    savedSeq_2,
                    kseq1,
                    kseq2);
            std::cout << difftime(time(nullptr), start) << " s" << std::endl;

            // 2) Sort k-mers
            start = time(nullptr);
            std::cout << "Query k-mer sorting    : " << std::flush;
            SORT_PARALLEL(queryKmerBuffer.buffer, queryKmerBuffer.buffer + queryKmerBuffer.startIndexOfReserve, Kmer::compareQueryKmer);
            std::cout << difftime(time(nullptr), start) << " s" << std::endl;

            // for (size_t i = 0; i < queryKmerBuffer.startIndexOfReserve; i++) {
            //     if (queryKmerBuffer.buffer[i].qInfo.sequenceID == 1003) {
            //         cout << "Sorted Kmer: ";
            //         metamerPattern->printDNA(queryKmerBuffer.buffer[i].value); cout << endl;
            //     }
            // }
            // 3) Match k-mers
            start = time(nullptr);
            bool searchComplete = kmerMatcher->matchKmersWithPos(&queryKmerBuffer, &matchBuffer, dbDir);
            if (searchComplete) {
                std::cout << "K-mer match count      : " << kmerMatcher->getTotalMatchCnt() << std::endl;

                // 4) Sort matches
                time_t beforeSortMatches = time(nullptr);
                std::cout << "K-mer match sorting    : " << std::flush;
                SORT_PARALLEL(matchBuffer.buffer,
                              matchBuffer.buffer + matchBuffer.startIndexOfReserve,
                              MatchWithPos::compare);
                std::cout << double(time(nullptr) - beforeSortMatches) << " s" << std::endl;

                // 5) Assign taxonomy
                assignTaxonomy<MatchWithPos>(matchBuffer.buffer, matchBuffer.startIndexOfReserve, queryList, par);
                for (const Query &query : queryList) {
                    if (query.classification != 0) {
                        addScoreToClade(taxonomy, query.classification, query.idScore, cladeScoreSums);
                    }
                }

                // 6) Write classification results
                start = time(nullptr);
                std::cout << "Writing results        : " << std::flush;
                reporter->writeReadClassification(queryList);
                processedReadCnt += seqCnt;                    
                std::cout << difftime(time(nullptr), start) << " s" << std::endl;
                std::cout << "Processed read count   : " << processedReadCnt << std::endl;
            } else {
                matchPerKmer *= 2;
                moreReads = true;
                std::cout << "--match-per-kmer was increased to " << matchPerKmer << " and searching again..." << std::endl;
                break;
            }
            std::cout << "--------------------" << std::endl;

        }
        delete kseq1;
        if (par.seqMode == 2) {
            delete kseq2;
        }
        if (!moreReads) {
            complete = true;
        }
        std::cout << "--------------------" << std::endl;
    }
    delete savedSeq_1;
    delete savedSeq_2;
    std::cout << "Total k-mer match count: " << kmerMatcher->getTotalMatchCnt() << std::endl;    
    reporter->closeReadClassificationFile();

    parseSp2GenomeSize();
    std::unordered_map<TaxID, std::vector<TaxID>> parentToChildren = taxonomy->getParentToChildren();
    unordered_map<TaxID, TaxonCounts> cladeCounts = taxonomy->getCladeCounts(taxCounts, parentToChildren);

    for (const auto& [spId, bins] : sp2coverage_global) {
        auto countIt = cladeCounts.find(spId);
        if (countIt == cladeCounts.end() || countIt->second.cladeCount == 0) {
            continue;
        }

        sp2covMetric[spId] = calCovMetrics(
            bins,
            countIt->second.cladeCount,
            sp2totalReadLength[spId],
            sp2genomeSize[spId]
        );
    }

    rollUpCoverageMetrics(parentToChildren, cladeCounts, sp2covMetric, 1);
    std::unordered_map<TaxID, double> avgScores;
    addAverageScoresToMap(cladeScoreSums, cladeCounts, avgScores);
    
    reporter->filterClassificationFile(
        reporter->getClassificationFileName(), 
        reporter->getClassificationFileName() + ".filtered", 
        sp2covMetric, 
        avgScores,
        0.4
    );
    reporter->writeReportFile(processedReadCnt, cladeCounts, sp2covMetric, avgScores, ReportType::Default);

    std::cout << "Taxonomic classification completed." << std::endl;

}

template <typename MatchType>
void Classifier::assignTaxonomy(const MatchType *matchList, 
                               size_t numOfMatches,
                               std::vector<Query> &queryList,
                               const LocalParameters &par) {
    time_t beforeAnalyze = time(nullptr);
    std::cout << "K-mer match analysis   : " << std::flush;
    // Divide matches into blocks for multi threading
    size_t seqNum = queryList.size();
    std::vector<MatchBlock> matchBlocks;
    matchBlocks.reserve(std::min(seqNum, numOfMatches));
    size_t matchIdx = 0;
    uint32_t currentQuery;
    if (par.printLog==1) {
#ifdef OPENMP
        omp_set_num_threads(1);
#endif
    }
    while (matchIdx < numOfMatches) {
        currentQuery = matchList[matchIdx].qKmer.qInfo.sequenceID;
        MatchBlock block;
        block.id = currentQuery;
        block.start = matchIdx;
        while (matchIdx < numOfMatches && currentQuery == matchList[matchIdx].qKmer.qInfo.sequenceID) {
            ++matchIdx;
        }
        block.end = matchIdx - 1;
        matchBlocks.push_back(block);
    }
    // Process each block
    if (!matchBlocks.empty()) {
#pragma omp parallel default(none), shared(matchBlocks, matchList, queryList, par)
        {
            Taxonomer<MatchType> taxonomer(par, taxonomy, metamerPattern);
            #pragma omp for schedule(guided, 8)
            for (size_t i = 0; i < matchBlocks.size(); ++i) {
                taxonomer.chooseBestTaxon(
                                matchBlocks[i].id - 1,
                                matchBlocks[i].start,
                                matchBlocks[i].end,
                                matchList,
                                queryList);
            }

            if (par.storeKmerPos) {
                for (auto& [spId, localArray] : taxonomer.sp2coverage) {
                    #pragma omp critical(sp2coverage_merge)
                    {
                        sp2totalReadLength[spId] += taxonomer.sp2totalReadLength[spId];

                        // Ensure the global map has this species
                        if (sp2coverage_global.find(spId) == sp2coverage_global.end()) {
                            sp2coverage_global[spId].resize(65536, 0);
                        }

                        uint8_t* globalBins = sp2coverage_global[spId].data();
                        uint8_t* threadBins = localArray.data();
                        saturatingAddBins(globalBins, threadBins, 65536);
                    }
                }
            }

        }
    }

    if (par.printLog) {
#ifdef OPENMP
        omp_set_num_threads(par.threads);
#endif
    }

    for (size_t i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
    }
    cout << double(time(nullptr) - beforeAnalyze) << " s" << endl;

}

void Classifier::em(size_t totalQueryCnt) {
    loadMappings(reporter->getMappingFileName());
    int maxTaxId = taxonomy->getMaxTaxID();
    vector<uint32_t> sp2uniqKmerCnt(maxTaxId + 1, 0);
    countUniqueKmerPerSpecies(sp2uniqKmerCnt);
    vector<double> sp2lengthFactor(maxTaxId + 1, 0.0);
    for (size_t i = 0; i < sp2uniqKmerCnt.size(); i++) {
        if (sp2uniqKmerCnt[i] > 0) {
            sp2lengthFactor[i] = 1.0 / log(sp2uniqKmerCnt[i]);
        } else {
            sp2lengthFactor[i] = 0.0;
        }
    }
    
    time_t st = time(nullptr);
    cout << "EM algorithm for taxonomic assignment: " << std::endl;
    size_t idx = 0;
    std::vector<std::pair<size_t, size_t>> queryRanges;
    while (idx < mappingResListSize) {
        TaxID currentQuery = mappingResList[idx].queryId;
        size_t startIdx = idx;
        while (idx < mappingResListSize && mappingResList[idx].queryId == currentQuery) {
            idx++;
        }
        queryRanges.emplace_back(startIdx, idx);
    }

    cout << "Species count: " << topSpeciesSet.size() << std::endl;
    std::unordered_map<TaxID, double> Fnew;
    std::vector<TaxID> spList;
    for (const auto & sp : topSpeciesSet) {
        taxProbs[sp] = 1.0 / topSpeciesSet.size();
        Fnew[sp] = 0.0;
        spList.push_back(sp);
    }

    std::atomic<int> converged{0};
    size_t queryCount = 0;
    #pragma omp parallel default(none) shared(queryRanges, queryCount, converged, mappingResList, taxProbs, Fnew, sp2lengthFactor, spList, std::cout)
    {
    std::unordered_map<TaxID, double> localFnew;
    size_t localQueryCount = 0;

    for (size_t iter = 0; iter < 1000; ++iter) {
        if (converged.load(std::memory_order_acquire))
            continue;

        // Init local variables
        for (auto & sp : localFnew) { sp.second = 0.0; }
        localQueryCount = 0;

        // Init shared variables
        #pragma omp single
        {
            for (auto & sp : Fnew) { sp.second = 0.0; }
            queryCount = 0;
            cout << "Iteration " << iter + 1 << ": " << std::endl;
        }

        // Calculate new probabilities for each species
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < queryRanges.size(); ++i) {
            double denom = 0.0;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                denom += mappingResList[j].score * taxProbs[mappingResList[j].speciesId] * sp2lengthFactor[mappingResList[j].speciesId];
            }
            if (denom == 0.0) { continue; }
            localQueryCount++;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                localFnew[mappingResList[j].speciesId] += (mappingResList[j].score * taxProbs[mappingResList[j].speciesId] * sp2lengthFactor[mappingResList[j].speciesId]) / denom;
            }
        }
    
        // Accumulate results from all threads
        #pragma omp critical
        {
            for (const auto & sp : localFnew) {
                Fnew[sp.first] += sp.second; 
            }
            queryCount += localQueryCount;
        }
        #pragma omp barrier

        // Update probabilities and check convergence
        #pragma omp single
        {
            for (size_t i = 0; i < spList.size(); i++) {
                Fnew[spList[i]] /= queryCount;
            }
            double delta = 0.0;
            for (size_t i = 0; i < spList.size(); i++) {
                delta += fabs(Fnew[spList[i]] - taxProbs[spList[i]]);
                if (iter > 10 && Fnew[spList[i]] < 1e-5) {
                    Fnew[spList[i]] = 0.0;
                }
            }
            cout << "Delta: " << delta << endl;
            taxProbs.swap(Fnew);    
            if (delta < 1e-6) {
                converged.fetch_add(1, std::memory_order_relaxed);
            }
        }
    }
    } // end of parallel region

    size_t explainedQueryCount = 0;
    for (size_t i = 0; i < spList.size(); i++) {
        emTaxCounts[spList[i]] = taxProbs[spList[i]] * queryCount;
        explainedQueryCount += emTaxCounts[spList[i]];
    }
    emTaxCounts[0] = totalQueryCnt - explainedQueryCount; // Unclassified queries
    cout << double(time(nullptr) - st) << " s" << endl;
    reclassify(queryRanges, mappingResList, sp2lengthFactor, totalQueryCnt);
}


void Classifier::reclassify(
    const std::vector<std::pair<size_t, size_t>> & queryRanges,
    const MappingRes * mappingResList,
    const vector<double> & sp2lengthFactor,
    size_t totalQueryCnt) 
{
    #pragma omp parallel default(none) shared(queryRanges, mappingResList, sp2lengthFactor, std::cout)
    {
        std::vector<std::pair<TaxID, double>> sp2prob;
        std::vector<TaxID> candidateSpecies;
        std::unordered_map<TaxID, unsigned int> localTaxCounts;

        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < queryRanges.size(); ++i) {
            double denom = 0.0;
            candidateSpecies.clear();
            sp2prob.clear();

            TaxID currentQuery = mappingResList[queryRanges[i].first].queryId;
            for (size_t j = queryRanges[i].first; j < queryRanges[i].second; ++j) {
                double score = taxProbs[mappingResList[j].speciesId] * mappingResList[j].score * sp2lengthFactor[mappingResList[j].speciesId];
                denom += score;
                sp2prob.emplace_back(mappingResList[j].speciesId, score);
            }
            if (denom == 0.0) {
                emResults[currentQuery].taxId = 0;
                emResults[currentQuery].score = 0.0;
                continue;
            }

            for (auto & sp : sp2prob) {
                sp.second /= denom; // Normalize scores
            }

            std::sort(sp2prob.begin(), sp2prob.end(),
                      [](const std::pair<TaxID, double> &a, const std::pair<TaxID, double> &b) {
                          return a.second > b.second;
                      });

            double sum = 0.0;
            for (size_t j = 0; j < sp2prob.size() && sum < 0.5; ++j) {
                sum += sp2prob[j].second;
                candidateSpecies.push_back(sp2prob[j].first);
            }
            emResults[currentQuery].taxId = taxonomy->LCA(candidateSpecies)->taxId;
            emResults[currentQuery].score = sum;
            localTaxCounts[emResults[currentQuery].taxId] += 1;
        }
        #pragma omp critical
        {
            for (const auto & taxCount : localTaxCounts) {
                reclassifyTaxCounts[taxCount.first] += taxCount.second;
            }
        }
    }
    size_t explainedQueryCnt = 0;
    for (const auto & taxCount : reclassifyTaxCounts) {
        if (taxCount.first != 0) { // Skip unclassified
            explainedQueryCnt += taxCount.second;
        }
    }
//    std::cout << double(time(nullptr) - st) << " s" << endl;
}

void Classifier::countUniqueKmerPerSpecies(
    vector<uint32_t> & sp2uniqKmerCnt) 
{
    string sp2uniqKmerCntFileName = dbDir + "/sp2uniqKmerCnt";
    if (FileUtil::fileExists(sp2uniqKmerCntFileName.c_str())) {
        cout << "Loading unique k-mer count per species from file: " << sp2uniqKmerCntFileName << endl;
        ifstream sp2uniqKmerCntFile(sp2uniqKmerCntFileName);
        if (!sp2uniqKmerCntFile.is_open()) {
            cerr << "Error: Could not open file " << sp2uniqKmerCntFileName << endl;
            return;
        }
        TaxID taxId;
        uint32_t count;
        while (sp2uniqKmerCntFile >> taxId >> count) {
            if (taxId < sp2uniqKmerCnt.size()) {
                sp2uniqKmerCnt[taxId] = count;
            } else {
                cerr << "Warning: TaxID " << taxId << " exceeds the size of sp2uniqKmerCnt vector." << endl;
            }
        }
        sp2uniqKmerCntFile.close();
        return;
    } else {
        string infoFileName = dbDir + "/info";
        ReadBuffer<TaxID> readBuffer(infoFileName, 1000000);
        TaxID taxId;
        time_t st = time(nullptr);
        unordered_map<TaxID, TaxID> & taxid2speciesId = kmerMatcher->getTaxId2SpeciesId();
        cout << "Count unique k-mers per species: " << std::flush;
        while((taxId = readBuffer.getNext()) != 0) {
            TaxID speciesTaxId = taxid2speciesId[taxId];
            if (speciesTaxId) {
                ++sp2uniqKmerCnt[speciesTaxId];
            }
        }
        cout << double(time(nullptr) - st) << " s" << endl;
    
        // Save the unique k-mer count per species to a file
        ofstream sp2uniqKmerCntFile(sp2uniqKmerCntFileName);
        if (!sp2uniqKmerCntFile.is_open()) {
            cerr << "Error: Could not open file " << sp2uniqKmerCntFileName << " for writing." << endl;
            return;
        }
        for (size_t i = 0; i < sp2uniqKmerCnt.size(); ++i) {
            if (sp2uniqKmerCnt[i] > 0) {
                sp2uniqKmerCntFile << i << " " << sp2uniqKmerCnt[i] << "\n";
            }
        }
        sp2uniqKmerCntFile.close();
    }
}

void Classifier::loadMappings(const string & mappingResFileName) {
    size_t fileSize = FileUtil::getFileSize(mappingResFileName);
    mappingResListSize = fileSize / sizeof(MappingRes);
    mappingResList = new MappingRes[mappingResListSize];
    FILE * mappingResFile = fopen(mappingResFileName.c_str(), "rb");
    if (mappingResFile == nullptr) {
        cerr << "Error: Could not open mapping results file " << mappingResFileName << endl;
        return;
    }
    size_t readCount = fread(mappingResList, sizeof(MappingRes), mappingResListSize, mappingResFile);
    if (readCount != mappingResListSize) {
        cerr << "Error: Could not read all mapping results from file " << mappingResFileName << endl;
        fclose(mappingResFile);
        return;
    }
    fclose(mappingResFile);
}

void Classifier::loadOriginalResults(
    const string & classificationFileName,
    size_t seqNum) 
{
    ifstream classificationFile(classificationFileName);
    if (!classificationFile.is_open()) {
        cerr << "Error: Could not open classification results file " << classificationFileName << endl;
        return;
    }
    emResults.reserve(seqNum);
    string line;
    while (getline(classificationFile, line)) {
        if (line.empty() || line[0] == '#') {
            continue; // Skip empty lines and comments
        }
        vector<string> columns = TaxonomyWrapper::splitByDelimiter(line, "\t", 4);
        if (columns.size() < 4) {
            cerr << "Error: Invalid line format in classification results file: " << line << endl;
            continue;
        }
        string queryName = columns[1];
        int length = atoi(columns[3].c_str());
        emResults.emplace_back(queryName, length);
    }
}

void Classifier::rollUpScore(
    std::unordered_map<TaxID, double> &sp2score_global,
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren,
    TaxID currentTaxID)
{
    auto childrenIt = parentToChildren.find(currentTaxID);
    if (childrenIt == parentToChildren.end() || childrenIt->second.empty()) {
        return;
    }

    double childrenSum = 0.0;
    for (TaxID childID : childrenIt->second) {
        rollUpScore(sp2score_global, parentToChildren, childID);

        auto scoreIt = sp2score_global.find(childID);
        if (scoreIt != sp2score_global.end()) {
            childrenSum += scoreIt->second;
        }
    }

    sp2score_global[currentTaxID] += childrenSum;
}

CovMetric Classifier::calCovMetrics(
    const std::vector<uint8_t>& bins,
    int readCnt,
    uint64_t totalReadLength,
    uint64_t genomeSize) 
{
    uint64_t totalCount = 0;
    double sum_c_log_c = 0.0;
    int occupiedBins = 0;

    // Macro-bin tracking
    int occupiedMacroBins = 0;
    bool macro_seen[256] = {false};

    size_t binLimit = std::min(genomeSize + 1, static_cast<uint64_t>(65536));
    for (size_t i = 1; i < binLimit; ++i) {
        uint8_t count = bins[i];
        if (count == 0) continue; // Skip empty bins

        totalCount += count;
        occupiedBins++;
        sum_c_log_c += C_LOG2_C[count];

        size_t macro_bin_idx = i >> 8; // Equivalent to i / 256
        if (!macro_seen[macro_bin_idx]) {
            macro_seen[macro_bin_idx] = true;
            occupiedMacroBins++;
        }
    }

    if (totalCount == 0) {
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }

    double effective_bins = std::min(65535.0, static_cast<double>(genomeSize));
    if (effective_bins <= 1.0) {
        return {0.0, 0.0, 0.0, 0.0, 0.0};
    }

    // Calculate Macro-bin Coverage
    double max_macro_bins = std::ceil(effective_bins / 256.0);
    max_macro_bins = std::min(max_macro_bins, 256.0);
    double macro_coverage = static_cast<double>(occupiedMacroBins) / max_macro_bins;
    macro_coverage = std::min(1.0, macro_coverage);

    // 3. Calculate Coverage (Breadth)
    double coverage = static_cast<double>(occupiedBins) / effective_bins;
    coverage = std::min(1.0, coverage); // Clamped just in case

    // 4. Calculate Observed Shannon Entropy (H_obs)
    double H_obs = std::log2(static_cast<double>(totalCount)) - (sum_c_log_c / totalCount);
    H_obs = std::max(0.0, H_obs);
    
    // 5. Calculate Standard Evenness
    // Normalized strictly against the maximum capacity of the genome
    double max_H_standard = std::log2(effective_bins);
    double evenness = (max_H_standard > 0.0001) ? std::min(1.0, H_obs / max_H_standard) : 0.0;

    // 6. Calculate Expected Occupied Bins under random uniform distribution
    double read_length = static_cast<double>(totalReadLength) / readCnt;
    double bin_size_bp = static_cast<double>(genomeSize) / effective_bins;
    double bins_per_read = 1.0 + (read_length / bin_size_bp);
    double expected_occupied = effective_bins * (1.0 - std::exp(-(readCnt * bins_per_read) / effective_bins));

    // 7. Calculate Expected Maximum Entropy (Poisson adjusted)
    double expected_max_H = std::log2(std::max(1.0, expected_occupied));

    // 8. Calculate Adjusted Evenness
    double adjustedEvenness = 0.0;
    if (expected_max_H >= 0.0001) {
        adjustedEvenness = std::min(1.0, H_obs / expected_max_H);
    }


    double unified_score = std::pow(2.0, H_obs) / effective_bins;

    // 9. Return the populated struct
    return {evenness, coverage, adjustedEvenness, unified_score, macro_coverage};
}


void Classifier::rollUpCoverageMetrics(
    const std::unordered_map<TaxID, std::vector<TaxID>>& parentToChildren,
    const std::unordered_map<TaxID, TaxonCounts>& cladeCounts,
    std::unordered_map<TaxID, CovMetric>& allMetrics, // Starts with species, gets filled with all ranks
    TaxID currentTaxID) 
{
// NEW BASE CASE: If this node already has a metric calculated (i.e., it's a Species),
    // it is the base of our roll-up. Stop recursing down to subspecies and return immediately.
    if (allMetrics.find(currentTaxID) != allMetrics.end()) {
        return;
    }

    // Secondary Base Case: True leaf node (no children at all) but no pre-calculated metric
    auto childrenIt = parentToChildren.find(currentTaxID);
    if (childrenIt == parentToChildren.end() || childrenIt->second.empty()) {
        return; 
    }

    double weightedEvenness = 0.0;
    double weightedCoverage = 0.0;
    double weightedAdjEvenness = 0.0;
    double weightedUnifiedScore = 0.0;
    double weightedMacroCoverage = 0.0;
    uint64_t totalChildReads = 0;

    // 1. Recursively process all children FIRST (Post-order traversal)
    for (TaxID childID : childrenIt->second) {
        rollUpCoverageMetrics(parentToChildren, cladeCounts, allMetrics, childID);
        
        // 2. Gather data for the weighted average
        auto countIt = cladeCounts.find(childID);
        auto metricIt = allMetrics.find(childID);

        if (countIt != cladeCounts.end() && metricIt != allMetrics.end()) {
            uint64_t reads = countIt->second.cladeCount;
            
            // Only include children that actually have metrics and reads
            if (reads > 0) {
                weightedEvenness += metricIt->second.evenness * reads;
                weightedCoverage += metricIt->second.coverage * reads;
                weightedAdjEvenness += metricIt->second.adjustedEvenness * reads;
                weightedUnifiedScore += metricIt->second.unifiedScore * reads;
                weightedMacroCoverage += metricIt->second.macroCoverage * reads; // Assuming macro coverage is the same as coverage for weighting

                totalChildReads += reads;
            }
        }
    }

    // 3. Calculate the final weighted metrics for THIS parent node
    if (totalChildReads > 0) {
        CovMetric parentMetric;
        parentMetric.evenness = weightedEvenness / static_cast<double>(totalChildReads);
        parentMetric.coverage = weightedCoverage / static_cast<double>(totalChildReads);
        parentMetric.adjustedEvenness = weightedAdjEvenness / static_cast<double>(totalChildReads);
        parentMetric.unifiedScore = weightedUnifiedScore / static_cast<double>(totalChildReads);
        parentMetric.macroCoverage = weightedMacroCoverage / static_cast<double>(totalChildReads);
        
        // Save it to the map
        allMetrics[currentTaxID] = parentMetric;
    } 
}
