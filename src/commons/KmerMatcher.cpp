#include "KmerMatcher.h"
#include "BitManipulateMacros.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "Mmap.h"
#include <chrono>
#include <ostream>
#include <vector>

constexpr uint16_t KmerMatcher::HAMMING_LUT0[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT1[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT2[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT3[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT4[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT5[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT6[64];
constexpr uint16_t KmerMatcher::HAMMING_LUT7[64];


KmerMatcher::KmerMatcher(
    const LocalParameters & par,
    TaxonomyWrapper * taxonomy,
    const MetamerPattern *metamerPattern) 
    : par(par), taxonomy(taxonomy), metamerPattern(metamerPattern)
{                        
    // Parameters
    threads = par.threads;
    dbDir = par.filenames[1 + (par.seqMode == 2)];
    hammingMargin = par.hammingMargin;
    DNA_MASK = ~metamerPattern->dnaMask; // metamer & DNA_MASK -> AA part 
    AA_MASK = metamerPattern->dnaMask;   // metamer & AA_MASK  -> DNA part
    totalMatchCnt = 0;
    kmerLen = metamerPattern->kmerLen;
    kmerFormat = 2;
    loadTaxIdList(par);
}

KmerMatcher::KmerMatcher(
    const LocalParameters &par,
    int kmerFormat) : par(par),
                      kmerFormat(kmerFormat)
{
    // if (par.reducedAA) {
    //     geneticCode = new ReducedGeneticCode();
    // } else {
    geneticCode = new RegularGeneticCode();
    // }
    dbDir = par.filenames[1];
    // targetDiffIdxFileName = dbDir + "/diffIdx";
    // targetInfoFileName = dbDir + "/info";
    // diffIdxSplitFileName = dbDir + "/split";
    totalMatchCnt = 0;
}

KmerMatcher::~KmerMatcher() {
    delete geneticCode;
}

void KmerMatcher::loadTaxIdList(const LocalParameters & par) {
    // cout << "Loading the list for taxonomy IDs ... " << std::flush;
    if (par.contamList != "") {
        vector<string> contams = Util::split(par.contamList, ",");
        for (auto &contam : contams) {
            string fileName = dbDir + "/" + contam + "/taxID_list";
            if (!FileUtil::fileExists(fileName.c_str())) {
                std::cout << "TaxID list file for " << contam << " does not exist." << std::endl;
                continue;
            }
            std::ifstream in{fileName};
            if (!in.is_open()) {
                std::cout << "Cannot open the taxID list file." << std::endl;
                return;
            }
            std::string line;
            while (std::getline(in, line)) {
                if (line.empty()) continue;                   
                TaxID taxId = static_cast<TaxID>(std::stoul(line));
                TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
                TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
                const TaxonNode* taxon = taxonomy->taxonNode(taxId);            
                if (taxId != taxon->taxId) {
                    taxId2speciesId[taxId] = speciesTaxID;
                    taxId2genusId[taxId] = genusTaxID;
                }
                while (taxon->taxId != speciesTaxID) {
                    taxId2speciesId[taxon->taxId] = speciesTaxID;
                    taxId2genusId[taxon->taxId] = genusTaxID;
                    taxon = taxonomy->taxonNode(taxon->parentTaxId);
                }
                taxId2speciesId[speciesTaxID] = speciesTaxID;
                taxId2genusId[speciesTaxID] = genusTaxID;
            }
            in.close();
        }
    } else {
        std::ifstream in{dbDir + "/taxID_list"};
        if (!in.is_open()) {
            std::cout << "Cannot open the taxID list file." << std::endl;
            return;
        }
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;                   
            TaxID taxId = static_cast<TaxID>(std::stoul(line));
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
            const TaxonNode* taxon = taxonomy->taxonNode(taxId);            
            if (taxId != taxon->taxId) {
                taxId2speciesId[taxId] = speciesTaxID;
                taxId2genusId[taxId] = genusTaxID;
            }
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxId2genusId[taxon->taxId] = genusTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2genusId[speciesTaxID] = genusTaxID;
        }
        in.close();
    }
    // cout << "Done" << endl;
}


bool KmerMatcher::matchKmers(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<Match> * totalMatches,
    const string & db)
{
    const string targetDiffIdxFileName = db + "/diffIdx";
    const string targetInfoFileName = db + "/info";
    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer, db);
    const Kmer * qKmers = queryKmerBuffer->buffer;

    std::vector<std::atomic<bool>> splitCheckList(querySplits.size());
    std::atomic<int> hasOverflow{0};

    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31)); 

    std::cout << "Reference k-mer search : " << flush;
    time_t beforeSearch = time(nullptr);
    #pragma omp parallel default(none), shared(splitCheckList, hasOverflow, db, \
    querySplits, qKmers, totalMatches, cout, mask)
    {
        Buffer<Match> localMatches(1024 * 1024);
        DeltaIdxReader * deltaIdxReaders = new DeltaIdxReader(db + "/diffIdx", db + "/info", 1024 * 1024 * 4, 1024 * 1024 * 4);
        std::vector<Kmer> candidates;
        std::vector<Match> filteredMatches;
        std::vector<uint8_t> hammings;
        bool localHasOverflow = false;
    #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (splitCheckList[i].exchange(true, std::memory_order_acq_rel))
                continue; 

            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            Kmer tKmer = deltaIdxReaders->next();
            Kmer qKmer(UINT64_MAX, 0);
            uint64_t qKmerAA = UINT64_MAX;
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the filtered matches if queries are exactly identical
                if ((qKmer.value == qKmers[j].value) && (qKmer.qInfo.frame/3 == qKmers[j].qInfo.frame/3)) {
                    size_t filteredMatchCnt = filteredMatches.size();
                    if (unlikely(!localMatches.afford(filteredMatchCnt))) {
                        if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                            hasOverflow.fetch_add(1, std::memory_order_relaxed);
                            localHasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(filteredMatchCnt);
                    memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                           sizeof(Match) * filteredMatches.size());
                    for (size_t k = 0; k < filteredMatchCnt; k++) {
                        localMatches.buffer[posToWrite + k].qKmer.qInfo = qKmers[j].qInfo;
                    }
                    continue;
                }
                filteredMatches.clear();

                // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                if (qKmerAA == AMINO_ACID_PART(qKmers[j].value)) {
                    filterCandidates(qKmers[j], candidates, filteredMatches, hammings);
                    size_t filteredMatchCnt = filteredMatches.size();
                    if (unlikely(!localMatches.afford(filteredMatchCnt))) {
                        if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                            hasOverflow.fetch_add(1, std::memory_order_relaxed);
                            localHasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(filteredMatchCnt);
                    memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                           sizeof(Match) * filteredMatchCnt);
                    qKmer = qKmers[j];
                    qKmerAA = AMINO_ACID_PART(qKmer.value);
                    continue;
                }
                candidates.clear();

                // Get next query, and start to find
                qKmer = qKmers[j];
                qKmerAA = AMINO_ACID_PART(qKmer.value);

                // Skip target k-mers lexiocographically smaller at amino acid level
                while (!deltaIdxReaders->isCompleted() 
                        && qKmerAA > AMINO_ACID_PART(tKmer.value)) {
                    tKmer = deltaIdxReaders->next();
                }

                // No match found - skip to the next query
                if (qKmerAA != AMINO_ACID_PART(tKmer.value)) {
                    continue;
                } 

                // Match found - load target k-mers matching at amino acid level
                while (!deltaIdxReaders->isCompleted() 
                        && qKmerAA == AMINO_ACID_PART(tKmer.value)) {
                    candidates.emplace_back(tKmer.value, tKmer.id & mask);
                    tKmer = deltaIdxReaders->next();   
                }

                filterCandidates(qKmer, candidates, filteredMatches, hammings);
                if (unlikely(!localMatches.afford(filteredMatches.size()))) {
                    if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                        hasOverflow.fetch_add(1, std::memory_order_relaxed);
                        localHasOverflow = true;
                        break;
                    }
                }
                size_t posToWrite = localMatches.reserveMemory(filteredMatches.size());
                memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                       sizeof(Match) * filteredMatches.size());
            } // End of one split


            // Move matches in the local buffer to the shared buffer
            if (!Buffer<Match>::moveSmallToLarge(&localMatches, totalMatches)) {
                hasOverflow.fetch_add(1, std::memory_order_relaxed);
                localHasOverflow = true;
            }

            // Check whether current split is completed or not
            if (!localHasOverflow) {
                splitCheckList[i] = true;
            }
        } // End of omp for (Iterating for splits)
        delete deltaIdxReaders;
    } // End of omp parallel
        
    if (hasOverflow.load(std::memory_order_acquire) > 0) {
        return false;
    }
    std::cout << double(time(nullptr) - beforeSearch) << " s" << std::endl;
    totalMatchCnt += totalMatches->startIndexOfReserve;
    return true;
}

bool KmerMatcher::matchKmersWithPos(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<MatchWithPos> * totalMatches,
    const string & db)
{
    const string targetDiffIdxFileName = db + "/diffIdx";
    const string targetInfoFileName = db + "/info";
    const string targetKmerPosFileName = db + "/kmerpos";

    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer, db);
    const Kmer * qKmers = queryKmerBuffer->buffer;

    std::vector<std::atomic<bool>> splitCheckList(querySplits.size());
    std::atomic<int> hasOverflow{0};

    std::cout << "Reference k-mer search : " << flush;
    time_t beforeSearch = time(nullptr);
    #pragma omp parallel default(none), shared(splitCheckList, hasOverflow, db, \
    querySplits, qKmers, totalMatches, cout, \
    targetDiffIdxFileName, targetInfoFileName, targetKmerPosFileName)
    {
        Buffer<MatchWithPos> localMatches(1024 * 1024);
        DeltaIdxReader * deltaIdxReaders = new DeltaIdxReader(
                targetDiffIdxFileName,
                targetInfoFileName,
                targetKmerPosFileName,
                1024 * 1024 * 4,
                1024 * 1024 * 4);
        std::vector<Kmer> candidates;
        std::vector<MatchWithPos> filteredMatches;
        std::vector<uint8_t> hammings;
        bool localHasOverflow = false;
    #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (splitCheckList[i].exchange(true, std::memory_order_acq_rel))
                continue; 

            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            Kmer tKmer = deltaIdxReaders->next();
            Kmer qKmer(UINT64_MAX, 0);
            uint64_t qKmerAA = UINT64_MAX;
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the filtered matches if queries are exactly identical
                if ((qKmer.value == qKmers[j].value) && (qKmer.qInfo.frame/3 == qKmers[j].qInfo.frame/3)) {
                    size_t filteredMatchCnt = filteredMatches.size();
                    if (unlikely(!localMatches.afford(filteredMatchCnt))) {
                        if (!Buffer<MatchWithPos>::moveSmallToLarge(&localMatches, totalMatches)) {
                            hasOverflow.fetch_add(1, std::memory_order_relaxed);
                            localHasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(filteredMatchCnt);
                    memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                           sizeof(MatchWithPos) * filteredMatches.size());
                    for (size_t k = 0; k < filteredMatchCnt; k++) {
                        localMatches.buffer[posToWrite + k].qKmer.qInfo = qKmers[j].qInfo;
                    }
                    continue;
                }
                filteredMatches.clear();

                // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                if (qKmerAA == AMINO_ACID_PART(qKmers[j].value)) {
                    filterCandidates(qKmers[j], candidates, filteredMatches, hammings);
                    size_t filteredMatchCnt = filteredMatches.size();
                    if (unlikely(!localMatches.afford(filteredMatchCnt))) {
                        if (!Buffer<MatchWithPos>::moveSmallToLarge(&localMatches, totalMatches)) {
                            hasOverflow.fetch_add(1, std::memory_order_relaxed);
                            localHasOverflow = true;
                            break;
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(filteredMatchCnt);
                    memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                           sizeof(MatchWithPos) * filteredMatchCnt);
                    qKmer = qKmers[j];
                    qKmerAA = AMINO_ACID_PART(qKmer.value);
                    continue;
                }
                candidates.clear();

                // Get next query, and start to find
                qKmer = qKmers[j];
                qKmerAA = AMINO_ACID_PART(qKmer.value);

                // Skip target k-mers lexiocographically smaller at amino acid level
                while (!deltaIdxReaders->isCompleted() 
                        && qKmerAA > AMINO_ACID_PART(tKmer.value)) {
                    tKmer = deltaIdxReaders->next();
                }

                // No match found - skip to the next query
                if (qKmerAA != AMINO_ACID_PART(tKmer.value)) {
                    continue;
                } 

                // Match found - load target k-mers matching at amino acid level
                while (!deltaIdxReaders->isCompleted() 
                        && qKmerAA == AMINO_ACID_PART(tKmer.value)) {
                    candidates.push_back(tKmer);
                    tKmer = deltaIdxReaders->next();   
                }

                filterCandidates(qKmer, candidates, filteredMatches, hammings);
                if (unlikely(!localMatches.afford(filteredMatches.size()))) {
                    if (!Buffer<MatchWithPos>::moveSmallToLarge(&localMatches, totalMatches)) {
                        hasOverflow.fetch_add(1, std::memory_order_relaxed);
                        localHasOverflow = true;
                        break;
                    }
                }
                size_t posToWrite = localMatches.reserveMemory(filteredMatches.size());
                memcpy(localMatches.buffer + posToWrite, filteredMatches.data(),
                       sizeof(MatchWithPos) * filteredMatches.size());
            } // End of one split


            // Move matches in the local buffer to the shared buffer
            if (!Buffer<MatchWithPos>::moveSmallToLarge(&localMatches, totalMatches)) {
                hasOverflow.fetch_add(1, std::memory_order_relaxed);
                localHasOverflow = true;
            }

            // Check whether current split is completed or not
            if (!localHasOverflow) {
                splitCheckList[i] = true;
            }
        } // End of omp for (Iterating for splits)
        delete deltaIdxReaders;
    } // End of omp parallel
        
    if (hasOverflow.load(std::memory_order_acquire) > 0) {
        return false;
    }
    std::cout << double(time(nullptr) - beforeSearch) << " s" << std::endl;
    totalMatchCnt += totalMatches->startIndexOfReserve;
    return true;
}

std::vector<QueryKmerSplit> KmerMatcher::makeQueryKmerSplits(
    const Buffer<Kmer> * qKmers,
    const string & dbDir) 
{
    auto beforeSplit = std::chrono::steady_clock::now();
    auto printSplitTime = [&beforeSplit]() {
        auto afterSplit = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = afterSplit - beforeSplit;
        std::cout << "Query k-mer splitting  : " << elapsed.count() << " s" << std::endl;
    };

    size_t blankCnt = std::find_if(qKmers->buffer,
                                   qKmers->buffer + qKmers->startIndexOfReserve, 
                                   [](const auto& kmer) { return kmer.id != 0;}
                                  ) - qKmers->buffer;
    size_t queryKmerNum = qKmers->startIndexOfReserve - blankCnt;
    std::cout << "Query k-mer number     : " << queryKmerNum << endl;

    // Filter out meaningless target splits
    string diffIdxSplitFileName = dbDir + "/split";
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);

    std::vector<DiffIdxSplit> targetSplits;
    targetSplits.reserve(numOfDiffIdxSplits);
    for (size_t i = 0; i < numOfDiffIdxSplits; i++) {
        if (i == 0 || (diffIdxSplits.data[i].ADkmer != 0 && diffIdxSplits.data[i].ADkmer != UINT64_MAX)) {
            targetSplits.push_back(diffIdxSplits.data[i]);
        }
    }

    std::vector<QueryKmerSplit> querySplits;
    if (queryKmerNum == 0 || targetSplits.empty()) {
        munmap(diffIdxSplits.data, diffIdxSplits.fileSize);
        printSplitTime();
        return querySplits;
    }

    auto targetStartInfo = [](const DiffIdxSplit & split) {
        return split.infoIdxOffset - (split.ADkmer != 0);
    };

    struct QueryTargetBin {
        size_t start;
        size_t end;
        size_t targetIdx;
        size_t startInfo;
        size_t endInfo;
    };

    const size_t totalTargetKmerNum = FileUtil::getFileSize((dbDir + "/info").c_str()) / sizeof(TaxID);
    std::vector<QueryTargetBin> queryBins;
    queryBins.reserve(targetSplits.size());

    auto firstQueryAfterAA = [&](uint64_t aa) {
        size_t left = blankCnt;
        size_t right = qKmers->startIndexOfReserve;
        while (left < right) {
            const size_t mid = left + (right - left) / 2;
            if (AminoAcidPart(qKmers->buffer[mid].value) <= aa) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }
        return left;
    };

    for (size_t splitIdx = 0; splitIdx < targetSplits.size(); splitIdx++) {
        const size_t queryStart = (splitIdx == 0)
                                ? blankCnt
                                : firstQueryAfterAA(AminoAcidPart(targetSplits[splitIdx].ADkmer));
        const size_t queryEnd = (splitIdx + 1 < targetSplits.size())
                              ? firstQueryAfterAA(AminoAcidPart(targetSplits[splitIdx + 1].ADkmer))
                              : qKmers->startIndexOfReserve;
        if (queryStart >= queryEnd) {
            continue;
        }

        const size_t startInfo = targetStartInfo(targetSplits[splitIdx]);
        const size_t endInfo = (splitIdx + 1 < targetSplits.size())
                             ? targetStartInfo(targetSplits[splitIdx + 1])
                             : totalTargetKmerNum;
        queryBins.push_back({queryStart, queryEnd - 1, splitIdx, startInfo, endInfo});
    }

    // Divide query k-mer list into blocks with similar estimated DB scan work.
    if (queryBins.empty()) {
        munmap(diffIdxSplits.data, diffIdxSplits.fileSize);
        printSplitTime();
        return querySplits;
    }

    size_t totalEstimatedWork = 0;
    for (const auto & bin : queryBins) {
        totalEstimatedWork += std::max<size_t>(bin.endInfo - bin.startInfo, 1);
    }

    const size_t targetTaskCnt = std::max<size_t>(
        1, std::min(queryBins.size(), static_cast<size_t>(std::max(par.threads, 1)) * 8));
    const size_t workPerTask = std::max<size_t>(
        1, (totalEstimatedWork + targetTaskCnt - 1) / targetTaskCnt);

    size_t binStart = 0;
    size_t chunkStartInfo = queryBins.front().startInfo;
    size_t chunkEndInfo = queryBins.front().endInfo;
    for (size_t i = 0; i < queryBins.size(); i++) {
        const size_t newChunkEndInfo = std::max(chunkEndInfo, queryBins[i].endInfo);
        const size_t newChunkWork = std::max<size_t>(newChunkEndInfo - chunkStartInfo, 1);

        if (i > binStart && querySplits.size() + 1 < targetTaskCnt && newChunkWork > workPerTask) {
            const auto & firstBin = queryBins[binStart];
            const auto & lastBin = queryBins[i - 1];
            querySplits.emplace_back(firstBin.start, lastBin.end, targetSplits[firstBin.targetIdx]);

            binStart = i;
            chunkStartInfo = queryBins[i].startInfo;
            chunkEndInfo = queryBins[i].endInfo;
        } else {
            chunkEndInfo = newChunkEndInfo;
        }
    }

    const auto & firstBin = queryBins[binStart];
    const auto & lastBin = queryBins.back();
    querySplits.emplace_back(firstBin.start, lastBin.end, targetSplits[firstBin.targetIdx]);

    // for (size_t i = 0; i < querySplits.size(); i++) {
    //     std::cout << "Split " << i << ": Query k-mer index " << querySplits[i].start << " - " << querySplits[i].end
    //               << ", Target k-mer offset " << querySplits[i].diffIdxSplit.ADkmer
    //               << " (diffIdx: " << querySplits[i].diffIdxSplit.diffIdxOffset
    //               << ", infoIdx: " << querySplits[i].diffIdxSplit.infoIdxOffset << ")" << std::endl;
    // }

    munmap(diffIdxSplits.data, diffIdxSplits.fileSize);
    printSplitTime();
    return querySplits;
}

bool KmerMatcher::matchKmers_AA(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<Match_AA> * totalMatches,
    const string & db)
{
    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer, db);
    const Kmer * qKmers = queryKmerBuffer->buffer;
    std::atomic<int> totalOverFlowCnt{0};
    #pragma omp parallel default(none), shared(totalOverFlowCnt, db, querySplits, qKmers, totalMatches, cout)
    {
        Buffer<Match_AA> localMatches(1024 * 1024 * 2);  // 16 Mb
        DeltaIdxReader * deltaIdxReaders 
            = new DeltaIdxReader(db + "/diffIdx",
                                 db + "/info", 
                                 1024 * 1024, 1024 * 1024);
        std::vector<Match_AA> tempMatches;  
    
        #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            if (totalOverFlowCnt.load() > 0) { continue; }
            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            Kmer tKmer = deltaIdxReaders->next();
            Kmer qKmer(UINT64_MAX, 0);
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                // Reuse the AA matches if queries are identical
                if (qKmer.value == qKmers[j].value) {
                    if (unlikely(!localMatches.afford(tempMatches.size()))) {
                        #pragma omp critical
                        {
                        if (!Buffer<Match_AA>::moveSmallToLarge(&localMatches, totalMatches)) {
                            totalOverFlowCnt++;
                            // break;
                        }
                        }
                    }
                    size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                    memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                           sizeof(Match_AA) * tempMatches.size());
                    for (size_t k = 0; k < tempMatches.size(); k++) {
                        localMatches.buffer[posToWrite + k].queryId = qKmers[j].qInfo.sequenceID;
                        localMatches.buffer[posToWrite + k].pos = qKmers[j].qInfo.pos;
                    }
                    continue;
                }
                tempMatches.clear();
                // Get next query, and start to find
                qKmer = qKmers[j];

                // Skip target k-mers lexiocographically smaller
                while (!deltaIdxReaders->isCompleted() && qKmer.value > tKmer.value) {
                    tKmer = deltaIdxReaders->next();
                }

                // No match found - skip to the next query
                if (qKmer.value != tKmer.value) { continue; } 

                // Match found - load target k-mers matching at amino acid level
                while (!deltaIdxReaders->isCompleted() && qKmer.value == tKmer.value) {
                    tempMatches.emplace_back((uint32_t) qKmer.qInfo.sequenceID, tKmer.id, (uint32_t) qKmer.qInfo.pos, tKmer.value);
                    tKmer = deltaIdxReaders->next();                                      
                }

                if (unlikely(!localMatches.afford(tempMatches.size()))) {
                                            #pragma omp critical
                        {
                    if (!Buffer<Match_AA>::moveSmallToLarge(&localMatches, totalMatches)) {
                        totalOverFlowCnt++;
                        // break;
                    }
                }
                }

                size_t posToWrite = localMatches.reserveMemory(tempMatches.size());
                memcpy(localMatches.buffer + posToWrite, tempMatches.data(),
                       sizeof(Match_AA) * tempMatches.size());
            } // End of one split

            // Move matches in the local buffer to the shared buffer
                                    #pragma omp critical
                        {
            if (!Buffer<Match_AA>::moveSmallToLarge(&localMatches, totalMatches)) {
                totalOverFlowCnt++;
            }
            }
        } // End of omp for (Iterating for splits)
    } // End of omp parallel
        
    if (totalOverFlowCnt.load() > 0) {
        return false;
    }
    this->totalMatchCnt += totalMatches->startIndexOfReserve;
    return true;
}

void KmerMatcher::sortMatches(Buffer<Match> * matchBuffer) {
    time_t beforeSortMatches = time(nullptr);
    std::cout << "K-mer match sorting    : " << std::flush;
    SORT_PARALLEL(matchBuffer->buffer,
                  matchBuffer->buffer + matchBuffer->startIndexOfReserve,
                  Match::compare);
    std::cout << double(time(nullptr) - beforeSortMatches) << " s" << std::endl;
}


void KmerMatcher::filterCandidates(
    Kmer qKmer,
    const std::vector<Kmer> &candidates,
    std::vector<Match> &filteredMatches,
    std::vector<uint8_t> & hammings
) {
    const size_t numCandidates = candidates.size();
    if (numCandidates == 0) return;

    hammings.clear();
    if (hammings.capacity() < numCandidates) {
        hammings.reserve(numCandidates);
    }

    const uint64_t qVal = qKmer.value;
    const Kmer* candPtr = candidates.data();

    if (par.useAllMatches) {
        for (size_t i = 0; i < numCandidates; i++) {
            filteredMatches.emplace_back(qKmer, candPtr[i]);
            filteredMatches.back().tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
        }
        return;
    }

    // Use only exact matches if any
    bool hasExactMatch = false;
    for (size_t i = 0; i < numCandidates; i++) {
        if (qVal == candPtr[i].value) {
            filteredMatches.emplace_back(qKmer, candPtr[i]);
            filteredMatches.back().tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
            hasExactMatch = true;
        } else if (hasExactMatch) {
            break;
        }
    }
    if (hasExactMatch) {
        return;
    }

    // Calculate hamming distances
    uint8_t minDist = UINT8_MAX;    
    for (size_t i = 0; i < numCandidates; i++) {
        uint8_t dist = metamerPattern->hammingDistSum(qVal, candPtr[i].value);
        hammings.push_back(dist);
        if (dist < minDist) {
            minDist = dist;
        }
    }

    uint8_t hDistCutoff = static_cast<uint8_t>(min(minDist * 2, kmerLen - 1));
    const uint8_t* hamPtr = hammings.data();
    for (size_t h = 0; h < numCandidates; h++) {
        if (hamPtr[h] <= hDistCutoff) {
            filteredMatches.emplace_back(qKmer, candPtr[h]);
            filteredMatches.back().tKmer.tInfo.speciesId = taxId2speciesId[candPtr[h].id];
        }
    }
}


void KmerMatcher::filterCandidates(
    Kmer qKmer,
    const std::vector<Kmer> &candidates,
    std::vector<MatchWithPos> &filteredMatches,
    std::vector<uint8_t> & hammings
) {
    const size_t numCandidates = candidates.size();
    if (numCandidates == 0) return;

    hammings.clear();
    if (hammings.capacity() < numCandidates) {
        hammings.reserve(numCandidates);
    }

    const uint64_t qVal = qKmer.value;
    const Kmer* candPtr = candidates.data();

    if (par.useAllMatches) {
        for (size_t i = 0; i < numCandidates; i++) {
            filteredMatches.emplace_back(qKmer, candPtr[i], static_cast<uint16_t>(candPtr[i].tInfo.pos));
            filteredMatches.back().tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
        }
        return;
    }

    // Use only exact matches if any
    bool hasExactMatch = false;
    for (size_t i = 0; i < numCandidates; i++) {
        if (qVal == candPtr[i].value) {
            filteredMatches.emplace_back(qKmer, candPtr[i], static_cast<uint16_t>(candPtr[i].tInfo.pos));
            filteredMatches.back().tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
            hasExactMatch = true;
        } else if (hasExactMatch) {
            break;
        }
    }
    if (hasExactMatch) {
        return;
    }

    // Calculate hamming distances
    uint8_t minDist = UINT8_MAX;    
    for (size_t i = 0; i < numCandidates; i++) {
        uint8_t dist = metamerPattern->hammingDistSum(qVal, candPtr[i].value);
        hammings.push_back(dist);
        if (dist < minDist) {
            minDist = dist;
        }
    }

    uint8_t hDistCutoff = static_cast<uint8_t>(min(minDist * 2, kmerLen - 1));
    const uint8_t* hamPtr = hammings.data();
    for (size_t h = 0; h < numCandidates; h++) {
        if (hamPtr[h] <= hDistCutoff) {
            filteredMatches.emplace_back(qKmer, candPtr[h], static_cast<uint16_t>(candPtr[h].tInfo.pos));
            filteredMatches.back().tKmer.tInfo.speciesId = taxId2speciesId[candPtr[h].id];
        }
    }
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void KmerMatcher::compareDna(uint64_t query,
                             std::vector<uint64_t> &targetKmersToCompare,
                             std::vector<uint8_t> & hammingDists,
                             std::vector<size_t> &selectedMatches,
                             std::vector<uint8_t> &selectedHammingSum,
                             std::vector<uint16_t> &selectedHammings,
                             size_t & selectedMatchIdx,
                             uint8_t frame) {
    hammingDists.resize(targetKmersToCompare.size());
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < targetKmersToCompare.size(); i++) {
        hammingDists[i] = getHammingDistanceSum(query, targetKmersToCompare[i]);
        minHammingSum = min(minHammingSum, hammingDists[i]);
    }

    // Select target k-mers that passed hamming criteria
    selectedMatchIdx = 0;
    uint8_t maxHamming = min(minHammingSum * 2, 7);
    for (size_t h = 0; h < targetKmersToCompare.size(); h++) {
        if (hammingDists[h] <= maxHamming) {
            selectedHammingSum[selectedMatchIdx] = hammingDists[h];
            selectedHammings[selectedMatchIdx] = !((frame < 3) ^ (kmerFormat == 2))
                ? getHammings(query, targetKmersToCompare[h])
                : getHammings_reverse(query, targetKmersToCompare[h]);
            selectedMatches[selectedMatchIdx++] = h;
        }
    }
}
