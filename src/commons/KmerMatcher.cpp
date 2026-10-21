#include "KmerMatcher.h"
#include "BitManipulateMacros.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "Mmap.h"
#include <algorithm>
#include <ostream>
#include <vector>

// DeltaIdxReader decodes the delta stream into valueBuffer (16-byte Kmers),
// which the join loop then reads back one at a time. When the database is
// resident in RAM (the common server case), a buffer larger than the private
// L2 forces that write->read round-trip through DRAM, and the search is
// DRAM-bandwidth-bound. Sizing both buffers to stay within one core's L2
// (512 KB on EPYC Milan; similar elsewhere) keeps the round-trip in cache.
// Candidate groups are NOT limited by this: they accumulate in a separate
// std::vector and next() refills transparently across buffer boundaries.
static constexpr size_t SEARCH_VALUE_BUFFER = 1024 * 16;  // 16K Kmers = 256 KB
static constexpr size_t SEARCH_DELTA_BUFFER = 1024 * 64;  // 64K uint16 = 128 KB

// METABULI_IO_PROF=1: per-search breakdown of time spent blocked in freads
// (summed over threads) vs the parallel search wall time, plus how many
// chunks the adaptive read-ahead decided to fetch in the background.
static void printIoProf(double wallSec, int threads,
                        uint64_t ioNanos0, uint64_t ioBytes0, uint64_t pfChunks0) {
    if (!ioProfEnabled()) {
        return;
    }
    const double waitSec = (ioProfNanos().load() - ioNanos0) / 1e9;
    const double gib = (ioProfBytes().load() - ioBytes0) / 1073741824.0;
    const double share = (wallSec > 0 && threads > 0)
                             ? waitSec / (wallSec * threads) * 100.0 : 0.0;
    std::cout << "I/O wait (thread sum)  : " << waitSec << " s, " << gib
              << " GiB read, " << share << "% of search thread-time, "
              << (ioPrefetchChunks().load() - pfChunks0) << " prefetched chunks"
              << std::endl;
}

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
    const MetamerPattern *metamerPattern,
    const std::vector<std::string> & taxIdListDirs)
    : par(par), taxonomy(taxonomy), metamerPattern(metamerPattern)
{
    // Parameters
    threads = par.threads;
    dbDir = par.filenames[1 + par.pairedFileInput()];
    hammingMargin = par.hammingMargin;
    DNA_MASK = ~metamerPattern->dnaMask; // metamer & DNA_MASK -> AA part
    AA_MASK = metamerPattern->dnaMask;   // metamer & AA_MASK  -> DNA part
    totalMatchCnt = 0;
    kmerLen = metamerPattern->kmerLen;
    kmerFormat = 2;
    if (!taxIdListDirs.empty() && par.contamList == "") {
        loadTaxIdLists(taxIdListDirs);
    } else {
        loadTaxIdList(par);
    }
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

void KmerMatcher::addTaxIdList(const std::string & fileName) {
    std::ifstream in{fileName};
    if (!in.is_open()) {
        std::cout << "Cannot open the taxID list file " << fileName << std::endl;
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
            addTaxIdList(fileName);
        }
    } else {
        addTaxIdList(dbDir + "/taxID_list");
    }
    // cout << "Done" << endl;
}

void KmerMatcher::loadTaxIdLists(const std::vector<std::string> & dbDirs) {
    for (const std::string & dir : dbDirs) {
        addTaxIdList(dir + "/taxID_list");
    }
}


bool KmerMatcher::matchKmers(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<Match> * totalMatches,
    const string & db)
{
    const string targetDiffIdxFileName = db + "/diffIdx";
    const string targetInfoFileName = db + "/info";
    // The caller may search several sub-DBs into the same buffer; count only
    // the matches this call appends.
    const size_t matchCntBefore = totalMatches->startIndexOfReserve;
    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer, db);
    const Kmer * qKmers = queryKmerBuffer->buffer;

    std::vector<std::atomic<bool>> splitCheckList(querySplits.size());
    std::atomic<int> hasOverflow{0};

    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31)); 

    if (par.lazyInfo == 0) {
        std::cout << "Note: --lazy-info 0 is deprecated and ignored; "
                     "lazy info fetching is always used." << std::endl;
    }
    std::cout << "Reference k-mer search : " << flush;
    time_t beforeSearch = time(nullptr);
    const uint64_t ioNanos0 = ioProfEnabled() ? ioProfNanos().load() : 0;
    const uint64_t ioBytes0 = ioProfEnabled() ? ioProfBytes().load() : 0;
    const uint64_t pfChunks0 = ioProfEnabled() ? ioPrefetchChunks().load() : 0;
    #pragma omp parallel default(none), shared(splitCheckList, hasOverflow, db, \
    querySplits, qKmers, totalMatches, cout, mask)
    {
        Buffer<Match> localMatches(1024 * 1024);
        // Lazy info: tax IDs are fetched per AA-matched candidate group
        // instead of being streamed alongside the metamers, so info-file
        // regions that never match a query are never read.
        DeltaIdxReader * deltaIdxReaders = new DeltaIdxReader(
            db + "/diffIdx",
            db + "/info",
            SEARCH_VALUE_BUFFER,
            SEARCH_DELTA_BUFFER,
            true /* lazyInfo */,
            true /* deferInitialFill: setReadPosition() seeks before any read */,
            true /* prefetch-eligible (active with METABULI_PREFETCH=1) */);
        std::vector<Kmer> candidates;
        std::vector<Match> filteredMatches;
        std::vector<uint8_t> hammings;
        std::vector<uint32_t> idScratch;
        candidates.reserve(1024);
        filteredMatches.reserve(1024);
        hammings.reserve(1024);
        idScratch.reserve(1024);
        size_t groupStartOrdinal = 0;
        bool localHasOverflow = false;
    #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;

            if (splitCheckList[i].exchange(true, std::memory_order_acq_rel))
                continue;

            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            // Pointer-based scan over the decoded window of bare metamer
            // values; tax IDs arrive later via fetchInfos().
            const uint64_t * tCur;
            const uint64_t * tEnd;
            deltaIdxReaders->currentWindow(tCur, tEnd);
            bool streamDone = (tCur == tEnd) && !deltaIdxReaders->nextWindow(tCur, tEnd);
            Kmer qKmer(UINT64_MAX, 0);
            uint64_t qKmerAA = UINT64_MAX;
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                const Kmer& currentQKmer = qKmers[j];
                const uint64_t currentQKmerAA = AminoAcidPart(currentQKmer.value);
                // Reuse the filtered matches if queries are exactly identical
                if ((qKmer.value == currentQKmer.value) && (qKmer.qInfo.frame/3 == currentQKmer.qInfo.frame/3)) {
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
                        localMatches.buffer[posToWrite + k].qKmer.qInfo = currentQKmer.qInfo;
                    }
                    continue;
                }
                filteredMatches.clear();

                // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                if (qKmerAA == currentQKmerAA) {
                    filterCandidates(currentQKmer, candidates, filteredMatches, hammings);
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
                    qKmer = currentQKmer;
                    qKmerAA = currentQKmerAA;
                    continue;
                }
                candidates.clear();

                // Get next query, and start to find
                qKmer = currentQKmer;
                qKmerAA = currentQKmerAA;

                // Skip target k-mers lexicographically smaller at amino acid
                // level. AA bits are the high bits and qKmerAA has the DNA
                // bits cleared, so AA(v) >= qKmerAA <=> v >= qKmerAA and the
                // sorted window can be searched directly: a short linear probe
                // for the common small gap, then a binary search for long ones.
                while (!streamDone) {
                    const uint64_t * probeEnd = (tEnd - tCur > 32) ? tCur + 32 : tEnd;
                    while (tCur < probeEnd && *tCur < qKmerAA) {
                        ++tCur;
                    }
                    if (tCur < tEnd && *tCur < qKmerAA) {
                        tCur = std::lower_bound(tCur, tEnd, qKmerAA);
                    }
                    if (tCur < tEnd) {
                        break;
                    }
                    streamDone = !deltaIdxReaders->nextWindow(tCur, tEnd);
                }

                // No match found - skip to the next query
                if (streamDone || AminoAcidPart(*tCur) != qKmerAA) {
                    continue;
                }

                // Match found - load target k-mers matching at amino acid
                // level. The group is contiguous but may span window refills.
                groupStartOrdinal = deltaIdxReaders->ordinalOf(tCur);
                while (!streamDone) {
                    while (tCur < tEnd && AminoAcidPart(*tCur) == qKmerAA) {
                        candidates.emplace_back(*tCur, 0u);
                        ++tCur;
                    }
                    if (tCur < tEnd) {
                        break;
                    }
                    streamDone = !deltaIdxReaders->nextWindow(tCur, tEnd);
                }

                // Candidates are consecutive k-mers of the sorted list; fetch
                // their tax IDs as one contiguous range.
                idScratch.resize(candidates.size());
                deltaIdxReaders->fetchInfos(groupStartOrdinal, candidates.size(), idScratch.data());
                for (size_t c = 0; c < candidates.size(); ++c) {
                    candidates[c].id = idScratch[c] & mask;
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
    printIoProf(double(time(nullptr) - beforeSearch), par.threads, ioNanos0, ioBytes0, pfChunks0);
    lastMatchCnt = totalMatches->startIndexOfReserve - matchCntBefore;
    totalMatchCnt += lastMatchCnt;
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
    const size_t matchCntBefore = totalMatches->startIndexOfReserve;

    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer, db);
    const Kmer * qKmers = queryKmerBuffer->buffer;

    std::vector<std::atomic<bool>> splitCheckList(querySplits.size());
    std::atomic<int> hasOverflow{0};

    std::cout << "Reference k-mer search : " << flush;
    time_t beforeSearch = time(nullptr);
    const uint64_t ioNanos0 = ioProfEnabled() ? ioProfNanos().load() : 0;
    const uint64_t ioBytes0 = ioProfEnabled() ? ioProfBytes().load() : 0;
    const uint64_t pfChunks0 = ioProfEnabled() ? ioPrefetchChunks().load() : 0;
    #pragma omp parallel default(none), shared(splitCheckList, hasOverflow, db, \
    querySplits, qKmers, totalMatches, cout, \
    targetDiffIdxFileName, targetInfoFileName, targetKmerPosFileName)
    {
        Buffer<MatchWithPos> localMatches(1024 * 1024);
        DeltaIdxReader * deltaIdxReaders = new DeltaIdxReader(
                targetDiffIdxFileName,
                targetInfoFileName,
                targetKmerPosFileName,
                SEARCH_VALUE_BUFFER,
                SEARCH_DELTA_BUFFER,
                true /* deferInitialFill: setReadPosition() seeks before any read */,
                true /* prefetch-eligible (active with METABULI_PREFETCH=1) */);
        std::vector<Kmer> candidates;
        std::vector<MatchWithPos> filteredMatches;
        std::vector<uint8_t> hammings;
        candidates.reserve(1024);
        filteredMatches.reserve(1024);
        hammings.reserve(1024);
        bool localHasOverflow = false;
    #pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < querySplits.size(); i++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (splitCheckList[i].exchange(true, std::memory_order_acq_rel))
                continue; 

            deltaIdxReaders->setReadPosition(querySplits[i].diffIdxSplit);
            Kmer tKmer = deltaIdxReaders->next();
            uint64_t tKmerAA = AminoAcidPart(tKmer.value);
            Kmer qKmer(UINT64_MAX, 0);
            uint64_t qKmerAA = UINT64_MAX;
            for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                const Kmer& currentQKmer = qKmers[j];
                const uint64_t currentQKmerAA = AminoAcidPart(currentQKmer.value);
                // Reuse the filtered matches if queries are exactly identical
                if ((qKmer.value == currentQKmer.value) && (qKmer.qInfo.frame/3 == currentQKmer.qInfo.frame/3)) {
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
                        localMatches.buffer[posToWrite + k].qKmer.qInfo = currentQKmer.qInfo;
                    }
                    continue;
                }
                filteredMatches.clear();

                // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                if (qKmerAA == currentQKmerAA) {
                    filterCandidates(currentQKmer, candidates, filteredMatches, hammings);
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
                    qKmer = currentQKmer;
                    qKmerAA = currentQKmerAA;
                    continue;
                }
                candidates.clear();

                // Get next query, and start to find
                qKmer = currentQKmer;
                qKmerAA = currentQKmerAA;

                // Skip target k-mers lexiocographically smaller at amino acid level
                while (!deltaIdxReaders->isCompleted() 
                        && qKmerAA > tKmerAA) {
                    tKmer = deltaIdxReaders->next();
                    tKmerAA = AminoAcidPart(tKmer.value);
                }

                // No match found - skip to the next query
                if (deltaIdxReaders->isCompleted() || qKmerAA != tKmerAA) {
                    continue;
                } 

                // Match found - load target k-mers matching at amino acid level
                while (!deltaIdxReaders->isCompleted() 
                        && qKmerAA == tKmerAA) {
                    candidates.push_back(tKmer);
                    tKmer = deltaIdxReaders->next();
                    tKmerAA = AminoAcidPart(tKmer.value);
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
    printIoProf(double(time(nullptr) - beforeSearch), par.threads, ioNanos0, ioBytes0, pfChunks0);
    lastMatchCnt = totalMatches->startIndexOfReserve - matchCntBefore;
    totalMatchCnt += lastMatchCnt;
    return true;
}

size_t KmerMatcher::countQueryKmerBlanks(const Buffer<Kmer> * qKmers) {
    return std::find_if(qKmers->buffer,
                        qKmers->buffer + qKmers->startIndexOfReserve,
                        [](const auto& kmer) { return kmer.id != 0;}
                       ) - qKmers->buffer;
}

std::vector<QueryKmerSplit> KmerMatcher::makeQueryKmerSplits(
    const Buffer<Kmer> * qKmers,
    const string & dbDir)
{
    size_t blankCnt = countQueryKmerBlanks(qKmers);
    size_t queryKmerNum = qKmers->startIndexOfReserve - blankCnt;

    // Filter out meaningless target splits
    string diffIdxSplitFileName = dbDir + "/split";
    MmapedData<DiffIdxSplit> diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading. Oversplit
    // relative to the thread count so schedule(dynamic) can balance
    // match-dense regions and the number of concurrent I/O streams does not
    // decay to the few threads holding the slowest splits. Each extra split
    // rescans from its anchor, so the overhead is bounded by the anchor count.
    size_t numSplits = std::min(static_cast<size_t>(par.threads) * 4, numOfDiffIdxSplits_use);
    numSplits = std::max(numSplits, static_cast<size_t>(par.threads)); // at least one split per thread
    numSplits = std::max<size_t>(std::min(numSplits, queryKmerNum), 1);
    std::vector<QueryKmerSplit> querySplits;
    size_t quotient = queryKmerNum / numSplits;
    size_t remainder = queryKmerNum % numSplits;
    size_t startIdx = blankCnt;
    size_t endIdx = 0; // endIdx is inclusive
    for (size_t i = 0; i < numSplits; i++) {
        endIdx = startIdx + quotient - 1;
        if (remainder > 0) {
            endIdx++;
            remainder--;
        }
        bool needLastTargetBlock = true;
        uint64_t queryAA = AminoAcidPart(qKmers->buffer[startIdx].value);
        for (size_t j = 0; j < numOfDiffIdxSplits_use; j ++) {
            if (queryAA <= AminoAcidPart(diffIdxSplits.data[j].ADkmer)) {
                querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[j - (j != 0)]);
                needLastTargetBlock = false;
                break;
            }
        }
        if (needLastTargetBlock) {
            querySplits.emplace_back(startIdx, endIdx, diffIdxSplits.data[numOfDiffIdxSplits_use - 1]);
        }
        startIdx = endIdx + 1;
    }

    // for (size_t i = 0; i < querySplits.size(); i++) {
    //     std::cout << "Split " << i << ": Query k-mer index " << querySplits[i].start << " - " << querySplits[i].end
    //               << ", Target k-mer offset " << querySplits[i].diffIdxSplit.ADkmer
    //               << " (diffIdx: " << querySplits[i].diffIdxSplit.diffIdxOffset
    //               << ", infoIdx: " << querySplits[i].diffIdxSplit.infoIdxOffset << ")" << std::endl;
    // }

    munmap(diffIdxSplits.data, diffIdxSplits.fileSize);
    return querySplits;
}

bool KmerMatcher::matchKmers_AA(
    const Buffer<Kmer> * queryKmerBuffer,
    Buffer<Match_AA> * totalMatches,
    const string & db)
{
    const size_t matchCntBefore = totalMatches->startIndexOfReserve;
    std::vector<QueryKmerSplit> querySplits = makeQueryKmerSplits(queryKmerBuffer, db);
    const Kmer * qKmers = queryKmerBuffer->buffer;
    std::atomic<int> totalOverFlowCnt{0};
    #pragma omp parallel default(none), shared(totalOverFlowCnt, db, querySplits, qKmers, totalMatches, cout)
    {
        Buffer<Match_AA> localMatches(1024 * 1024 * 2);  // 16 Mb
        DeltaIdxReader * deltaIdxReaders
            = new DeltaIdxReader(db + "/diffIdx",
                                 db + "/info",
                                 SEARCH_VALUE_BUFFER, SEARCH_DELTA_BUFFER, false,
                                 true /* deferInitialFill: setReadPosition() seeks before any read */,
                                 true /* prefetch-eligible (active with METABULI_PREFETCH=1) */);
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
    this->lastMatchCnt = totalMatches->startIndexOfReserve - matchCntBefore;
    this->totalMatchCnt += lastMatchCnt;
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

    const uint64_t qVal = qKmer.value;
    const Kmer* candPtr = candidates.data();
    if (filteredMatches.capacity() < numCandidates) {
        filteredMatches.reserve(numCandidates);
    }

    if (par.useAllMatches) {
        for (size_t i = 0; i < numCandidates; i++) {
            Kmer tKmer = candPtr[i];
            tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
            filteredMatches.emplace_back(qKmer, tKmer);
        }
        return;
    }

    if (par.maxHdist >= 0) {
        hammings.resize(numCandidates);
        uint8_t* hamPtr = hammings.data();
        metamerPattern->hammingDistSumBatch(qVal, candPtr, numCandidates, hamPtr);
        for (size_t h = 0; h < numCandidates; h++) {
            if (static_cast<int>(hamPtr[h]) <= par.maxHdist) {
                Kmer tKmer = candPtr[h];
                tKmer.tInfo.speciesId = taxId2speciesId[candPtr[h].id];
                filteredMatches.emplace_back(qKmer, tKmer);
            }
        }
        return;
    }

    // Use only exact matches if any
    bool hasExactMatch = false;
    for (size_t i = 0; i < numCandidates; i++) {
        if (qVal == candPtr[i].value) {
            Kmer tKmer = candPtr[i];
            tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
            filteredMatches.emplace_back(qKmer, tKmer);
            hasExactMatch = true;
        } else if (hasExactMatch) {
            break;
        }
    }
    if (hasExactMatch) {
        return;
    }

    // Calculate hamming distances
    hammings.resize(numCandidates);
    uint8_t* hamPtr = hammings.data();
    metamerPattern->hammingDistSumBatch(qVal, candPtr, numCandidates, hamPtr);
    uint8_t minDist = UINT8_MAX;
    for (size_t i = 0; i < numCandidates; i++) {
        const uint8_t dist = hamPtr[i];
        if (dist < minDist) {
            minDist = dist;
        }
    }

    const uint8_t hDistCutoff = static_cast<uint8_t>(min(static_cast<int>(minDist) * 2, kmerLen - 1));
    for (size_t h = 0; h < numCandidates; h++) {
        if (hamPtr[h] <= hDistCutoff) {
            Kmer tKmer = candPtr[h];
            tKmer.tInfo.speciesId = taxId2speciesId[candPtr[h].id];
            filteredMatches.emplace_back(qKmer, tKmer);
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

    const uint64_t qVal = qKmer.value;
    const Kmer* candPtr = candidates.data();
    if (filteredMatches.capacity() < numCandidates) {
        filteredMatches.reserve(numCandidates);
    }

    if (par.useAllMatches) {
        for (size_t i = 0; i < numCandidates; i++) {
            const uint16_t posId = static_cast<uint16_t>(candPtr[i].tInfo.pos);
            Kmer tKmer = candPtr[i];
            tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
            filteredMatches.emplace_back(qKmer, tKmer, posId);
        }
        return;
    }

    if (par.maxHdist >= 0) {
        hammings.resize(numCandidates);
        uint8_t* hamPtr = hammings.data();
        metamerPattern->hammingDistSumBatch(qVal, candPtr, numCandidates, hamPtr);
        for (size_t h = 0; h < numCandidates; h++) {
            if (static_cast<int>(hamPtr[h]) <= par.maxHdist) {
                const uint16_t posId = static_cast<uint16_t>(candPtr[h].tInfo.pos);
                Kmer tKmer = candPtr[h];
                tKmer.tInfo.speciesId = taxId2speciesId[candPtr[h].id];
                filteredMatches.emplace_back(qKmer, tKmer, posId);
            }
        }
        return;
    }

    // Use only exact matches if any
    bool hasExactMatch = false;
    for (size_t i = 0; i < numCandidates; i++) {
        if (qVal == candPtr[i].value) {
            const uint16_t posId = static_cast<uint16_t>(candPtr[i].tInfo.pos);
            Kmer tKmer = candPtr[i];
            tKmer.tInfo.speciesId = taxId2speciesId[candPtr[i].id];
            filteredMatches.emplace_back(qKmer, tKmer, posId);
            hasExactMatch = true;
        } else if (hasExactMatch) {
            break;
        }
    }
    if (hasExactMatch) {
        return;
    }

    // Calculate hamming distances
    hammings.resize(numCandidates);
    uint8_t* hamPtr = hammings.data();
    metamerPattern->hammingDistSumBatch(qVal, candPtr, numCandidates, hamPtr);
    uint8_t minDist = UINT8_MAX;
    for (size_t i = 0; i < numCandidates; i++) {
        const uint8_t dist = hamPtr[i];
        if (dist < minDist) {
            minDist = dist;
        }
    }

    const uint8_t hDistCutoff = static_cast<uint8_t>(min(static_cast<int>(minDist) * 2, kmerLen - 1));
    for (size_t h = 0; h < numCandidates; h++) {
        if (hamPtr[h] <= hDistCutoff) {
            const uint16_t posId = static_cast<uint16_t>(candPtr[h].tInfo.pos);
            Kmer tKmer = candPtr[h];
            tKmer.tInfo.speciesId = taxId2speciesId[candPtr[h].id];
            filteredMatches.emplace_back(qKmer, tKmer, posId);
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
