#include "KmerMatcher.h"

KmerMatcher::KmerMatcher(const LocalParameters & par,
                         NcbiTaxonomy * taxonomy) {
    threads = par.threads;
    std::string dbDir = par.filenames[1 + (par.seqMode == 2)];
    targetDiffIdxFileName = dbDir + "/diffIdx";
    targetInfoFileName = dbDir + "/info";
    diffIdxSplitFileName = dbDir + "/split";

    diffIdxSplits = mmapData<DiffIdxSplit>(diffIdxSplitFileName.c_str(), 3);

    MARKER = 16777215;
    MARKER = ~ MARKER;
    hammingMargin = par.hammingMargin;
    totalMatchCnt = 0;

    // Load the taxonomy ID list
    FILE * taxIdFile;
    if((taxIdFile = fopen((dbDir + "/taxID_list").c_str(),"r")) == NULL){
        std::cout<<"Cannot open the taxID list file."<<std::endl;
        return;
    }
    char taxID[100];

    while(feof(taxIdFile) == 0) {
        fscanf(taxIdFile,"%s",taxID);
        TaxID taxId = atol(taxID);
        TaxonNode const * taxon = taxonomy->taxonNode(taxId);
        if (taxId == taxon->taxId) {
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxId2genusId[taxon->taxId] = genusTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2genusId[speciesTaxID] = genusTaxID;
        } else {
            TaxID speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
            TaxID genusTaxID = taxonomy->getTaxIdAtRank(taxId, "genus");
            while (taxon->taxId != speciesTaxID) {
                taxId2speciesId[taxon->taxId] = speciesTaxID;
                taxId2genusId[taxon->taxId] = genusTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[speciesTaxID] = speciesTaxID;
            taxId2genusId[speciesTaxID] = genusTaxID;
            taxId2speciesId[taxId] = speciesTaxID;
            taxId2genusId[taxId] = genusTaxID;
        }
    }
    fclose(taxIdFile);
}

KmerMatcher::~KmerMatcher() {
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);
}

int KmerMatcher::matchKmers(QueryKmerBuffer * queryKmerBuffer, Buffer<Match> * matchBuffer) {
    size_t queryKmerNum = queryKmerBuffer->startIndexOfReserve;
    QueryKmer *queryKmerList = queryKmerBuffer->buffer;

    size_t numOfDiffIdx = FileUtil::getFileSize(targetDiffIdxFileName) / sizeof(uint16_t);

    std::cout << "Comparing query and reference metamers..." << std::endl;

    // Find the first index of garbage query k-mer (UINT64_MAX) and discard from there
    for (size_t checkN = queryKmerNum - 1; checkN > 0; checkN--) {
        if (queryKmerList[checkN].ADkmer != UINT64_MAX) {
            queryKmerNum = checkN + 1;
            break;
        }
    }

    // Filter out meaningless target splits
    size_t numOfDiffIdxSplits = diffIdxSplits.fileSize / sizeof(DiffIdxSplit);
    size_t numOfDiffIdxSplits_use = numOfDiffIdxSplits;
    for (size_t i = 1; i < numOfDiffIdxSplits; i++) {
        if (diffIdxSplits.data[i].ADkmer == 0 || diffIdxSplits.data[i].ADkmer == UINT64_MAX) {
            diffIdxSplits.data[i] = {UINT64_MAX, UINT64_MAX, UINT64_MAX};
            numOfDiffIdxSplits_use--;
        }
    }

    // Divide query k-mer list into blocks for multi threading.
    // Each split has start and end points of query list + proper offset point of target k-mer list
    std::vector<QueryKmerSplit> querySplits;
    uint64_t queryAA;
    std::vector<int> targetSplitIdxs;

    if (threads == 1) { //Single thread
        querySplits.emplace_back(0, queryKmerNum - 1, queryKmerNum, diffIdxSplits.data[0]);
    } else if (threads == 2) { //Two threads
        size_t splitWidth = queryKmerNum / 2;
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (size_t tSplitCnt = 0; tSplitCnt < numOfDiffIdxSplits_use; tSplitCnt++) {
            queryAA = AminoAcidPart(queryKmerList[splitWidth].ADkmer);
            if (queryAA <= AminoAcidPart(diffIdxSplits.data[tSplitCnt].ADkmer)) {
                tSplitCnt = tSplitCnt - (tSplitCnt != 0);
                querySplits.emplace_back(splitWidth, queryKmerNum - 1, queryKmerNum - splitWidth,
                                         diffIdxSplits.data[tSplitCnt]);
                break;
            }
        }
    } else { //More than two threads
        // Devide query k-mers into blocks
        size_t splitWidth = queryKmerNum / (threads - 1);
        querySplits.emplace_back(0, splitWidth - 1, splitWidth, diffIdxSplits.data[0]);
        for (size_t i = 1; i < threads; i++) {
            queryAA = AminoAcidPart(queryKmerList[splitWidth * i].ADkmer);
            bool needLastTargetBlock = true;
            for (size_t j = 0; j < numOfDiffIdxSplits_use; j++) {
                if (queryAA <= AminoAcidPart(diffIdxSplits.data[j].ADkmer)) {
                    j = j - (j != 0);
                    if (i != threads - 1) {
                        querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                                 diffIdxSplits.data[j]);
                    } else {
                        querySplits.emplace_back(splitWidth * i, queryKmerNum - 1, queryKmerNum - splitWidth * i,
                                                 diffIdxSplits.data[j]);
                    }
                    targetSplitIdxs.emplace_back(j);
                    needLastTargetBlock = false;
                    break;
                }
            }
            if (needLastTargetBlock) {
                if (i != threads - 1) { // If it is not the last split
                    querySplits.emplace_back(splitWidth * i, splitWidth * (i + 1) - 1, splitWidth,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
                    targetSplitIdxs.emplace_back(numOfDiffIdxSplits_use - 2);
                } else {
                    querySplits.emplace_back(splitWidth * i, queryKmerNum - 1, queryKmerNum - splitWidth * i,
                                             diffIdxSplits.data[numOfDiffIdxSplits_use - 2]);
                    targetSplitIdxs.emplace_back(numOfDiffIdxSplits_use - 2);
                }
            }
        }
    }

    bool *splitCheckList = (bool *) malloc(sizeof(bool) * threads);
    std::fill_n(splitCheckList, threads, false);
    size_t completedSplitCnt = 0;

    time_t beforeSearch = time(nullptr);

    while (completedSplitCnt < threads) {
        bool hasOverflow = false;
#pragma omp parallel default(none), shared(completedSplitCnt, splitCheckList, hasOverflow, \
querySplits, queryKmerList, matchBuffer, cout, targetDiffIdxFileName, numOfDiffIdx, targetInfoFileName, targetSplitIdxs)
        {
            // FILE
            FILE * diffIdxFp = fopen(targetDiffIdxFileName.c_str(), "rb");
            FILE * kmerInfoFp = fopen(targetInfoFileName.c_str(), "rb");

            // Target K-mer buffer
            uint16_t * diffIdxBuffer = (uint16_t *) malloc(sizeof(uint16_t) * (BufferSize + 1)); // size = 32 Mb
            TargetKmerInfo * kmerInfoBuffer = (TargetKmerInfo *) malloc(sizeof(TargetKmerInfo) * (BufferSize+1)); // 64 Mb
            size_t kmerInfoBufferIdx = 0;
            size_t diffIdxBufferIdx = 0;

            //query variables
            uint64_t currentQuery = UINT64_MAX;
            uint64_t currentQueryAA = UINT64_MAX;
            QueryKmerInfo currentQueryInfo;

            //target variables
            size_t diffIdxPos = 0;
            std::vector<uint64_t> candidateTargetKmers; //vector for candidate target k-mer, some of which are selected after based on hamming distance
            std::vector<TargetKmerInfo> candidateKmerInfos;
            uint64_t currentTargetKmer;

            //Match buffer for each thread
            int localBufferSize = 2'000'000; // 32 Mb
            auto *matches = new Match[localBufferSize]; // 16 * 2'000'000 = 32 Mb
            int matchCnt = 0;

            // For debug
//            SeqIterator seqIterator(par);

            //vectors for selected target k-mers
            std::vector<uint8_t> selectedHammingSum;
            std::vector<size_t> selectedMatches;
            std::vector<uint16_t> selectedHammings;
            size_t posToWrite;

            int currMatchNum;
            size_t idx;
#pragma omp for schedule(dynamic, 1)
            for (size_t i = 0; i < querySplits.size(); i++) {
                if (hasOverflow || splitCheckList[i]) {
                    continue;
                }

                currentTargetKmer = querySplits[i].diffIdxSplit.ADkmer;
                diffIdxBufferIdx = querySplits[i].diffIdxSplit.diffIdxOffset;
                kmerInfoBufferIdx = querySplits[i].diffIdxSplit.infoIdxOffset
                                    - (querySplits[i].diffIdxSplit.ADkmer != 0);
                diffIdxPos = querySplits[i].diffIdxSplit.diffIdxOffset;

                fseek(kmerInfoFp, 4 * (long)(kmerInfoBufferIdx), SEEK_SET);
                loadBuffer(kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx, BufferSize);
                fseek(diffIdxFp, 2 * (long) (diffIdxBufferIdx), SEEK_SET);
                loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize);

                if (i == 0) {
                    currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                          diffIdxBufferIdx, diffIdxPos);
                }
                currentQuery = UINT64_MAX;
                currentQueryAA = UINT64_MAX;

                size_t lastMovedQueryIdx = 0;
                for (size_t j = querySplits[i].start; j < querySplits[i].end + 1; j++) {
                    querySplits[i].start++;

                    // Reuse the comparison data if queries are exactly identical
                    if (currentQuery == queryKmerList[j].ADkmer
                        && (currentQueryInfo.frame/3 == queryKmerList[j].info.frame/3)) {
                        currMatchNum = selectedMatches.size();
                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > localBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer->reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >= matchBuffer->bufferSize) {
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1;
                                __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }
                        for (int k = 0; k < currMatchNum; k++) {
                            idx = selectedMatches[k];
                            matches[matchCnt] = {queryKmerList[j].info,
                                                 candidateKmerInfos[idx].sequenceID,
                                                 taxId2genusId[candidateKmerInfos[idx].sequenceID],
                                                 taxId2speciesId[candidateKmerInfos[idx].sequenceID],
                                                 selectedHammings[k],
                                                 selectedHammingSum[k],
                                                 (bool) candidateKmerInfos[idx].redundancy};
                            matchCnt++;
                        }
                        continue;
                    }
                    selectedMatches.clear();
                    selectedHammingSum.clear();
                    selectedHammings.clear();

                    // Reuse the candidate target k-mers to compare in DNA level if queries are the same at amino acid level but not at DNA level
                    if (currentQueryAA == AminoAcidPart(queryKmerList[j].ADkmer)) {
                        compareDna(queryKmerList[j].ADkmer, candidateTargetKmers, selectedMatches,
                                   selectedHammingSum, selectedHammings,queryKmerList[j].info.frame);
                        currMatchNum = selectedMatches.size();

                        // If local buffer is full, copy them to the shared buffer.
                        if (matchCnt + currMatchNum > localBufferSize) {
                            // Check if the shared buffer is full.
                            posToWrite = matchBuffer->reserveMemory(matchCnt);
                            if (posToWrite + matchCnt >= matchBuffer->bufferSize) {
                                hasOverflow = true;
                                querySplits[i].start = lastMovedQueryIdx + 1;
                                __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                                break;
                            } else { // not full -> copy matches to the shared buffer
                                moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                                lastMovedQueryIdx = j;
                            }
                        }
                        for (int k = 0; k < currMatchNum; k++) {
                            idx = selectedMatches[k];
                            matches[matchCnt] = {queryKmerList[j].info,
                                                 candidateKmerInfos[idx].sequenceID,
                                                 taxId2genusId[candidateKmerInfos[idx].sequenceID],
                                                 taxId2speciesId[candidateKmerInfos[idx].sequenceID],
                                                 selectedHammings[k],
                                                 selectedHammingSum[k],
                                                 (bool) candidateKmerInfos[idx].redundancy};
                            matchCnt++;
                        }
                        currentQuery = queryKmerList[j].ADkmer;
                        currentQueryAA = AminoAcidPart(currentQuery);
                        currentQueryInfo = queryKmerList[j].info;
                        continue;
                    }
                    candidateTargetKmers.clear();
                    candidateKmerInfos.clear();

                    // Get next query, and start to find
                    currentQuery = queryKmerList[j].ADkmer;
                    currentQueryAA = AminoAcidPart(currentQuery);
                    currentQueryInfo = queryKmerList[j].info;

                    // Skip target k-mers that are not matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx
                           && (currentQueryAA > AminoAcidPart(currentTargetKmer))) {
                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx, BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }
                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos);
                        kmerInfoBufferIdx ++;
                    }

                    if (currentQueryAA != AminoAcidPart(currentTargetKmer)) // Move to next query k-mer if there isn't any match.
                        continue;

                    // Load target k-mers that are matched in amino acid level
                    while (diffIdxPos != numOfDiffIdx &&
                           currentQueryAA == AminoAcidPart(currentTargetKmer)) {
                        candidateTargetKmers.push_back(currentTargetKmer);
                        candidateKmerInfos.push_back(getKmerInfo(BufferSize, kmerInfoFp, kmerInfoBuffer, kmerInfoBufferIdx));
                        // Print the target k-mer
//                        if (par.printLog == 1) {
//                            cout << queryKmerList[j].info.sequenceID << "\t" << queryKmerList[j].info.pos << "\t"
//                                 << (int) queryKmerList[j].info.frame << endl;
//                            cout << "Query  k-mer: ";
//                            print_binary64(64, currentQuery);
//                            cout << "\t";
//                            seqIterator.printKmerInDNAsequence(currentQuery);
//                            cout << endl;
//                            cout << "Target k-mer: ";
//                            print_binary64(64, currentTargetKmer);
//                            cout << "\t";
//                            seqIterator.printKmerInDNAsequence(currentTargetKmer);
//                            cout << "\t" << kmerInfoBuffer[kmerInfoBufferIdx].sequenceID
//                                 << "\t" << taxId2speciesId[kmerInfoBuffer[kmerInfoBufferIdx].sequenceID] << endl;
//                            cout << (int) getHammingDistanceSum(currentQuery, currentTargetKmer) << "\t";
//                            print_binary16(16, getHammings(currentQuery, currentTargetKmer)); cout << endl;
//                        }

                        if (unlikely(BufferSize < diffIdxBufferIdx + 7)){
                            loadBuffer(diffIdxFp, diffIdxBuffer, diffIdxBufferIdx,
                                       BufferSize, ((int)(BufferSize - diffIdxBufferIdx)) * -1 );
                        }

                        currentTargetKmer = getNextTargetKmer(currentTargetKmer, diffIdxBuffer,
                                                              diffIdxBufferIdx, diffIdxPos);
                        kmerInfoBufferIdx ++;
                    }

                    // Compare the current query and the loaded target k-mers and select
                    compareDna(currentQuery, candidateTargetKmers, selectedMatches, selectedHammingSum,
                               selectedHammings, queryKmerList[j].info.frame);

                    // If local buffer is full, copy them to the shared buffer.
                    currMatchNum = selectedMatches.size();
                    if (matchCnt + currMatchNum > localBufferSize) {
                        // Check if the shared buffer is full.
                        posToWrite = matchBuffer->reserveMemory(matchCnt);
                        if (posToWrite + matchCnt >= matchBuffer->bufferSize) { // full -> write matches to file first
                            hasOverflow = true;
                            querySplits[i].start = lastMovedQueryIdx + 1;
                            __sync_fetch_and_sub(&matchBuffer->startIndexOfReserve, matchCnt);
                            break;
                        } else { // not full -> copy matches to the shared buffer
                            moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                            lastMovedQueryIdx = j;
                        }
                    }

                    for (int k = 0; k < currMatchNum; k++) {
                        idx = selectedMatches[k];
                        matches[matchCnt] = {queryKmerList[j].info,
                                             candidateKmerInfos[idx].sequenceID,
                                             taxId2genusId[candidateKmerInfos[idx].sequenceID],
                                             taxId2speciesId[candidateKmerInfos[idx].sequenceID],
                                             selectedHammings[k],
                                             selectedHammingSum[k],
                                             (bool) candidateKmerInfos[idx].redundancy};
                        matchCnt++;
                    }
                } // End of one split

                // Move matches in the local buffer to the shared buffer
                posToWrite = matchBuffer->reserveMemory(matchCnt);
                if (posToWrite + matchCnt >= matchBuffer->bufferSize) {
                    hasOverflow = true;
                    querySplits[i].start = lastMovedQueryIdx + 1;
                    __sync_fetch_and_sub(& matchBuffer->startIndexOfReserve, matchCnt);
                } else {
                    moveMatches(matchBuffer->buffer + posToWrite, matches, matchCnt);
                }

                // Check whether current split is completed or not
                if (querySplits[i].start - 1 == querySplits[i].end) {
                    splitCheckList[i] = true;
                    __sync_fetch_and_add(&completedSplitCnt, 1);
                }
            } // End of omp for (Iterating for splits)
            delete[] matches;
            fclose(diffIdxFp);
            fclose(kmerInfoFp);
            free(diffIdxBuffer);
            free(kmerInfoBuffer);
        } // End of omp parallel
        if (hasOverflow) {
            std::cout << "overflow!!!" << std::endl;
            return 2;
        }
    } // end of while(completeSplitCnt < threadNum)
    std::cout << "Time spent for the comparison: " << double(time(nullptr) - beforeSearch) << std::endl;
    munmap(diffIdxSplits.data, diffIdxSplits.fileSize + 1);
    free(splitCheckList);
    queryKmerNum = 0;

#ifdef OPENMP
    omp_set_num_threads(threads);
#endif

    // Sort matches
    time_t beforeSortMatches = time(nullptr);
    totalMatchCnt += matchBuffer->startIndexOfReserve;
    std::cout << "Sorting matches ..." << std::endl;
    SORT_PARALLEL(matchBuffer->buffer, matchBuffer->buffer + matchBuffer->startIndexOfReserve,
                  sortMatch());
    std::cout << "Time spent for sorting matches: " << double(time(nullptr) - beforeSortMatches) << std::endl;

    return 1;
}

void KmerMatcher::moveMatches(Match *dest, Match *src, int &matchNum) {
    memcpy(dest, src, sizeof(Match) * matchNum);
    matchNum = 0;
}

// It compares query k-mers to target k-mers.
// If a query has matches, the matches with the smallest hamming distance will be selected
void KmerMatcher::compareDna(uint64_t query,
                             std::vector<uint64_t> &targetKmersToCompare,
                             std::vector<size_t> &selectedMatches,
                             std::vector<uint8_t> &selectedHammingSum,
                             std::vector<uint16_t> &selectedHammings, uint8_t frame) {

    size_t size = targetKmersToCompare.size();
    auto *hammingSums = new uint8_t[size + 1];
    uint8_t currentHammingSum;
    uint8_t minHammingSum = UINT8_MAX;

    // Calculate hamming distance
    for (size_t i = 0; i < size; i++) {
        currentHammingSum = getHammingDistanceSum(query, targetKmersToCompare[i]);
        if (currentHammingSum < minHammingSum) {
            minHammingSum = currentHammingSum;
        }
        hammingSums[i] = currentHammingSum;
    }

    // Select target k-mers that passed hamming criteria
    for (size_t h = 0; h < size; h++) {
        if (hammingSums[h] <= minHammingSum + hammingMargin) {
            selectedMatches.push_back(h);
            selectedHammingSum.push_back(hammingSums[h]);
            if (frame < 3) {
                selectedHammings.push_back(getHammings(query, targetKmersToCompare[h]));
            } else {
                selectedHammings.push_back(getHammings_reverse(query, targetKmersToCompare[h]));
            }
        }
    }
    delete[] hammingSums;
}