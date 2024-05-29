#include "FuncIndexer.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "KmerMatcher.h"
#include "Match.h"
#include "ProteinDbIndexer.h"
#include <cstdint>
#include <cstring>

FuncIndexer::FuncIndexer(const LocalParameters &par) : IndexCreator(par), par(par) {
    protDBFileName = par.proteinDB + "/protein.mtbl";
    protDbSplitFileName = par.proteinDB + "/prot_split.mtbl";
    protIdMapFileName = par.proteinDB + "/prtIdMap.mtbl";
    loadProtIdMap();
    kmerMatcher = new KmerMatcher(par, taxonomy);

    queryCdsList.resize(par.bufferSize / 50);
    
}

FuncIndexer::~FuncIndexer() {
    
}

void FuncIndexer::createIndex() {
    loadCdsInfo(par.cdsInfo);
    
    // Read through FASTA files and make blocks of sequences to be processed by each thread
    if (par.accessionLevel) {
        makeBlocksForParallelProcessing_accession_level();
    } else {
        makeBlocksForParallelProcessing();
    }
    writeTaxIdList();

    // Process the splits until all are processed
    bool * tempChecker = new bool[fnaSplits.size()];
    bool * completionChecker = new bool[fnaSplits.size()];
    memset(tempChecker, 0, fnaSplits.size() * sizeof(bool));
    memset(completionChecker, 0, fnaSplits.size() * sizeof(bool));
    // fill_n(completionChecker, fnaSplits.size(), false);
    size_t totalProcessedSplitCnt = 0;
    Buffer<TargetMetamerF> metamerBuffer(par.bufferSize);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(totalProcessedSplitCnt < fnaSplits.size()) {
        queryCdsList.clear();
        queryCdsList.resize(par.bufferSize / 50);

        size_t processedSplitCnt = 0;
        memset(metamerBuffer.buffer, 0, metamerBuffer.bufferSize * sizeof(TargetMetamerF));
        fillTargetKmerBuffer(metamerBuffer, tempChecker, processedSplitCnt);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(metamerBuffer.buffer,
                      metamerBuffer.buffer + metamerBuffer.startIndexOfReserve, 
                      [] (const TargetMetamerF & a, const TargetMetamerF & b) {
                                return a.metamerF.metamer < b.metamerF.metamer;});
        time_t sort = time(nullptr);
        cout << "Sort time: " << sort - start << endl;
        
        // Search protein database to label the k-mers
        if (getProteinId(metamerBuffer)) {
            memcpy(completionChecker, tempChecker, fnaSplits.size() * sizeof(bool));
            totalProcessedSplitCnt += processedSplitCnt;
        } else {
            memcpy(tempChecker, completionChecker, fnaSplits.size() * sizeof(bool));
            metamerBuffer.reallocateMemory(metamerBuffer.bufferSize / 2);
            continue;
        }

        // Reduce redundancy



        // // Reduce redundancy
        // auto * uniqKmerIdx = new size_t[metamerBuffer.startIndexOfReserve + 1];
        // size_t uniqKmerCnt = 0;
        // reduceRedundancy(metamerBuffer, uniqKmerIdx, uniqKmerCnt, par);
        // time_t reduction = time(nullptr);
        // cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        // if(processedSplitCnt == fnaSplits.size() && numOfFlush == 0){
        //     writeTargetFilesAndSplits(metamerBuffer.buffer, metamerBuffer.startIndexOfReserve, par, uniqKmerIdx, uniqKmerCnt);
        // } else {
        //     writeTargetFiles(metamerBuffer.buffer, metamerBuffer.startIndexOfReserve, par,uniqKmerIdx, uniqKmerCnt);
        // }
        // delete[] uniqKmerIdx;
    }
    delete[] tempChecker;
    delete[] completionChecker;
    writeTaxonomyDB();
    writeDbParameters();
}


size_t FuncIndexer::fillTargetKmerBuffer(Buffer<TargetMetamerF> &kmerBuffer,
                                         bool * tempChecker,
                                         size_t &processedSplitCnt) {
    uint32_t cdsCnt = 1;
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(kmerBuffer, tempChecker, processedSplitCnt, hasOverflow, cout, cdsCnt)
    {
        ProbabilityMatrix probMatrix(*subMat);
        SeqIterator seqIterator(par);
        size_t posToWrite;
        priority_queue<uint64_t> observedNonCDSKmerHashes;
        char *reverseComplement;
        string reverseComplementStr;
        kseq_buffer_t buffer;
        kseq_t *seq;
        vector<uint64_t> observedNonCDSKmers;
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < fnaSplits.size(); i++) {
            if (!tempChecker[i] && !hasOverflow) {
                tempChecker[i] = true;
                observedNonCDSKmers.clear();

                // Estimate the number of k-mers to be extracted from current split
                size_t totalLength = 0;
                for (size_t p = 0; p < fnaSplits[i].cnt; p++) {
                    totalLength += fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + p].length;
                }
                size_t estimatedKmerCnt = (totalLength + totalLength / 10) / 3;

                // Process current split if buffer has enough space.
                posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                    struct MmapedData<char> fastaFile = mmapData<char>(fastaList[fnaSplits[i].file_idx].path.c_str());
                    vector<string> cds;
                    vector<string> nonCds;
                    // ** Use provided CDS information **
                    for (size_t s_cnt = 0; s_cnt < fnaSplits[i].cnt; ++s_cnt) {
                        buffer = {const_cast<char *>(&fastaFile.data[fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + s_cnt].start]),
                                  static_cast<size_t>(fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + s_cnt].length)};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);
                        
                        // Mask low complexity regions
                        char *maskedSeq = nullptr;
                        if (par.maskMode) {
                            maskedSeq = new char[seq->seq.l + 1];
                            SeqIterator::maskLowComplexityRegions(seq->seq.s, maskedSeq, probMatrix, par.maskProb, subMat);
                        } else {
                            maskedSeq = seq->seq.s;
                        }

                        cout << "Processing " << seq->name.s << "\t" << seq->seq.l << "\t" << posToWrite << endl;
                        if (cdsInfoMap.find(string(seq->name.s)) != cdsInfoMap.end()) { // CDS provided
                            cds.clear();
                            nonCds.clear();
                            seqIterator.devideToCdsAndNonCds(
                                maskedSeq, seq->seq.l,
                                cdsInfoMap[string(seq->name.s)], cds, nonCds);
                            
                            // cout << "CDS: " << endl;
                            // for (size_t j = 0; j < cds.size(); j++) {
                            //     cout << cdsInfoMap[string(seq->name.s)][j].proteinId << endl;
                            //     cout << cds[j] << endl << endl;
                            // }

                            // cout << "Non-CDS: " << endl;
                            // for (size_t j = 0; j < nonCds.size(); j++) {
                            //     cout << nonCds[j] << endl << endl;
                            // }

                            // Translate CDS and extract metamers
                            for (size_t j = 0; j < cds.size(); j++) {
                                // if (cdsCnt == 1) {
                                //     cout << "First CDS: " << cds[j] << endl;
                                // }
                                seqIterator.translate(cds[j]);
                                // if (cdsCnt == 1) {
                                //     for (size_t pp = 0; pp < seqIterator.aaFrames[0].size() ; pp++) {
                                //         cout << seqIterator.aaFrames[0][pp] << " ";
                                //     }
                                //     cout << endl;
                                // }
                                uint32_t cdsId = __sync_fetch_and_add(&cdsCnt, 1);
                                queryCdsList[cdsId] = QueryCDS(cdsInfoMap[string(seq->name.s)][j].proteinId, cdsId, cds[j].size());
                                // queryCdsList.emplace_back(cdsInfoMap[string(seq->name.s)][j].proteinId, cdsId, cds[j].size());
                                seqIterator.computeMetamerF(cds[j].c_str(), 0, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, cdsId);
                            }

                            // Translate non-CDS and extract metamers
                            for (size_t j = 0; j < nonCds.size(); j++) {
                                if (observedNonCDSKmers.empty()) {
                                    seqIterator.translate(nonCds[j]);
                                    seqIterator.computeMetamerF(nonCds[j].c_str(), 0, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, 0);
                                    seqIterator.getMinHashList(observedNonCDSKmerHashes, nonCds[j].c_str());
                                } else {
                                    int frame = selectReadingFrame(observedNonCDSKmerHashes, nonCds[j].c_str(), seqIterator);
                                    if (frame < 3) { // Forward
                                        seqIterator.translate(nonCds[j], frame);
                                        seqIterator.computeMetamerF(nonCds[j].c_str(), frame, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, 0);
                                    } else { // Reverse complement
                                        reverseComplementStr.clear();
                                        reverseComplementStr = seqIterator.reverseComplement(nonCds[j]);
                                        seqIterator.translate(reverseComplementStr, frame - 3);
                                        seqIterator.computeMetamerF(reverseComplementStr.c_str(), frame - 3, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, 0);
                                    }
                                }
                            }
                        } else { // CDS not provided
                            // Translate the whole sequence
                            if (observedNonCDSKmers.empty()) {
                                seqIterator.translate(seq->seq.s);
                                seqIterator.getMinHashList(observedNonCDSKmerHashes, seq->seq.s);
                                seqIterator.computeMetamerF(seq->seq.s, 0, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, 0);                                    
                            } else {
                                int frame = selectReadingFrame(observedNonCDSKmerHashes, seq->seq.s, seqIterator);
                                if (frame < 3) { // Forward
                                    seqIterator.translate(seq->seq.s, frame);
                                    seqIterator.computeMetamerF(seq->seq.s, frame, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, 0);
                                } else { // Reverse complement
                                    reverseComplement = seqIterator.reverseComplement(seq->seq.s, seq->seq.l);
                                    seqIterator.translate(reverseComplement,  frame - 3);
                                    seqIterator.computeMetamerF(reverseComplement, frame - 3, kmerBuffer, posToWrite,int(processedSeqCnt[fnaSplits[i].file_idx] + fnaSplits[i].offset + s_cnt), fnaSplits[i].speciesID, 0);
                                    free(reverseComplement);
                                }
                            }                     
                        }
                        if (par.maskMode) { delete[] maskedSeq;}
                        kseq_destroy(seq);
                    }
                    __sync_fetch_and_add(&processedSplitCnt, 1);
#ifdef OPENMP
                    cout << omp_get_thread_num() << " Processed " << i << "th splits (" << processedSplitCnt << ")" << endl;
#endif
                    munmap(fastaFile.data, fastaFile.fileSize + 1);
                } else {
                    // Withdraw the reservation if the buffer is full.
                    tempChecker[i] = false;
                    __sync_fetch_and_add(&hasOverflow, 1);
                    __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
                }
            }
        }
    }
    return 0;
}

int FuncIndexer::getProteinId(Buffer<TargetMetamerF> &kmerBuffer) {
    Buffer<ProtMatch> matchBuffer(par.bufferSize);
    if (kmerMatcher->matchAAKmers(& kmerBuffer, & matchBuffer, dbDir)) {
        SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, 
        [] (const ProtMatch & a, const ProtMatch & b) {
            if (a.cdsId != b.cdsId) {
                return a.cdsId < b.cdsId;
            } else if (a.protId != b.protId) {
                return a.protId < b.protId;
            } else {
                return a.cdsPos < b.cdsPos;
            }
        });
        // print matchBuffer
        // uint32_t lastPost = 0;
        // for (size_t i = 0; i < matchBuffer.startIndexOfReserve; i++) {
        //     cout << matchBuffer.buffer[i].cdsId << " " 
        //          << matchBuffer.buffer[i].protId << " " 
        //          << matchBuffer.buffer[i].cdsPos << " "
        //          << queryCdsList[matchBuffer.buffer[i].cdsId].cdsId << " "
        //          << protIdMap[matchBuffer.buffer[i].protId] << " ";
        //     if (matchBuffer.buffer[i].cdsPos == lastPost) {
        //         cout << "Duplicated";
        //     } 
        //     cout << endl; 
        //     lastPost = matchBuffer.buffer[i].cdsPos;
        // }

        unordered_map<uint32_t, uint32_t> cdsId2protId;
        labelCdsWithProtId(matchBuffer, cdsId2protId);
        // // print cdsId2protId
        // cout << "print cdsId2protId" << endl;
        // for (auto & it : cdsId2protId) {
        //     cout << it.first << " " << it.second << " " <<  queryCdsList[it.first].cdsId << " " << protIdMap[it.second] << endl;
        // }
        // cout << "print queryCdsList" << endl;
        // for (size_t i = 0; i < queryCdsList.size(); i++) {
        //     cout << queryCdsList[i].cdsId << " " << queryCdsList[i].proteinId << endl;
        // }
        labelKmerWithProtId(kmerBuffer, cdsId2protId);
        cout << "K-mers from CDS are labeled with protein IDs" << endl;
        return 1;
    } else {
        return 0;
    }
}

void FuncIndexer::labelCdsWithProtId(Buffer<ProtMatch> & protMatch, unordered_map<uint32_t, uint32_t> & cdsId2protId) {
    // Divide ProtMatches into blocks for multi threading
    size_t seqNum = queryCdsList.size();
    // ProtMatchBlock *matchBlocks = new ProtMatchBlock[seqNum];
    vector<ProtMatchBlock> matchBlocks;
    size_t matchIdx = 0;
    // size_t blockIdx = 0;
    uint32_t currentCDS;
    uint32_t currentProt;
    size_t numOfMatch = protMatch.startIndexOfReserve;
    while (matchIdx < numOfMatch) {
        currentCDS = protMatch.buffer[matchIdx].cdsId;        
        size_t protCnt = 0;
        while((matchIdx < numOfMatch) && (currentCDS == protMatch.buffer[matchIdx].cdsId)) {
            currentProt = protMatch.buffer[matchIdx].protId;
            // matchBlocks[blockIdx].cdsId = currentCDS;
            // matchBlocks[blockIdx].protId = protCnt;
            // matchBlocks[blockIdx].start = matchIdx;
            size_t start = matchIdx;
            while ((matchIdx < numOfMatch) &&
                   (currentCDS == protMatch.buffer[matchIdx].cdsId) &&
                   (currentProt == protMatch.buffer[matchIdx].protId)) {
                ++matchIdx;
            }
            // matchBlocks[blockIdx++].end = matchIdx - 1;
            matchBlocks.emplace_back(start, matchIdx - 1, currentCDS, protCnt++);
            cds2protScoreMap[currentCDS].emplace_back(0, 0);
        }
    }
    // How about start a thread when a block is ready? Then, we can avoid the overhead of creating threads.
    for (size_t i = 0; i < matchBlocks.size(); ++i) {
            cout << i << " " << matchBlocks[i].cdsId << " " << matchBlocks[i].protId << " " <<  queryCdsList[matchBlocks[i].cdsId].cdsId << " " << protIdMap[protMatch.buffer[matchBlocks[i].start].protId] << " " << matchBlocks[i].start << " " << matchBlocks[i].end << endl;
    }

    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, protMatch, seqNum, cdsId2protId)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < matchBlocks.size(); ++i) {
            cds2protScoreMap[matchBlocks[i].cdsId][matchBlocks[i].protId] = scoreProteinMatches(matchBlocks[i].start, matchBlocks[i].end + 1, protMatch);
        }
    }
    
    for (auto & it : cds2protScoreMap) {
        float maxScore = -1;
        uint32_t bestProt = 0;
        for (size_t i = 0; i < it.second.size(); i++) {
            if (it.second[i].score > maxScore) {
                maxScore = it.second[i].score;
                bestProt = it.second[i].protIdx;
            }
        }
        cdsId2protId[it.first] = bestProt;
        // For debugging
        queryCdsList[it.first].protIdx = bestProt;
        queryCdsList[it.first].proteinId = protIdMap[bestProt];
        queryCdsList[it.first].score = maxScore;
    }
}

void FuncIndexer::labelKmerWithProtId(Buffer<TargetMetamerF> &kmerBuffer, const unordered_map<uint32_t, uint32_t> & cdsId2protId) {
#pragma omp parallel default(none), shared(kmerBuffer, cdsId2protId, cout)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < kmerBuffer.startIndexOfReserve; i++) {
            if (kmerBuffer.buffer[i].metamerF.protId != 0) {
                kmerBuffer.buffer[i].metamerF.protId = cdsId2protId.at(kmerBuffer.buffer[i].metamerF.protId);
            }
        }
    }
}

ProtScore FuncIndexer::scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch) {
    uint32_t protId = protMatch.buffer[start].protId;
    float score = 0;
    float len = (float) queryCdsList[protMatch.buffer[start].cdsId].cdsLength;
    vector<const ProtMatch *> consecutiveMatches;
    consecutiveMatches.reserve(protMatch.startIndexOfReserve);
    for (size_t i = start; i + 1 < end; i++) {
        size_t consecutive = 1;
        // MUST: protMatch.buffer[i].cdsPos != protMatch.buffer[i + 1].cdsPos
        while ((i + 1 < end) && protMatch.buffer[i].cdsPos + 3 == protMatch.buffer[i + 1].cdsPos) {
            ++consecutive;
            ++i;
        }
        if (consecutive >= 4) {
            for (size_t j = i - consecutive + 1; j <= i; j++) {
                consecutiveMatches.push_back(&protMatch.buffer[j]);
            }
        }
    }

    // Compute score
    int aminoAcidNum = (int) len / 3;
    auto *checker = new bool[aminoAcidNum + 1];
    memset(checker, 0, (aminoAcidNum + 1));
    for (size_t i = 0; i < consecutiveMatches.size(); i++) {
        uint32_t currPos = consecutiveMatches[i]->cdsPos / 3;
        checker[currPos] = true;
        checker[currPos + 1] = true;
        checker[currPos + 2] = true;
        checker[currPos + 3] = true;
        checker[currPos + 4] = true;
        checker[currPos + 5] = true;
        checker[currPos + 6] = true;
        checker[currPos + 7] = true;
    }

    for (int i = 0; i < aminoAcidNum; i++) {
        score += 3.0f * checker[i];
    }
    delete[] checker;    
    return ProtScore(protId, score / len);
}

bool FuncIndexer::sortTargetMetamerF(const TargetMetamerF & a, const TargetMetamerF & b){
    return a.metamerF.metamer < b.metamerF.metamer;
}

void FuncIndexer::loadProtIdMap() {
    ifstream protIdMapFile(protIdMapFileName);
    if (!protIdMapFile.is_open()) {
        cerr << "Cannot open " << protIdMapFileName << endl;
        exit(1);
    }
    string line;
    while (getline(protIdMapFile, line)) {
        istringstream iss(line);
        uint32_t protIdx;
        string protId;
        iss >> protId >> protIdx;
        protIdMap[protIdx] = protId;
    }
    protIdMapFile.close();
}