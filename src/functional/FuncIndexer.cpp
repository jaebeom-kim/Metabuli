#include "FuncIndexer.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "KmerMatcher.h"
#include "LocalUtil.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include "ProteinDbIndexer.h"
#include "SeqIterator.h"
#include <cstdint>
#include <cstring>
#include <iterator>

FuncIndexer::FuncIndexer(const LocalParameters &par) : IndexCreator(par), par(par) {
    protDBFileName = par.proteinDB + "/protein.mtbl";
    protDbSplitFileName = par.proteinDB + "/prot_split.mtbl";
    protIdMapFileName = par.proteinDB + "/prtIdMap.mtbl";
    protIdx2taxIdFileName = dbDir + "/id2taxonomy.mtbl";
    protIdx2unirefIdFileName = dbDir + "/id2uniref.mtbl";
    protIdx2taxIdFileName_debug = dbDir + "/id2taxonomy_debug.mtbl";
    protIdx2unirefIdFileName_debug = dbDir + "/id2uniref_debug.mtbl";

    // loadProtIdMap();
    kmerMatcher = new KmerMatcher(par, taxonomy);

    queryCdsList.resize(par.bufferSize / 50);
    
}

FuncIndexer::~FuncIndexer() {
    
}

void FuncIndexer::createIndex() {
    // Read through FASTA files and make blocks of sequences to be processed by each thread
    if (par.accessionLevel) {
        makeBlocksForParallelProcessing_accession_level();
    } else {
        makeBlocksForParallelProcessing();
    }
    writeTaxIdList();
    
    loadCdsInfo(par.cdsInfo);
    uint32_t nonCdsIdx = lastProtIdx + 1;
    
    // Process the splits until all are processed
    bool * tempChecker = new bool[fnaSplits.size()];
    bool * completionChecker = new bool[fnaSplits.size()];
    memset(tempChecker, 0, fnaSplits.size() * sizeof(bool));
    memset(completionChecker, 0, fnaSplits.size() * sizeof(bool));
    // fill_n(completionChecker, fnaSplits.size(), false);
    size_t totalProcessedSplitCnt = 0;
    Buffer<ExtractedMetamer> metamerBuffer(par.bufferSize);
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(totalProcessedSplitCnt < fnaSplits.size()) {
        queryCdsList.clear();
        queryCdsList.resize(par.bufferSize / 50);

        size_t processedSplitCnt = 0;
        memset(metamerBuffer.buffer, 0, metamerBuffer.bufferSize * sizeof(ExtractedMetamer));
        fillTargetKmerBuffer(metamerBuffer, tempChecker, processedSplitCnt, nonCdsIdx);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(metamerBuffer.buffer,
                      metamerBuffer.buffer + metamerBuffer.startIndexOfReserve, 
                      [] (const ExtractedMetamer & a, const ExtractedMetamer & b) {
                                return a.metamer.metamer < b.metamer.metamer;});
                              
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

        SORT_PARALLEL(metamerBuffer.buffer,
                      metamerBuffer.buffer + metamerBuffer.startIndexOfReserve,
                      sortExtractedMetamer);
        
        // Reduce redundancy
        time_t reduction = time(nullptr);
        auto * uniqKmerIdx = new size_t[metamerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
        cout<<"Time spent for reducing redundancy: "<<(double) (reduction - sort) << endl;
        
        // Write the k-mers to the database
        if(processedSplitCnt == fnaSplits.size() && numOfFlush == 0){
            writeTargetFilesAndSplits(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
        }
        delete[] uniqKmerIdx;
    }
    delete[] tempChecker;
    delete[] completionChecker;
    writeTaxonomyDB();
    writeDbParameters();
    LocalUtil::writeMappingFile(this->protIdx2taxId, this->protIdx2taxIdFileName);
    LocalUtil::writeMappingFile(this->protIdx2unirefId, this->protIdx2unirefIdFileName);
    LocalUtil::writeMappingFile_text(this->protIdx2taxId, this->protIdx2taxIdFileName_debug);
    LocalUtil::writeMappingFile_text(this->protIdx2unirefId, this->protIdx2unirefIdFileName_debug);
}


size_t FuncIndexer::fillTargetKmerBuffer(Buffer<ExtractedMetamer> &kmerBuffer,
                                         bool * tempChecker,
                                         size_t &processedSplitCnt,
                                         uint32_t & nonCdsIdx) {
    uint32_t cdsCnt = 1;
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(kmerBuffer, tempChecker, processedSplitCnt, hasOverflow, cout, cdsCnt, nonCdsIdx)
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

                        // cout << "Processing " << seq->name.s << "\t" << seq->seq.l << "\t" << posToWrite << endl;
                        TaxID currentTaxId = foundAcc2taxid[string(seq->name.s)];
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
                                queryCdsList[cdsId] = QueryCDS(cdsInfoMap[string(seq->name.s)][j].protId, cdsId, cds[j].size());
                                // queryCdsList.emplace_back(cdsInfoMap[string(seq->name.s)][j].proteinId, cdsId, cds[j].size());
                                seqIterator.computeMetamerF(
                                    cds[j].c_str(), // DNA sequence string
                                    0,            // frame
                                    kmerBuffer,       // extracted k-mer buffer
                                    posToWrite,       // position to write
                                    fnaSplits[i].speciesID,
                                    cdsInfoMap[string(seq->name.s)][j].protId);
                            }

                            // Translate non-CDS and extract metamers
                            for (size_t j = 0; j < nonCds.size(); j++) {
                                uint32_t nonCdsId = __sync_fetch_and_add(&nonCdsIdx, 1);
                                protIdx2taxId[nonCdsId] = currentTaxId;
                                if (observedNonCDSKmers.empty()) {
                                    seqIterator.translate(nonCds[j]);
                                    seqIterator.computeMetamerF(nonCds[j].c_str(), 0, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId);
                                    seqIterator.getMinHashList(observedNonCDSKmerHashes, nonCds[j].c_str());
                                } else {
                                    int frame = selectReadingFrame(observedNonCDSKmerHashes, nonCds[j].c_str(), seqIterator);
                                    if (frame < 3) { // Forward
                                        seqIterator.translate(nonCds[j], frame);
                                        seqIterator.computeMetamerF(nonCds[j].c_str(), frame, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId);
                                    } else { // Reverse complement
                                        reverseComplementStr.clear();
                                        reverseComplementStr = seqIterator.reverseComplement(nonCds[j]);
                                        seqIterator.translate(reverseComplementStr, frame - 3);
                                        seqIterator.computeMetamerF(reverseComplementStr.c_str(), frame - 3, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId);
                                    }
                                }
                            }
                        } else { // CDS not provided
                            // Translate the whole sequence
                            uint32_t nonCdsId = __sync_fetch_and_add(&nonCdsIdx, 1);
                            protIdx2taxId[nonCdsId] = currentTaxId;
                            if (observedNonCDSKmers.empty()) {
                                seqIterator.translate(seq->seq.s);
                                seqIterator.getMinHashList(observedNonCDSKmerHashes, seq->seq.s);
                                seqIterator.computeMetamerF(seq->seq.s, 0, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId);                                    
                            } else {
                                int frame = selectReadingFrame(observedNonCDSKmerHashes, seq->seq.s, seqIterator);
                                if (frame < 3) { // Forward
                                    seqIterator.translate(seq->seq.s, frame);
                                    seqIterator.computeMetamerF(seq->seq.s, frame, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId);
                                } else { // Reverse complement
                                    reverseComplement = seqIterator.reverseComplement(seq->seq.s, seq->seq.l);
                                    seqIterator.translate(reverseComplement,  frame - 3);
                                    seqIterator.computeMetamerF(reverseComplement, frame - 3, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId);
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

int FuncIndexer::getProteinId(Buffer<ExtractedMetamer> &kmerBuffer) {
    Buffer<ProtMatch> matchBuffer(par.bufferSize);
    if (kmerMatcher->matchAAKmers(& kmerBuffer, & matchBuffer, lastProtIdx, dbDir)) {
        cout << "Number of matches: " << matchBuffer.startIndexOfReserve << endl;
        SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, 
        [] (const ProtMatch & a, const ProtMatch & b) {
            if (a.queryProtId != b.queryProtId) {
                return a.queryProtId < b.queryProtId;
            } else if (a.targetProtId != b.targetProtId) {
                return a.targetProtId < b.targetProtId;
            } else {
                return a.queryProtPos < b.queryProtPos;
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

        // unordered_map<uint32_t, uint32_t> cdsId2protId;
        mapQueryProt2TargetProt(matchBuffer);
        labelExtractedMetamerWithUniRefId(kmerBuffer);
        // // print cdsId2protId
        // cout << "print cdsId2protId" << endl;
        // for (auto & it : cdsId2protId) {
        //     cout << it.first << " " << it.second << " " <<  queryCdsList[it.first].cdsId << " " << protIdMap[it.second] << endl;
        // }
        // cout << "print cdsId2protId done" << endl;
        // labelKmerWithProtId(kmerBuffer, cdsId2protId);
        cout << "K-mers from CDS are labeled with protein IDs" << endl;
        return 1;
    } else {
        return 0;
    }
}

// int FuncIndexer::getProteinId(Buffer<TargetMetamerF> &kmerBuffer) {
//     Buffer<ProtMatch> matchBuffer(par.bufferSize);
//     if (kmerMatcher->matchAAKmers(& kmerBuffer, & matchBuffer, dbDir)) {
//         cout << "Number of matches: " << matchBuffer.startIndexOfReserve << endl;
//         SORT_PARALLEL(matchBuffer.buffer, matchBuffer.buffer + matchBuffer.startIndexOfReserve, 
//         [] (const ProtMatch & a, const ProtMatch & b) {
//             if (a.queryProtId != b.queryProtId) {
//                 return a.queryProtId < b.queryProtId;
//             } else if (a.targetProtId != b.targetProtId) {
//                 return a.targetProtId < b.targetProtId;
//             } else {
//                 return a.queryProtPos < b.queryProtPos;
//             }
//         });
//         // print matchBuffer
//         // uint32_t lastPost = 0;
//         // for (size_t i = 0; i < matchBuffer.startIndexOfReserve; i++) {
//         //     cout << matchBuffer.buffer[i].cdsId << " " 
//         //          << matchBuffer.buffer[i].protId << " " 
//         //          << matchBuffer.buffer[i].cdsPos << " "
//         //          << queryCdsList[matchBuffer.buffer[i].cdsId].cdsId << " "
//         //          << protIdMap[matchBuffer.buffer[i].protId] << " ";
//         //     if (matchBuffer.buffer[i].cdsPos == lastPost) {
//         //         cout << "Duplicated";
//         //     } 
//         //     cout << endl; 
//         //     lastPost = matchBuffer.buffer[i].cdsPos;
//         // }

//         unordered_map<uint32_t, uint32_t> cdsId2protId;
//         mapQueryProt2TargetProt(matchBuffer, cdsId2protId);
//         // // print cdsId2protId
//         // cout << "print cdsId2protId" << endl;
//         // for (auto & it : cdsId2protId) {
//         //     cout << it.first << " " << it.second << " " <<  queryCdsList[it.first].cdsId << " " << protIdMap[it.second] << endl;
//         // }
//         // cout << "print cdsId2protId done" << endl;
//         labelKmerWithProtId(kmerBuffer, cdsId2protId);
//         cout << "K-mers from CDS are labeled with protein IDs" << endl;
//         return 1;
//     } else {
//         return 0;
//     }
// }

void FuncIndexer::mapQueryProt2TargetProt(Buffer<ProtMatch> & protMatch) {
    // Divide ProtMatches into blocks for multi threading
    size_t seqNum = queryCdsList.size();
    vector<ProtMatchBlock> matchBlocks;
    size_t matchIdx = 0;
    uint32_t currentQueryProtIdx;
    uint32_t currentProt;
    size_t numOfMatch = protMatch.startIndexOfReserve;
    while (matchIdx < numOfMatch) {
        currentQueryProtIdx = protMatch.buffer[matchIdx].queryProtId;        
        size_t protCnt = 0;
        while((matchIdx < numOfMatch) && (currentQueryProtIdx == protMatch.buffer[matchIdx].queryProtId)) {
            currentProt = protMatch.buffer[matchIdx].targetProtId;
            size_t start = matchIdx;
            while ((matchIdx < numOfMatch) &&
                   (currentQueryProtIdx == protMatch.buffer[matchIdx].queryProtId) &&
                   (currentProt == protMatch.buffer[matchIdx].targetProtId)) {
                ++matchIdx;
            }
            matchBlocks.emplace_back(start, matchIdx - 1, currentQueryProtIdx, protCnt++);
            query2targetProtScMap[currentQueryProtIdx].emplace_back(0, 0);
        }
    }
    // How about start a thread when a block is ready? Then, we can avoid the overhead of creating threads.
    // for (size_t i = 0; i < matchBlocks.size(); ++i) {
    //         cout << i << " " << matchBlocks[i].cdsId << " " << matchBlocks[i].protId << " " <<  queryCdsList[matchBlocks[i].cdsId].cdsId << " " << protIdMap[protMatch.buffer[matchBlocks[i].start].protId] << " " << matchBlocks[i].start << " " << matchBlocks[i].end << endl;
    // }

    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, protMatch, seqNum)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < matchBlocks.size(); ++i) {
            query2targetProtScMap[matchBlocks[i].queryProtId][matchBlocks[i].targetProtCnt] = scoreProteinMatches(matchBlocks[i].start, matchBlocks[i].end + 1, protMatch);
        }
    }
    
    for (auto & it : query2targetProtScMap) {
        float maxScore = -1;
        uint32_t bestProt = 0;
        for (size_t i = 0; i < it.second.size(); i++) {
            if (it.second[i].score > maxScore) {
                maxScore = it.second[i].score;
                bestProt = it.second[i].targetProtIdx;
            }
        }
        protIdx2unirefId[it.first] = bestProt;
        // cdsId2protId[it.first] = bestProt;
        // For debugging
        queryCdsList[it.first].assigend_uniref_idx = bestProt;
        // queryCdsList[it.first].assigend_uniref_id = uniRefIdMap[bestProt];
        queryCdsList[it.first].score = maxScore;
    }
}


void FuncIndexer::labelExtractedMetamerWithUniRefId(Buffer<ExtractedMetamer> &kmerBuffer) {
#pragma omp parallel default(none), shared(kmerBuffer, cout)
    {
#pragma omp for schedule(dynamic, 10000)
        for (size_t i = 0; i < kmerBuffer.startIndexOfReserve; i++) {
            if (kmerBuffer.buffer[i].metamer.id > lastProtIdx) {
                kmerBuffer.buffer[i].unirefId = protIdx2unirefId[kmerBuffer.buffer[i].metamer.id];
            }
        }
    }
    
}

void FuncIndexer::labelKmerWithProtId(Buffer<TargetMetamerF> &kmerBuffer, const unordered_map<uint32_t, uint32_t> & cdsId2protId) {
#pragma omp parallel default(none), shared(kmerBuffer, cdsId2protId, cout)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < kmerBuffer.startIndexOfReserve; i++) {
            // cdsKmerCnt[kmerBuffer.buffer[i].metamerF.protId]++;
            if (kmerBuffer.buffer[i].metamerF.protId != 0) {
                kmerBuffer.buffer[i].metamerF.protId = cdsId2protId.at(kmerBuffer.buffer[i].metamerF.protId);
            }
        }
    }
    
}

ProtScore FuncIndexer::scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch) {
    uint32_t protId = protMatch.buffer[start].targetProtId;
    float score = 0;
    float len = (float) queryCdsList[protMatch.buffer[start].queryProtId].cdsLength;
    vector<const ProtMatch *> consecutiveMatches;
    consecutiveMatches.reserve(protMatch.startIndexOfReserve);
    for (size_t i = start; i + 1 < end; i++) {
        size_t consecutive = 1;
        // MUST: protMatch.buffer[i].cdsPos != protMatch.buffer[i + 1].cdsPos
        while ((i + 1 < end) && protMatch.buffer[i].queryProtPos + 3 == protMatch.buffer[i + 1].queryProtPos) {
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
        uint32_t currPos = consecutiveMatches[i]->queryProtPos / 3;
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

inline bool FuncIndexer::sortTargetMetamerF(const TargetMetamerF & a, const TargetMetamerF & b){
    if (a.metamerF.metamer != b.metamerF.metamer) {
        return a.metamerF.metamer < b.metamerF.metamer;
    }
    if (a.speciesId != b.speciesId) {
        return a.speciesId < b.speciesId;
    }
    if (a.metamerF.protId != b.metamerF.protId) {
        return a.metamerF.protId < b.metamerF.protId;
    }
    return a.metamerF.seqId < b.metamerF.seqId;
}

inline bool FuncIndexer::sortExtractedMetamer(const ExtractedMetamer & a, const ExtractedMetamer & b){
    if (a.metamer.metamer != b.metamer.metamer) {
        return a.metamer.metamer < b.metamer.metamer;
    }
    if (a.speciesId != b.speciesId) {
        return a.speciesId < b.speciesId;
    }
    if (a.unirefId != b.unirefId) {
        return a.unirefId < b.unirefId;
    }
    return a.metamer.id < b.metamer.id;
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
        uniRefIdMap[protIdx] = protId;
    }
    protIdMapFile.close();
}

void FuncIndexer::reduceRedundancy(Buffer<ExtractedMetamer> &kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqKmerCnt) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t i = kmerBuffer.startIndexOfReserve - 1; i != 0; i--){
        if(kmerBuffer.buffer[i].metamer.metamer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = i + 1;
            break;
        }
    }

    SeqIterator seqIterator(par);
    
    // for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
    //     if(kmerBuffer.buffer[i].speciesId == 0){
    //         cout << i << " " << kmerBuffer.buffer[i].metamerF.metamer << " " << kmerBuffer.buffer[i].speciesId << " " << kmerBuffer.buffer[i].metamerF.protId << " " << kmerBuffer.buffer[i].metamerF.seqId << endl;
    //         seqIterator.printKmerInDNAsequence(kmerBuffer.buffer[i].metamerF.metamer); cout << endl;
    //     }
    // }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].speciesId != 0){
            startIdx = i;
            cout << "startIdx: " << startIdx << endl;
            break;
        }
    }

    // for (size_t i = startIdx; i < kmerBuffer.startIndexOfReserve; i++) {
    //     cout << i << " "; seqIterator.printKmerInDNAsequence(kmerBuffer.buffer[i].metamerF.metamer);
    //     cout << " " << kmerBuffer.buffer[i].speciesId << " " << kmerBuffer.buffer[i].metamerF.protId << " " << kmerBuffer.buffer[i].metamerF.seqId << endl;
    // }
    cout << "K-mer counts before reducing redundancy: " << kmerBuffer.startIndexOfReserve - startIdx << endl;

    ExtractedMetamer * lookingKmer = nullptr;
    size_t lookingIndex = 0;
    bool isEnd = false;
    vector<TaxID> taxIds;

    lookingKmer = & kmerBuffer.buffer[startIdx];
    lookingIndex = startIdx;
    for (size_t i = startIdx + 1; i < kmerBuffer.startIndexOfReserve; i++) {
        taxIds.clear();
        taxIds.push_back(protIdx2taxId[lookingKmer->metamer.id]);
        while((lookingKmer->metamer.metamer == kmerBuffer.buffer[i].metamer.metamer)
              && (lookingKmer->speciesId == kmerBuffer.buffer[i].speciesId)
              && (lookingKmer->unirefId == kmerBuffer.buffer[i].unirefId)){
            taxIds.push_back(protIdx2taxId[kmerBuffer.buffer[i].metamer.id]);
            i++;
            if (i == kmerBuffer.startIndexOfReserve) {
                isEnd = true;
                break;
            }
        }
        if(taxIds.size() > 1) {
            protIdx2taxId[lookingKmer->metamer.id] = taxonomy->LCA(taxIds)->taxId;
        } 
        uniqeKmerIdx[uniqKmerCnt] = lookingIndex;
        uniqKmerCnt ++;
        if (isEnd) break;
        lookingKmer = & kmerBuffer.buffer[i];
        lookingIndex = i;
    }

    if (!((kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 2].metamer.metamer == kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].metamer.metamer)
        && (kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 2].speciesId == kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].speciesId)
        && (kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 2].unirefId == kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].unirefId))) {
        uniqeKmerIdx[uniqKmerCnt] = kmerBuffer.startIndexOfReserve - 1;
        uniqKmerCnt ++;
    }
    cout << "K-mer counts after reducing redundancy: " << uniqKmerCnt << endl;
}

void FuncIndexer::writeTargetFilesAndSplits(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName = dbDir + "/diffIdx";
    string infoFileName = dbDir + "/info";
    string splitFileName = dbDir + "/split";

    // Make splits
    FILE * diffIdxSplitFile = fopen(splitFileName.c_str(), "wb");
    DiffIdxSplit splitList[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    size_t splitWidth = uniqKmerCnt / par.splitNum;
    size_t remainder = uniqKmerCnt % par.splitNum;
    size_t splitCnt = 1;
    size_t start = 0;
    for (size_t i = 1; i < (size_t) par.splitNum; i++) {
        start = start + splitWidth;
        if (remainder > 0) {
            start++;
            remainder--;
        }
        for (size_t j = start; j + 1 < start + splitWidth; j++) {
            if (AminoAcidPart(kmerBuffer.buffer[uniqKmerIdx[j]].metamer.metamer) 
                != AminoAcidPart(kmerBuffer.buffer[uniqKmerIdx[j + 1]].metamer.metamer)) {
                splitList[splitCnt].ADkmer = kmerBuffer.buffer[uniqKmerIdx[j + 1]].metamer.metamer;
                cout << splitList[splitCnt].ADkmer << endl;
                splitCnt++;
                break;
            }
        }
    }

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;


    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * size_t(par.bufferSize));
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer.buffer[uniqKmerIdx[i]].metamer.id, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[i]].metamer.metamer, diffIdxFile,
                   diffIdxBuffer, this->bufferSize, localBufIdx, totalDiffIdx);
        lastKmer = kmerBuffer.buffer[uniqKmerIdx[i]].metamer.metamer;
        if((splitIdx < splitCnt) && (lastKmer == splitList[splitIdx].ADkmer)){
            splitList[splitIdx].diffIdxOffset = totalDiffIdx;
            splitList[splitIdx].infoIdxOffset = write;
            splitIdx ++;
        }
    }
    
    cout<<"K-mer counts recorded on disk: "<< write << endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    printIndexSplitList(splitList);
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, diffIdxSplitFile);

    free(diffIdxBuffer);
    fclose(diffIdxSplitFile);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerBuffer.startIndexOfReserve = 0;
}

// This function sort the TargetKmerBuffer, do redundancy reducing task, write the differential index of them
void FuncIndexer::writeTargetFiles(Buffer<ExtractedMetamer> &kmerBuffer,
                                   const size_t * uniqeKmerIdx,
                                   size_t & uniqKmerCnt){
    string diffIdxFileName;
    string infoFileName;
    diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    FILE * infoFile = fopen(infoFileName.c_str(), "wb");
    if (diffIdxFile == nullptr || infoFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000'000);
    size_t localBufIdx = 0;
    uint64_t lastKmer = 0;
    size_t write = 0;

    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fwrite(& kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.id, sizeof (TargetKmerInfo), 1, infoFile);
        write++;
        getDiffIdx(lastKmer, kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.id, diffIdxFile, diffIdxBuffer, this->bufferSize, localBufIdx);
        lastKmer = kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.metamer;
    }
    cout<<"K-mer counts recorded on disk: "<< write << endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerBuffer.startIndexOfReserve = 0;
}

void FuncIndexer::reduceRedundancy(Buffer<TargetMetamerF> &kmerBuffer, size_t * uniqeKmerIdx, size_t & uniqKmerCnt) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t i = kmerBuffer.startIndexOfReserve - 1; i != 0; i--){
        if(kmerBuffer.buffer[i].metamerF.metamer != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = i + 1;
            break;
        }
    }

    SeqIterator seqIterator(par);
    cout << "Before reducing redundancy: " << kmerBuffer.startIndexOfReserve << endl;
    // for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
    //     if(kmerBuffer.buffer[i].speciesId == 0){
    //         cout << i << " " << kmerBuffer.buffer[i].metamerF.metamer << " " << kmerBuffer.buffer[i].speciesId << " " << kmerBuffer.buffer[i].metamerF.protId << " " << kmerBuffer.buffer[i].metamerF.seqId << endl;
    //         seqIterator.printKmerInDNAsequence(kmerBuffer.buffer[i].metamerF.metamer); cout << endl;
    //     }
    // }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++){
        if(kmerBuffer.buffer[i].speciesId != 0){
            startIdx = i;
            cout << "startIdx: " << startIdx << endl;
            break;
        }
    }

    // for (size_t i = startIdx; i < kmerBuffer.startIndexOfReserve; i++) {
    //     cout << i << " "; seqIterator.printKmerInDNAsequence(kmerBuffer.buffer[i].metamerF.metamer);
    //     cout << " " << kmerBuffer.buffer[i].speciesId << " " << kmerBuffer.buffer[i].metamerF.protId << " " << kmerBuffer.buffer[i].metamerF.seqId << endl;
    // }

    cout << kmerBuffer.startIndexOfReserve - startIdx << endl;

    TargetMetamerF * lookingKmer = nullptr;
    size_t lookingIndex = 0;
    bool isEnd = false;
    vector<TaxID> taxIds;

    lookingKmer = & kmerBuffer.buffer[startIdx];
    lookingIndex = startIdx;
    // cout << "Unique k-mers:" << endl;
    for (size_t i = startIdx + 1; i < kmerBuffer.startIndexOfReserve; i++) {
        taxIds.clear();
        taxIds.push_back(taxIdList[lookingKmer->metamerF.seqId]);
        while((lookingKmer->metamerF.metamer == kmerBuffer.buffer[i].metamerF.metamer)
              && (lookingKmer->speciesId == kmerBuffer.buffer[i].speciesId)
              && (lookingKmer->metamerF.protId == kmerBuffer.buffer[i].metamerF.protId)){
            taxIds.push_back(taxIdList[kmerBuffer.buffer[i].metamerF.seqId]);
            i++;
            if (i == kmerBuffer.startIndexOfReserve) {
                isEnd = true;
                break;
            }
        }
        if(taxIds.size() > 1){
            lookingKmer->metamerF.seqId = taxonomy->LCA(taxIds)->taxId;
        } else {
            lookingKmer->metamerF.seqId = taxIds[0];
        }
        uniqeKmerIdx[uniqKmerCnt] = lookingIndex;
        // cout << uniqKmerCnt << " " << lookingIndex << " "; seqIterator.printKmerInDNAsequence(lookingKmer->metamerF.metamer);
        // cout << " " << lookingKmer->speciesId << " " << lookingKmer->metamerF.protId << " " << lookingKmer->metamerF.seqId << endl;
        uniqKmerCnt ++;
        if (isEnd) break;
        lookingKmer = & kmerBuffer.buffer[i];
        lookingIndex = i;
    }

    if (!((kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 2].metamerF.metamer == kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].metamerF.metamer)
        && (kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 2].speciesId == kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].speciesId)
        && (kmerBuffer.buffer[kmerBuffer.startIndexOfReserve- 2].metamerF.protId == kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].metamerF.protId))) {
        kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].metamerF.seqId = taxIdList[kmerBuffer.buffer[kmerBuffer.startIndexOfReserve - 1].metamerF.seqId];
        uniqeKmerIdx[uniqKmerCnt] = kmerBuffer.startIndexOfReserve - 1;
        uniqKmerCnt ++;
    }
    cout << "After reducing redundancy: " << uniqKmerCnt << endl;
}

void FuncIndexer::loadCdsInfo(const string & cdsInfoFileList) {
    uint32_t prtId = 1;
    ifstream cdsInfoList(cdsInfoFileList);
    if (cdsInfoList.is_open()) {
        string cdsInfoFile;
        while (getline(cdsInfoList, cdsInfoFile)) { // Read each CDS info file
            ifstream cdsInfo(cdsInfoFile);
            if (cdsInfo.is_open()) {
                string line;
                while (getline(cdsInfo, line)) { // Read each line of the CDS info file
                    if (line[0] == '>') { // Check if the line starts with ">"
                        // Get the accession number between the '|' and '.'.
                        size_t start = line.find('|') + 1;
                        size_t end = line.find('.', start);
                        string accession = line.substr(start, end - start + 2);
                        TaxID taxonomyId = foundAcc2taxid[accession];
                        // cout << "Accession: " << accession << endl;
                        int frame = 1;
                        while (true) {
                            start = line.find('[', end) + 1;
                            end = line.find(']', start);
                            if (start == string::npos) { break;}
                            size_t equalPos = line.find('=', start);
                            string feature = line.substr(start, equalPos - start);
                            string value = line.substr(equalPos + 1, end - equalPos - 1);
                            if (feature == "pseudo") {
                                break;
                            } else if (feature == "protein" && value == "hypothetical protein") {
                                break;
                            } else if (feature == "frame") {
                                frame = stoi(value);
                            } else if (feature == "protein_id") {
                                // cout << "Protein ID: " << value << "\t" << frame << endl;
                                protIdMap[prtId] = value; // protein idx to string protein id
                                protIdx2taxId[prtId] = taxonomyId;
                                cdsInfoMap[accession].emplace_back(CDSinfo(prtId++, frame));
                            } else if (feature == "location") {
                                // cout << "Location: " << value << endl;
                                // Check if the location is complement
                                size_t complementPos = value.find('c');
                                bool isComplement = (complementPos != string::npos);
                                if (isComplement) {
                                    cdsInfoMap[accession].back().isComplement = true;
                                    value = value.substr(complementPos + 11, value.size() - complementPos - 12);
                                } else {
                                    cdsInfoMap[accession].back().isComplement = false;
                                }

                                // Check if spliced
                                size_t joinPos = value.find('j');
                                if (joinPos != string::npos) {
                                    value = value.substr(joinPos + 5, value.size() - joinPos - 6);
                                }
                                
                                // Load the locations
                                size_t commaPos = value.find(',');
                                size_t dotPos;
                                string locationBegin, locationEnd;
                                while (commaPos != string::npos) {
                                    dotPos = value.find('.');
                                    locationBegin = value.substr(0, dotPos);
                                    locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                                    
                                    // Check < and > signs
                                    if (locationBegin[0] == '<') {
                                        locationBegin = locationBegin.substr(1, locationBegin.size() - 1);
                                    }
                                    if (locationEnd[0] == '>') {
                                        locationEnd = locationEnd.substr(1, locationEnd.size() - 1);
                                    }

                                    cdsInfoMap[accession].back().loc.emplace_back(stoi(locationBegin), stoi(locationEnd));
                                    value = value.substr(commaPos + 1, value.size() - commaPos - 1);
                                    commaPos = value.find(',');
                                }
                                dotPos = value.find('.');
                                locationBegin = value.substr(0, dotPos);
                                locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                                // Check < and > signs
                                if (locationBegin[0] == '<') {
                                    locationBegin = locationBegin.substr(1, locationBegin.size() - 1);
                                }
                                if (locationEnd[0] == '>') {
                                    locationEnd = locationEnd.substr(1, locationEnd.size() - 1);
                                }
                                cdsInfoMap[accession].back().loc.emplace_back(stoi(locationBegin), stoi(locationEnd));

                                // Frame correction
                                if (frame != 1) {
                                    if (!isComplement) {
                                        cdsInfoMap[accession].back().loc[0].first += frame - 1;
                                    } else {
                                        cdsInfoMap[accession].back().loc.back().second -= frame - 1;
                                    }
                                }
                                break;
                            } 
                        }
                    }
                }
            } else {
                Debug(Debug::ERROR) << "Cannot open file " << cdsInfoFile << "\n";
                EXIT(EXIT_FAILURE);
            }
        }
    } else {
        Debug(Debug::ERROR) << "Cannot open file " << cdsInfoFileList << "\n";
        EXIT(EXIT_FAILURE);
    }
    lastProtIdx = prtId - 1;
}