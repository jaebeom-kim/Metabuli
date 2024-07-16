#include "FuncIndexer.h"
#include "KaijuWrapper.h"
#include "Kmer.h"
#include "KmerMatcher.h"
#include "LocalUtil.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include "ProdigalWrapper.h"
#include "ProteinDbIndexer.h"
#include "SeqIterator.h"
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iterator>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>

FuncIndexer::FuncIndexer(const LocalParameters &par) : IndexCreator(par) {
    // UniRef Idx
    unirefId2IdxFileName = par.proteinDB + "/unirefId2Idx.mtbl";
    unirefIdx2taxIdFileName = par.proteinDB + "/unirefIdx2taxId.mtbl";

    // Output files
    regionId2taxIdFileName = dbDir + "/id2taxonomy.mtbl";
    regionId2unirefIdFileName = dbDir + "/id2uniref.mtbl";
    regionId2taxIdFileName_debug = dbDir + "/id2taxonomy_debug.mtbl";
    regionId2unirefIdFileName_debug = dbDir + "/id2uniref_debug.mtbl";
    
    ncbi2gtdbFileName = par.ncbi2gtdb;

    // Load UniRef mapping files
    // loadUniRefId2Idx(); // UniRef Id -> UniRef Idx
    // LocalUtil::loadMappingFile<uint32_t, int>(unirefIdx2taxIdFileName, unirefIdx2taxId); // UniRef Idx -> Taxonomy ID
    
    kaiju = new KaijuWrapper(par);
}

FuncIndexer::~FuncIndexer() {
    delete kaiju;
}

void FuncIndexer::createIndex() {
    // Read through FASTA files and make blocks of sequences to be processed by each thread
    if (par.accessionLevel) {
        makeBlocksForParallelProcessing_accession_level();
    } else {
        makeBlocksForParallelProcessing();
    }
    writeTaxIdList();
    generateTaxId2SpeciesIdMap();
    if (!par.cdsInfo.empty()) {
        loadCdsInfo(par.cdsInfo);   
    }
    uint32_t nonCdsIdx = lastProtIdx + 1;
    
    // Process the splits until all are processed
    bool * tempChecker = new bool[fnaSplits.size()];
    bool * completionChecker = new bool[fnaSplits.size()];
    memset(tempChecker, 0, fnaSplits.size() * sizeof(bool));
    memset(completionChecker, 0, fnaSplits.size() * sizeof(bool));
    // fill_n(completionChecker, fnaSplits.size(), false);
    size_t totalProcessedSplitCnt = 0;
    Buffer<ExtractedMetamer> metamerBuffer(par.bufferSize);
    uint32_t regionId = 0;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(totalProcessedSplitCnt < fnaSplits.size()) {
        // queryCdsList.clear();
        // queryCdsList.resize(par.bufferSize / 50);

        size_t processedSplitCnt = 0;
        memset(metamerBuffer.buffer, 0, metamerBuffer.bufferSize * sizeof(ExtractedMetamer));
        if (!par.cdsInfo.empty()) {
            fillTargetKmerBuffer(metamerBuffer, tempChecker, processedSplitCnt, nonCdsIdx);
        } else {
            fillTargetKmerBufferUsingProdigal(metamerBuffer, tempChecker, processedSplitCnt, regionId);
        }
        totalProcessedSplitCnt += processedSplitCnt;
        
        SORT_PARALLEL(metamerBuffer.buffer,
                      metamerBuffer.buffer + metamerBuffer.startIndexOfReserve,
                      sortExtractedMetamer);
        
        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[metamerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundancy(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
        
        // Write the k-mers to the database
        if(processedSplitCnt == fnaSplits.size() && numOfFlush == 0){
            writeTargetFilesAndSplits(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
            writeTargetFilesAndSplits2(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
        } else {
            writeTargetFiles(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
            writeTargetFiles2(metamerBuffer, uniqKmerIdx, uniqKmerCnt);
        }
        delete[] uniqKmerIdx;
    }

    mergeDeltaIndexFiles();
    delete[] tempChecker;
    delete[] completionChecker;
    writeTaxonomyDB();
    writeDbParameters();
    LocalUtil::writeMappingFile(this->regionId2taxId, this->regionId2taxIdFileName);
    LocalUtil::writeMappingFile(this->regionId2unirefId, this->regionId2unirefIdFileName);
    // LocalUtil::writeMappingFile_text(this->protIdx2taxId, this->protIdx2taxIdFileName_debug);
    // LocalUtil::writeMappingFile_text(this->protIdx2unirefId, this->protIdx2unirefIdFileName_debug);
}


size_t FuncIndexer::fillTargetKmerBuffer(Buffer<ExtractedMetamer> &kmerBuffer,
                                         bool * tempChecker,
                                         size_t &processedSplitCnt,
                                         uint32_t & nonCdsIdx) {
    // uint32_t cdsCnt = 1;
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(kmerBuffer, tempChecker, processedSplitCnt, hasOverflow, cout, nonCdsIdx)
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
                                // uint32_t cdsId = __sync_fetch_and_add(&cdsCnt, 1);
                                queryCodingRegionMap[cdsInfoMap[string(seq->name.s)][j].protId] = QueryCodingRegionInfo(cdsInfoMap[string(seq->name.s)][j].protId, cds[j].size());
                                // queryCdsList[cdsId] = QueryCDS(cdsInfoMap[string(seq->name.s)][j].protId, cds[j].size());
                                // queryCdsList.emplace_back(cdsInfoMap[string(seq->name.s)][j].proteinId, cdsId, cds[j].size());
                                seqIterator.computeMetamerF(
                                    cds[j].c_str(), // DNA sequence string
                                    0,            // frame
                                    kmerBuffer,       // extracted k-mer buffer
                                    posToWrite,       // position to write
                                    fnaSplits[i].speciesID,
                                    cdsInfoMap[string(seq->name.s)][j].protId,
                                    1);
                            }

                            // Translate non-CDS and extract metamers
                            for (size_t j = 0; j < nonCds.size(); j++) {
                                uint32_t nonCdsId = __sync_fetch_and_add(&nonCdsIdx, 1);
                                regionId2taxId[nonCdsId] = currentTaxId;
                                if (observedNonCDSKmers.empty()) {
                                    seqIterator.translate(nonCds[j]);
                                    seqIterator.computeMetamerF(nonCds[j].c_str(), 0, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId, 0);
                                    seqIterator.getMinHashList(observedNonCDSKmerHashes, nonCds[j].c_str());
                                } else {
                                    int frame = selectReadingFrame(observedNonCDSKmerHashes, nonCds[j].c_str(), seqIterator);
                                    if (frame < 3) { // Forward
                                        seqIterator.translate(nonCds[j], frame);
                                        seqIterator.computeMetamerF(nonCds[j].c_str(), frame, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId, 0);
                                    } else { // Reverse complement
                                        reverseComplementStr.clear();
                                        reverseComplementStr = seqIterator.reverseComplement(nonCds[j]);
                                        seqIterator.translate(reverseComplementStr, frame - 3);
                                        seqIterator.computeMetamerF(reverseComplementStr.c_str(), frame - 3, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId, 0);
                                    }
                                }
                            }
                        } else { // CDS not provided
                            // Translate the whole sequence
                            uint32_t nonCdsId = __sync_fetch_and_add(&nonCdsIdx, 1);
                            regionId2taxId[nonCdsId] = currentTaxId;
                            if (observedNonCDSKmers.empty()) {
                                seqIterator.translate(seq->seq.s);
                                seqIterator.getMinHashList(observedNonCDSKmerHashes, seq->seq.s);
                                seqIterator.computeMetamerF(seq->seq.s, 0, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId, 0);                                    
                            } else {
                                int frame = selectReadingFrame(observedNonCDSKmerHashes, seq->seq.s, seqIterator);
                                if (frame < 3) { // Forward
                                    seqIterator.translate(seq->seq.s, frame);
                                    seqIterator.computeMetamerF(seq->seq.s, frame, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId, 0);
                                } else { // Reverse complement
                                    reverseComplement = seqIterator.reverseComplement(seq->seq.s, seq->seq.l);
                                    seqIterator.translate(reverseComplement,  frame - 3);
                                    seqIterator.computeMetamerF(reverseComplement, frame - 3, kmerBuffer, posToWrite, fnaSplits[i].speciesID, nonCdsId, 0);
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

size_t FuncIndexer::fillTargetKmerBufferUsingProdigal(Buffer<ExtractedMetamer> &kmerBuffer,
                                                      bool * tempChecker,
                                                      size_t &processedSplitCnt,
                                                      uint32_t & regionId) {
    cout << "Using Prodigal to extract CDS" << endl;
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(kmerBuffer, tempChecker, processedSplitCnt, hasOverflow, cout, regionId)
    {
        ProbabilityMatrix probMatrix(*subMat);
        SeqIterator seqIterator(par);
        size_t posToWrite;
        priority_queue<uint64_t> observedNonCDSKmerHashes;
        // char *reverseComplement;
        string reverseComplementStr;
        kseq_buffer_t buffer;
        kseq_t *seq;
        vector<uint64_t> observedNonCDSKmers;
        std::unordered_map<uint32_t, int> localRegionId2taxId;
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
                ProdigalWrapper * prodigal = new ProdigalWrapper();
                // Process current split if buffer has enough space.
                posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                    struct MmapedData<char> fastaFile = mmapData<char>(fastaList[fnaSplits[i].file_idx].path.c_str());
                    // Load sequence for prodigal training
                    buffer = {const_cast<char *>(&fastaFile.data[fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].training].start]),
                              static_cast<size_t>(fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].training].length)};
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    size_t lengthOfTrainingSeq = seq->seq.l;

                    // Train prodigal
                    prodigal->is_meta = 0;
                    if (lengthOfTrainingSeq < 100'000) {
                        prodigal->is_meta = 1;
                        prodigal->trainMeta(seq->seq.s);
                    } else {
                        prodigal->trainASpecies(seq->seq.s);
                    }
                    kseq_destroy(seq);

                    vector<string> cds;
                    vector<string> nonCds;
                    vector<string> masked_cds;
                    vector<string> masked_nonCds;
                    // ** Use predicted CDS information **
                    for (size_t s_cnt = 0; s_cnt < fnaSplits[i].cnt; ++s_cnt) {
                        cds.clear();
                        nonCds.clear();
                        masked_cds.clear();
                        masked_nonCds.clear();

                        buffer = {const_cast<char *>(&fastaFile.data[fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + s_cnt].start]),
                                  static_cast<size_t>(fastaList[fnaSplits[i].file_idx].sequences[fnaSplits[i].offset + s_cnt].length)};
                        seq = kseq_init(&buffer);
                        kseq_read(seq);

                        // Get predicted CDS
                        prodigal->getPredictedGenes(seq->seq.s);

                        // Mask low complexity regions
                        char *maskedSeq = nullptr;
                        if (par.maskMode) {
                            maskedSeq = new char[seq->seq.l + 1];
                            SeqIterator::maskLowComplexityRegions(seq->seq.s, maskedSeq, probMatrix, par.maskProb, subMat);
                        } else {
                            maskedSeq = seq->seq.s;
                        }

                        TaxID currentTaxId = foundAcc2taxid[string(seq->name.s)];

                        seqIterator.devideToCdsAndNonCds(
                            seq->seq.s, seq->seq.l,
                            prodigal, cds, nonCds);

                        seqIterator.devideToCdsAndNonCds(
                            maskedSeq, seq->seq.l,
                            prodigal, masked_cds, masked_nonCds);
                        
                        uint32_t localRegionId = __sync_fetch_and_add(&regionId, cds.size() + nonCds.size());

                        for (size_t j = 0; j < cds.size(); j++) {
                            // Get UniRef ID of the CDS
                            // cout << ">" << localRegionId << endl;
                            // cout << cds[j] << endl;
                            // cout << ">" << localRegionId << "_aa" <<endl;
                            // seqIterator.printTranslation(cds[j]);
                            uint32_t uniRefId = getUniRefIdx(cds[j], currentTaxId);
                            seqIterator.translate(masked_cds[j]);
                            localRegionId2taxId[localRegionId] = currentTaxId;
                            queryCodingRegionMap[localRegionId] = QueryCodingRegionInfo(localRegionId, cds[j].size());
                            seqIterator.computeMetamerF(cds[j].c_str(),         // DNA sequence string
                                                      0,                      // frame
                                                          kmerBuffer,             // extracted k-mer buffer
                                                          posToWrite,             // position to write
                                                fnaSplits[i].speciesID, 
                                                     localRegionId,
                                                     uniRefId);
                            localRegionId++;
                        }
                        
                        for (size_t j = 0; j < nonCds.size(); j++) {
                            // uint32_t nonCdsId = __sync_fetch_and_add(&nonCdsIdx, 1);
                            localRegionId2taxId[localRegionId] = currentTaxId;
                            if (observedNonCDSKmers.empty()) {
                                seqIterator.translate(masked_nonCds[j]);
                                seqIterator.computeMetamerF(masked_nonCds[j].c_str(), 0, kmerBuffer, posToWrite, fnaSplits[i].speciesID, localRegionId, 0);
                                seqIterator.getMinHashList(observedNonCDSKmerHashes, masked_nonCds[j].c_str());
                            } else {
                                int frame = selectReadingFrame(observedNonCDSKmerHashes, masked_nonCds[j].c_str(), seqIterator);
                                if (frame < 3) { // Forward
                                    seqIterator.translate(masked_nonCds[j], frame);
                                    seqIterator.computeMetamerF(masked_nonCds[j].c_str(), frame, kmerBuffer, posToWrite, fnaSplits[i].speciesID, localRegionId, 0);
                                } else { // Reverse complement
                                    reverseComplementStr.clear();
                                    reverseComplementStr = seqIterator.reverseComplement(masked_nonCds[j]);
                                    seqIterator.translate(reverseComplementStr, frame - 3);
                                    seqIterator.computeMetamerF(reverseComplementStr.c_str(), frame - 3, kmerBuffer, posToWrite, fnaSplits[i].speciesID, localRegionId, 0);
                                }
                            }
                            localRegionId++;
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
                delete prodigal;
            }
        }
        // Write the local regionId2taxId to the global regionId2taxId
#pragma omp critical
        {
            for (auto & it : localRegionId2taxId) {
                regionId2taxId[it.first] = it.second;
            }
        }
    }
    return 0;
}

uint32_t FuncIndexer::getUniRefIdx(const string & cds, TaxID taxId) {
    std::set<char *> match_ids;
    int len;
    kaiju->classify_length(cds, match_ids, len);
    // // print match_ids
    // cout <<">";
    // for (auto & it : match_ids) {
    //     cout << it << ",";
    // }
    // cout << len << endl;
    // string seq;
    // kaiju->translate(cds, seq);
    // cout << seq << endl;
    // return 0;

    if (match_ids.empty()) {
        return 0;
    } else if (match_ids.size() == 1) {
        return uniRefId2Idx[*match_ids.begin()];
    } else {
        // Find the closest (the minimum number of edges in a taxonomy tree) UniRef ID
        map<TaxID, unordered_set<uint32_t>> taxId2UniRefIdx;
        int minDist = INT32_MAX;
        for (auto & it : match_ids) {
            TaxID currUniRefTaxId = unirefIdx2taxId[uniRefId2Idx[it]];
            int dist = taxonomy->getDistance(currUniRefTaxId, taxId);
            if (dist < minDist) {
                minDist = dist;
                taxId2UniRefIdx.clear();
                taxId2UniRefIdx[taxId].insert(uniRefId2Idx[it]);
            } else if (dist == minDist) {
                taxId2UniRefIdx[taxId].insert(uniRefId2Idx[it]);
            } 
        }

        if (taxId2UniRefIdx.size() == 1) {
            return *(taxId2UniRefIdx.begin()->second.begin());
        } else {
            int level = INT32_MIN;
            TaxID mostSpecificTaxId;
            for (auto & it : taxId2UniRefIdx) {
                int currLevel = taxonomy->getLevel(it.first);
                if (currLevel > level) {
                    level = currLevel;
                    mostSpecificTaxId = it.first;
                }
            }
            return *(taxId2UniRefIdx[mostSpecificTaxId].begin());    
        }
    }
}




// ProtScore FuncIndexer::scoreProteinMatches(size_t start, size_t end, const Buffer<ProtMatch> & protMatch) {
//     uint32_t protId = protMatch.buffer[start].targetProtId;
//     float score = 0;
//     float len = (float) queryCodingRegionMap[protMatch.buffer[start].queryProtId].cdsLength;
//     vector<const ProtMatch *> consecutiveMatches;
//     consecutiveMatches.reserve(protMatch.startIndexOfReserve);
//     for (size_t i = start; i + 1 < end; i++) {
//         size_t consecutive = 1;
//         // MUST: protMatch.buffer[i].cdsPos != protMatch.buffer[i + 1].cdsPos
//         while ((i + 1 < end) && protMatch.buffer[i].queryProtPos + 3 == protMatch.buffer[i + 1].queryProtPos) {
//             ++consecutive;
//             ++i;
//         }
//         if (consecutive >= 4) {
//             for (size_t j = i - consecutive + 1; j <= i; j++) {
//                 consecutiveMatches.push_back(&protMatch.buffer[j]);
//             }
//         }
//     }

//     // Compute score
//     int aminoAcidNum = (int) len / 3;
//     auto *checker = new bool[aminoAcidNum + 1];
//     memset(checker, 0, (aminoAcidNum + 1));
//     for (size_t i = 0; i < consecutiveMatches.size(); i++) {
//         uint32_t currPos = consecutiveMatches[i]->queryProtPos / 3;
//         checker[currPos] = true;
//         checker[currPos + 1] = true;
//         checker[currPos + 2] = true;
//         checker[currPos + 3] = true;
//         checker[currPos + 4] = true;
//         checker[currPos + 5] = true;
//         checker[currPos + 6] = true;
//         checker[currPos + 7] = true;
//     }

//     for (int i = 0; i < aminoAcidNum; i++) {
//         score += 3.0f * checker[i];
//     }
//     delete[] checker;    
//     return ProtScore(protId, score / len);
// }

// inline bool FuncIndexer::sortTargetMetamerF(const TargetMetamerF & a, const TargetMetamerF & b){
//     if (a.metamerF.metamer != b.metamerF.metamer) {
//         return a.metamerF.metamer < b.metamerF.metamer;
//     }
//     if (a.speciesId != b.speciesId) {
//         return a.speciesId < b.speciesId;
//     }
//     if (a.metamerF.protId != b.metamerF.protId) {
//         return a.metamerF.protId < b.metamerF.protId;
//     }
//     return a.metamerF.seqId < b.metamerF.seqId;
// }

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




void FuncIndexer::loadUniRefId2Idx() {
    // Each line: UniRef90_K0JAJ8 187136235
    // string 2 uint32_t
    ifstream uniRefId2IdxFile(unirefId2IdxFileName);
    if (!uniRefId2IdxFile.is_open()) {
        cerr << "Cannot open " << unirefId2IdxFileName << endl;
        exit(1);
    }
    string line;
    while (getline(uniRefId2IdxFile, line)) {
        istringstream iss(line);
        string protId;
        uint32_t protIdx;
        iss >> protId >> protIdx;
        uniRefId2Idx[protId] = protIdx;
    }
    uniRefId2IdxFile.close();
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
        taxIds.push_back(regionId2taxId[lookingKmer->metamer.id]);
        while((lookingKmer->metamer.metamer == kmerBuffer.buffer[i].metamer.metamer)
              && (lookingKmer->speciesId == kmerBuffer.buffer[i].speciesId)
              && (lookingKmer->unirefId == kmerBuffer.buffer[i].unirefId)){
            taxIds.push_back(regionId2taxId[kmerBuffer.buffer[i].metamer.id]);
            i++;
            if (i == kmerBuffer.startIndexOfReserve) {
                isEnd = true;
                break;
            }
        }
        if(taxIds.size() > 1) {
            regionId2taxId[lookingKmer->metamer.id] = taxonomy->LCA(taxIds)->taxId;
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
                   diffIdxBuffer, par.bufferSize, localBufIdx, totalDiffIdx);
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

void FuncIndexer::writeTargetFilesAndSplits2(Buffer<ExtractedMetamer> &kmerBuffer, const size_t * uniqKmerIdx, size_t & uniqKmerCnt){
    string diffIdxFileName = dbDir + "/deltaIdx.mtbl";
    string splitFileName = dbDir + "/deltaIdxSplits.mtbl";

    // Make splits
    FILE * diffIdxSplitFile = fopen(splitFileName.c_str(), "wb");
    DeltaIdxOffset * offsetList = new DeltaIdxOffset[par.splitNum];
    memset(offsetList, 0, sizeof(DiffIdxSplit) * par.splitNum);
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
                offsetList[splitCnt].metamer = kmerBuffer.buffer[uniqKmerIdx[j + 1]].metamer;
                // cout << offsetList[splitCnt].ADkmer << endl;
                splitCnt++;
                break;
            }
        }
    }

    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    if (diffIdxFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;


    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * size_t(par.bufferSize));
    size_t localBufIdx = 0;
    Metamer previousMetamer;
    size_t splitIdx = 1;
    size_t totalDiffIdx = 0;
    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        fillDeltaIndexing(previousMetamer, kmerBuffer.buffer[uniqKmerIdx[i]].metamer, diffIdxFile,
                   diffIdxBuffer, par.bufferSize, localBufIdx, totalDiffIdx);
        previousMetamer = kmerBuffer.buffer[uniqKmerIdx[i]].metamer;
        if((splitIdx < splitCnt) && (previousMetamer.metamer == offsetList[splitIdx].metamer.metamer)){
            offsetList[splitIdx].offset = totalDiffIdx;
            // splitList[splitIdx].infoIdxOffset = write;
            splitIdx ++;
        }
    }
    
    cout<<"K-mer counts recorded on disk: "<< uniqKmerCnt << endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    // printIndexSplitList(splitList);
    fwrite(offsetList, sizeof(DeltaIdxOffset), par.splitNum, diffIdxSplitFile);

    free(diffIdxBuffer);
    free(offsetList);
    fclose(diffIdxSplitFile);
    fclose(diffIdxFile);
    kmerBuffer.startIndexOfReserve = 0;
}

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
        getDiffIdx(lastKmer, kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.id, diffIdxFile, diffIdxBuffer, par.bufferSize, localBufIdx);
        lastKmer = kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.metamer;
    }
    cout<<"K-mer counts recorded on disk: "<< write << endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    fclose(infoFile);
    kmerBuffer.startIndexOfReserve = 0;
}

void FuncIndexer::writeTargetFiles2(Buffer<ExtractedMetamer> &kmerBuffer,
                                    const size_t * uniqeKmerIdx,
                                    size_t & uniqKmerCnt){
    string diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_deltaIdx.mtbl";
    FILE * diffIdxFile = fopen(diffIdxFileName.c_str(), "wb");
    if (diffIdxFile == nullptr){
        cout<<"Cannot open the file for writing target DB"<<endl;
        return;
    }
    numOfFlush++;

    uint16_t *diffIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * 10'000'000'000);
    size_t localBufIdx = 0;
    size_t write = 0;
    Metamer previousMetamer = Metamer();

    for(size_t i = 0; i < uniqKmerCnt ; i++) {
        write++;
        fillDeltaIndexing(previousMetamer, kmerBuffer.buffer[uniqeKmerIdx[i]].metamer, diffIdxFile, diffIdxBuffer, par.bufferSize, localBufIdx);
        previousMetamer = kmerBuffer.buffer[uniqeKmerIdx[i]].metamer;
    }
    cout<<"K-mer counts recorded on disk: "<< write << endl;

    flushKmerBuf(diffIdxBuffer, diffIdxFile, localBufIdx);
    free(diffIdxBuffer);
    fclose(diffIdxFile);
    kmerBuffer.startIndexOfReserve = 0;
}

void FuncIndexer::fillDeltaIndexing(const Metamer & previousMetamer,
                                    const Metamer & currentMetamer,
                                    FILE* handleKmerTable,
                                    uint16_t * deltaIndexBuffer,
                                    size_t bufferSize,
                                    size_t & localBufIdx) {
    bitset<96> diff = Metamer::substract(currentMetamer, previousMetamer);                               
    // uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[7];
    // uint16_t buffer[5];
    int idx = 5;
    buffer[6] = SET_END_FLAG(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
    // buffer[6] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    diff >>= 15U;
    while (diff.any()) {
        uint16_t toWrite = GET_15_BITS(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
        diff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    writeDiffIdx(deltaIndexBuffer, handleKmerTable, (buffer + idx + 1), (6 - idx), localBufIdx, bufferSize);
}

void FuncIndexer::fillDeltaIndexing(const Metamer & previousMetamer,
                                    const Metamer & currentMetamer,
                                    FILE* handleKmerTable,
                                    uint16_t * deltaIndexBuffer,
                                    size_t bufferSize,
                                    size_t & localBufIdx,
                                    size_t & totalBufferIdx) {
    bitset<96> diff = Metamer::substract(currentMetamer, previousMetamer);                               
    uint16_t buffer[7];
    int idx = 5;
    buffer[6] = SET_END_FLAG(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
    diff >>= 15U;
    while (diff.any()) {
        uint16_t toWrite = GET_15_BITS(static_cast<uint16_t>((diff & bitset<96>(0x7FFF)).to_ulong()));
        diff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    totalBufferIdx += 6 - idx;
    writeDiffIdx(deltaIndexBuffer, handleKmerTable, (buffer + idx + 1), (6 - idx), localBufIdx, bufferSize);
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
                                // regionId2taxId[prtId] = taxonomyId;
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

// This function write the target DB in text format for debugging
void FuncIndexer::writeTargetFilesInText(Buffer<ExtractedMetamer> &kmerBuffer,
                                         const size_t * uniqeKmerIdx,
                                         size_t & uniqKmerCnt) {
    string textDbFileName = dbDir + "/textDb";
    ofstream textDb(textDbFileName);
    if (!textDb) {
        cout << "Cannot open the file for writing target DB" << endl;
        return;
    }
    for (size_t i = 0; i < uniqKmerCnt; i++) {
        textDb << bitset<64>(kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.metamer) << " "
               << kmerBuffer.buffer[uniqeKmerIdx[i]].metamer.id << endl;
    }
    textDb.close();
}

void FuncIndexer::testDeltaIndexing(const Buffer<ExtractedMetamer> & kmerBuffer,
                                    const size_t * uniqueKmerIdx,
                                    size_t uniqueKmerCnt) {
    string deltaFileName = dbDir + "/delta.debug";
    ofstream deltaFile(deltaFileName);
    if (!deltaFile) {
        cout << "Cannot open the file for writing delta" << endl;
        return;
    }
    Metamer previousMetamer(0, 0);
    for (size_t i = 0; i < uniqueKmerCnt; i++) {
        bitset<96> delta = Metamer::substract(kmerBuffer.buffer[uniqueKmerIdx[i]].metamer, previousMetamer);
        Metamer currentMetamer = previousMetamer.add(delta);
        deltaFile << bitset<64>(currentMetamer.metamer) << " "
                  << currentMetamer.id << endl;
    }
    deltaFile.close();
}

void FuncIndexer::mergeDeltaIndexFiles(){
    string mergedFileName = dbDir + "/deltaIdx.mtbl";
    string splitFileName = dbDir + "/deltaIdxSplits.mtbl";
    FILE * mergedFile = fopen(mergedFileName.c_str(), "wb");
    FILE * idxSplitFile = fopen(splitFileName.c_str(), "wb");

    // Buffer
    uint16_t * kmerBuffer = (uint16_t *)malloc(sizeof(uint16_t) * par.bufferSize);
    size_t kmerBufferIdx = 0;
    size_t totalBufferIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    Metamer * lookingMetamers = new Metamer[numOfFlush];
    auto * maxIdxOfEachFiles = new size_t[numOfFlush];
    struct MmapedData<uint16_t> * deltaIdxList = new struct MmapedData<uint16_t>[numOfFlush];
    // struct MmapedData<size_t> * kmerCntFileList = new struct MmapedData<size_t>[numOfFlush];
    size_t * currPositions = new size_t[numOfFlush];
    memset(currPositions, 0, numOfFlush * sizeof(size_t));

    for (size_t file = 0; file < numOfFlush; file++) {
        deltaIdxList[file] = mmapData<uint16_t>((dbDir + "/" + to_string(file) + "_deltaIdx.mtbl").c_str());
        // kmerCntFileList[file] = mmapData<size_t>((dbDir + "/" + to_string(file) + "_protein_kmer_num.mtbl").c_str());
        maxIdxOfEachFiles[file] = deltaIdxList[file].fileSize / sizeof(uint16_t);
        numOfKmerBeforeMerge += maxIdxOfEachFiles[file];
    }    
    cout << "Number of k-mers before merge: " << numOfKmerBeforeMerge << endl;
    
    // Temporary offsets
    size_t splitWidth = numOfKmerBeforeMerge / par.splitNum;
    size_t remainder = numOfKmerBeforeMerge % par.splitNum;
    size_t * tempOffsetList = new size_t[par.splitNum + 1];
    tempOffsetList[0] = 0;
    for (size_t i = 1; i < (size_t) par.splitNum; i++) {
        tempOffsetList[i] = tempOffsetList[i - 1] + splitWidth;
        if (remainder > 0) {
            tempOffsetList[i]++;
            remainder--;
        }
    }
    tempOffsetList[par.splitNum] = UINT64_MAX;

    // Index splits to fill in
    DeltaIdxOffset * offsetList = new DeltaIdxOffset[par.splitNum];
    // ProtIdxSplit splitList[par.splitNum];
    memset(offsetList, 0, sizeof(DeltaIdxOffset) * par.splitNum);

    // Get the first k-mer of each ProtSplitFile
    for(size_t file = 0; file < numOfFlush; file++){
        Metamer previousMetamer(0, 0);
        lookingMetamers[file] = KmerMatcher::getNextTargetKmer(previousMetamer, deltaIdxList[file].data, currPositions[file]);
    }
    size_t idxOfMin = getSmallestMetamer(lookingMetamers, numOfFlush);
    Metamer entry(0, 0);
    size_t remainedFiles = numOfFlush;
    size_t splitCnt = 0;
    Metamer lastWrittenMetamer(0, 0);
    size_t offsetListIdx = 0;
    bool check = false;
    bool end = false;
    size_t processedKmers = 0;
    vector<TaxID> taxIds;
    while(true) {
        // Update the entry
        entry = lookingMetamers[idxOfMin];

        // Update looking metamers
        if (currPositions[idxOfMin] == maxIdxOfEachFiles[idxOfMin]) {
            lookingMetamers[idxOfMin] = Metamer(UINT64_MAX, UINT32_MAX);
            remainedFiles--;
            if(remainedFiles == 0) break;
        } else {
            lookingMetamers[idxOfMin] = KmerMatcher::getNextTargetKmer(entry, deltaIdxList[idxOfMin].data, currPositions[idxOfMin]);
        }

        // Find the smallest k-mer
        idxOfMin = getSmallestMetamer(lookingMetamers, numOfFlush);
        // Metamer smallestMetamer = lookingMetamers[idxOfMin];

        // Scan redundancy
        taxIds.clear();
        taxIds.push_back(regionId2taxId[entry.id]);
        while(entry.metamer == lookingMetamers[idxOfMin].metamer
           && taxId2speciesId[regionId2taxId[entry.id]] == taxId2speciesId[regionId2taxId[lookingMetamers[idxOfMin].id]]
           && regionId2unirefId[entry.id] == regionId2unirefId[lookingMetamers[idxOfMin].id]) {
            taxIds.push_back(regionId2taxId[lookingMetamers[idxOfMin].id]);
            if (currPositions[idxOfMin] == maxIdxOfEachFiles[idxOfMin]) {
                lookingMetamers[idxOfMin] = Metamer(UINT64_MAX, UINT32_MAX);
                remainedFiles--;
                if(remainedFiles == 0) {
                    end = true;
                    break;
                }
            } else {
                lookingMetamers[idxOfMin] = KmerMatcher::getNextTargetKmer(lookingMetamers[idxOfMin], deltaIdxList[idxOfMin].data, currPositions[idxOfMin]);
            }
            idxOfMin = getSmallestMetamer(lookingMetamers, numOfFlush);
        }
        if(taxIds.size() > 1) {
            regionId2taxId[entry.id] = taxonomy->LCA(taxIds)->taxId;
        }

        fillDeltaIndexing(lastWrittenMetamer,
                          entry,
                          mergedFile, kmerBuffer,
                          par.bufferSize,
                          kmerBufferIdx,
                          totalBufferIdx);
        
        // Check if offset is reached
        if (processedKmers == tempOffsetList[offsetListIdx]) {
            if (processedKmers == 0) {
                offsetList[splitCnt++] = {entry, totalBufferIdx};
                offsetListIdx++;
            } else {
                check = true;
                offsetListIdx++;
            }
        }
        processedKmers++;

        // If new offset is reached and new AA is seen, it is time to split
        if (check && ((lastWrittenMetamer.metamer & DNA_MASK) != (entry.metamer & DNA_MASK))) {
            offsetList[splitCnt++] = {entry, totalBufferIdx};
            check = false;
        }

        // // Update looking k-mers
        // if (currPositions[idxOfMin] == maxIdxOfEachFiles[idxOfMin]) {
        //     lookingMetamers[idxOfMin] = Metamer(UINT64_MAX, UINT32_MAX);
        //     remainedFiles--;
        //     if(remainedFiles == 0) break;
        // } else {
        //     lookingMetamers[idxOfMin] = KmerMatcher::getNextTargetKmer(smallestMetamer, deltaIdxList[idxOfMin].data, currPositions[idxOfMin]);
        // }
        lastWrittenMetamer = entry;
    }

    // flush buffer
    IndexCreator::flushKmerBuf(kmerBuffer, mergedFile, kmerBufferIdx);

    // write split
    fwrite(offsetList, sizeof(DeltaIdxOffset), par.splitNum, idxSplitFile);
    free(kmerBuffer);
    fclose(mergedFile);
    fclose(idxSplitFile);
 
    for(size_t file = 0; file < numOfFlush; file++){
        munmap(deltaIdxList[file].data, deltaIdxList[file].fileSize + 1);
    }

    delete[] tempOffsetList;
    delete[] deltaIdxList;
    delete[] lookingMetamers;
    delete[] currPositions;
    delete[] maxIdxOfEachFiles;
}

size_t FuncIndexer::getSmallestMetamer(const Metamer * lookingMetamers, size_t numOfFiles) {
    size_t idxOfMin = 0;
    for (size_t i = 1; i < numOfFiles; i++) {
        if (lookingMetamers[i] < lookingMetamers[idxOfMin]) {
            idxOfMin = i;
        }
    }
    return idxOfMin;
}

void FuncIndexer::generateTaxId2SpeciesIdMap() {
    for (size_t i = 0; i < taxIdList.size(); i ++) {
        TaxonNode const * taxon = taxonomy->taxonNode(taxIdList[i]);
        TaxID speciesId = taxonomy->getTaxIdAtRank(taxIdList[i], "species");
        speciesTaxIds.insert(speciesId);
        while (taxon->taxId != speciesId) {
            taxId2speciesId[taxon->taxId] = speciesId;
            taxon = taxonomy->taxonNode(taxon->parentTaxId);
        }
        taxId2speciesId[taxon->taxId] = speciesId;
        taxId2speciesId[taxIdList[i]] = speciesId; 
    }

}

