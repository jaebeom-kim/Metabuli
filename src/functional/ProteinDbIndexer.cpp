#include "ProteinDbIndexer.h"
#include "FileMerger.h"
#include "IndexCreator.h"
#include "Kmer.h"
#include "LocalUtil.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "printBinary.h"
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <ostream>

ProteinDbIndexer::ProteinDbIndexer(const LocalParameters &par) : par(par){
    mask_getId = (1ULL << 28) - 1;
    mask_getAA = ~((1ULL << 28) - 1);
    dbDir = par.filenames[0];
    aaKmerBuffer = new Buffer<uint64_t>(par.bufferSize);
    proteinIndexSplitFileName = dbDir + "/prot_split.mtbl";
    prtId2taxIdFileName = dbDir + "/unirefIdx2taxId.mtbl";
    bufferSize = par.bufferSize;
    flushCnt = 0;

    // Read taxonomy
    taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    
}

ProteinDbIndexer::~ProteinDbIndexer() {
    delete aaKmerBuffer;
}

void ProteinDbIndexer::index() {
    splitFasta(par.proteinDB.c_str(), this->sequenceBlocks);
    // LocalUtil::writeMappingFile(proteinId2Index, dbDir + "/prtIdMap.mtbl");
    writePrtIdMap();
    LocalUtil::writeMappingFile<uint32_t, int>(proteinIndex2taxonomyId, prtId2taxIdFileName);
    LocalUtil::writeMappingFile_text(proteinIndex2taxonomyId, prtId2taxIdFileName + ".txt");
    
    size_t numOfSplits = sequenceBlocks.size();
    bool * splitChecker = new bool[numOfSplits];
    fill_n(splitChecker, numOfSplits, false);
    size_t processedSplitCnt = 0;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedSplitCnt < numOfSplits) {
        fillAAKmerBuffer(splitChecker, processedSplitCnt);
        SORT_PARALLEL(aaKmerBuffer->buffer, aaKmerBuffer->buffer + aaKmerBuffer->startIndexOfReserve);
        
        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[aaKmerBuffer->startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        reduceRedundantAAKmers(uniqKmerIdx, uniqKmerCnt);

        writeProteinDeltaIndexFile(uniqKmerIdx, uniqKmerCnt, processedSplitCnt == numOfSplits && flushCnt == 0);
        delete[] uniqKmerIdx;
    }
    if (flushCnt > 1) {
        // Merge
        // mergeProteinIndexFiles();
        mergeProteinDeltaIndexFiles();
    }
    delete [] splitChecker;
}

size_t ProteinDbIndexer::fillAAKmerBuffer(bool *checker, size_t & processedSplitCnt) {
    aaKmerBuffer->startIndexOfReserve = 0;
    int hasOverflow = 0;
#pragma omp parallel default(none), shared(aaKmerBuffer, checker, processedSplitCnt, hasOverflow, par, cout)
    {
        SeqIterator seqIterator(par);
        size_t posToWrite;
        kseq_buffer_t buffer;
        kseq_t *seq;
        struct MmapedData<char> fastaFile = mmapData<char>(par.proteinDB.c_str());
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < sequenceBlocks.size(); i++) {
            if (!checker[i] && !hasOverflow) {
                checker[i] = true;
                size_t estimatedKmerCnt = sequenceBlocks[i].seqLength - 7;
                if(estimatedKmerCnt < 4) {
                    __sync_fetch_and_add(&processedSplitCnt, 1);  
                    continue;
                }
                posToWrite = this->aaKmerBuffer->reserveMemory(estimatedKmerCnt);
                if (posToWrite + estimatedKmerCnt < this->aaKmerBuffer->bufferSize) {
                    buffer = {const_cast<char *>(&fastaFile.data[sequenceBlocks[i].start]),
                              static_cast<size_t>(sequenceBlocks[i].length)};
                    seq = kseq_init(&buffer);
                    kseq_read(seq);
                    // cout << "Processing " << seq->name.s << " " << this->proteinId2Index[seq->name.s] << endl;
                    seqIterator.computeAAKmer(this->aaKmerBuffer, seq->seq.s, seq->seq.l, posToWrite, this->proteinId2Index[seq->name.s]);
                    kseq_destroy(seq);
                    __sync_fetch_and_add(&processedSplitCnt, 1);  
                } else {
                    checker[i] = false;
                    __sync_fetch_and_add(&hasOverflow, 1);
                    __sync_fetch_and_sub(&aaKmerBuffer->startIndexOfReserve, estimatedKmerCnt);
                }
            }
        }
        munmap(fastaFile.data, fastaFile.fileSize + 1);
    }
    return 0;
}

void ProteinDbIndexer::reduceRedundantAAKmers(size_t * uniqeKmerIdx, size_t & uniqueKmerCnt) {
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = aaKmerBuffer->startIndexOfReserve - 1; checkN != 0; checkN--){
        if(aaKmerBuffer->buffer[checkN] != UINT64_MAX){
            aaKmerBuffer->startIndexOfReserve = checkN + 1;
            break;
        }
    }
    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for(size_t i = 0; i < aaKmerBuffer->startIndexOfReserve ; i++){
        if(aaKmerBuffer->buffer[i] != 0){
            startIdx = i;
            break;
        }
    }
    // Make splits
    vector<Split> splits;
    size_t splitWidth = (aaKmerBuffer->startIndexOfReserve - startIdx) / par.threads;
    for (int i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < aaKmerBuffer->startIndexOfReserve; j++) {
            if ((aaKmerBuffer->buffer[j] & mask_getId) != (aaKmerBuffer->buffer[j+1] & mask_getId)) {
                splits.emplace_back(startIdx, j);
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, aaKmerBuffer->startIndexOfReserve - 1);

    size_t ** idxOfEachSplit = new size_t * [splits.size()];
    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++) {
        idxOfEachSplit[i] = new size_t[splits[i].end - splits[i].offset + 2];
        cntOfEachSplit[i] = 0;
    }

#pragma omp parallel default(none), shared(aaKmerBuffer, idxOfEachSplit, cntOfEachSplit, splits, par)
    {
        uint64_t lookingKmer;
        size_t lookingIndex;
        int endFlag;
#pragma omp for schedule(dynamic, 1)
        for(size_t split = 0; split < splits.size(); split ++){
            lookingKmer = aaKmerBuffer->buffer[splits[split].offset];
            lookingIndex = splits[split].offset;
            endFlag = 0;
            for(size_t i = 1 + splits[split].offset; i < splits[split].end + 1 ; i++) {
                while((lookingKmer & mask_getId) == (aaKmerBuffer->buffer[i] & mask_getId)){
                    if ((lookingKmer & mask_getAA) != (aaKmerBuffer->buffer[i] & mask_getAA)) {
                        break;
                    }
                    i++;
                    if(i == splits[split].end + 1){
                        endFlag = 1;
                        break;
                    }
                }               
                idxOfEachSplit[split][cntOfEachSplit[split]++] = lookingIndex;
                if(endFlag == 1) break;
                lookingKmer = aaKmerBuffer->buffer[i];
                lookingIndex = i;
            }
            //For the end part
            if (aaKmerBuffer->buffer[splits[split].end - 1] != aaKmerBuffer->buffer[splits[split].end]) {
                idxOfEachSplit[split][cntOfEachSplit[split]++] = splits[split].end;
            }
        }
    }

    // Merge
    for(size_t i = 0; i < splits.size(); i++){
        memcpy(uniqeKmerIdx + uniqueKmerCnt, idxOfEachSplit[i], cntOfEachSplit[i] * sizeof(size_t));
        uniqueKmerCnt += cntOfEachSplit[i];
    }

    for(size_t i = 0; i < splits.size(); i++) { delete[] idxOfEachSplit[i];}
    delete[] idxOfEachSplit;
    delete[] cntOfEachSplit;
}


void ProteinDbIndexer::writeProteinIndexFile(const size_t * uniqKmerIdx, size_t uniqKmerCnt, bool completed) {
    string idxFileName;
    if (completed) {
        idxFileName = dbDir + "/protein_expanded.mtbl";
    } else {
        idxFileName = dbDir + "/" + to_string(flushCnt) + "_protein_expanded.mtbl";
    }
    FILE *fp = fopen(idxFileName.c_str(), "wb");
    flushCnt++;

    // Make splits
    FILE * idxSplitFile;
    ProtIdxSplit splitList[par.splitNum];    
    size_t splitCnt = 1;
    if (completed) {
        idxSplitFile = fopen(proteinIndexSplitFileName.c_str(), "wb");
        memset(splitList, 0, sizeof(ProtIdxSplit) * par.splitNum);
        size_t splitWidth = uniqKmerCnt / par.splitNum;
        size_t remainder = uniqKmerCnt % par.splitNum;
        size_t start = 0;
        splitList[0].kmer = aaKmerBuffer->buffer[uniqKmerIdx[0]];
        splitList[0].idxOffset = 0;
        for (size_t i = 1; i < (size_t) par.splitNum; i++) {
            start = start + splitWidth;
            if (remainder > 0) {
                start++;
                remainder--;
            }
            for (size_t j = start; j + 1 < start + splitWidth; j++) {
                if ((aaKmerBuffer->buffer[uniqKmerIdx[j]] & mask_getAA) != (aaKmerBuffer->buffer[uniqKmerIdx[j + 1]] & mask_getAA)) {
                    splitList[splitCnt].kmer = aaKmerBuffer->buffer[uniqKmerIdx[j + 1]];
                    splitCnt++;
                    break;
                }
            }
        }
    }
    
    if (fp == NULL) {
        cerr << "Error: Cannot open file " << idxFileName << endl;
        exit(1);
    }
    
    size_t splitIdx = 0;
    if (!completed) {
        for(size_t i = 0; i < uniqKmerCnt ; i++) {
            fwrite(&aaKmerBuffer->buffer[uniqKmerIdx[i]], sizeof (uint64_t), 1, fp);
        }
    } else {
        for(size_t i = 0; i < uniqKmerCnt ; i++) {
            fwrite(&this->aaKmerBuffer->buffer[uniqKmerIdx[i]], sizeof (uint64_t), 1, fp);
            if ((splitIdx < splitCnt) && (aaKmerBuffer->buffer[uniqKmerIdx[i]] == splitList[splitIdx].kmer)) {
                splitList[splitIdx].idxOffset = i;
                splitIdx++;
            }
        }
        fwrite(splitList, sizeof(ProtIdxSplit), par.splitNum, idxSplitFile);
    }
    // Write number of kmers
    fclose(fp);
}

void ProteinDbIndexer::writeProteinDeltaIndexFile(const size_t * uniqKmerIdx, size_t uniqKmerCnt, bool completed) {
    string deltaIdxFileName;
    string kmerNumFileName;
    FILE * kmerNum_fp;
    if (completed) {
        deltaIdxFileName = dbDir + "/protein.mtbl";
    } else {
        deltaIdxFileName = dbDir + "/" + to_string(flushCnt) + "_protein.mtbl";
        kmerNumFileName = dbDir + "/" + to_string(flushCnt) + "_protein_kmer_num.mtbl";
        kmerNum_fp = fopen(kmerNumFileName.c_str(), "wb");
        if (kmerNum_fp == NULL) {
            cerr << "Error: Cannot open file " << kmerNumFileName << endl;
            exit(1);
        }
    }
    FILE *delta_fp = fopen(deltaIdxFileName.c_str(), "wb");
    if (delta_fp == NULL) {
        cerr << "Error: Cannot open file " << deltaIdxFileName << endl;
        exit(1);
    }
    flushCnt++;

    // Temporary offsets
    size_t splitWidth = uniqKmerCnt / par.splitNum;
    size_t remainder = uniqKmerCnt % par.splitNum;
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

    // Offset list
    FILE * offsetList_fp;
    ProtIdxSplit * offsetList = new ProtIdxSplit[par.splitNum + 1];

    
    uint16_t *deltaIdxBuffer = (uint16_t *)malloc(sizeof(uint16_t) * this->bufferSize);
    size_t localBufIdx = 0;
    size_t totalBufferIdx = 0;
    uint64_t lastKmer = 0;
    if (!completed) {
        for(size_t i = 0; i < uniqKmerCnt ; i++) {
            IndexCreator::getDiffIdx(lastKmer,
             aaKmerBuffer->buffer[uniqKmerIdx[i]],
              delta_fp,
               deltaIdxBuffer,
                this->bufferSize,
                 localBufIdx);
            lastKmer = aaKmerBuffer->buffer[uniqKmerIdx[i]];
        }
    } else {
        offsetList_fp = fopen(proteinIndexSplitFileName.c_str(), "wb");
        if (offsetList_fp == NULL) {
            cerr << "Error: Cannot open file " << proteinIndexSplitFileName << endl;
            exit(1);
        }
        size_t processedKmers = 0;
        size_t offsetListIdx = 0;
        size_t processedOffset = 0;
        bool check = false;
        for(size_t i = 0; i < uniqKmerCnt ; i++) {
            const uint64_t & currentKmer = aaKmerBuffer->buffer[uniqKmerIdx[i]];
            IndexCreator::getDiffIdx(lastKmer,
             currentKmer,
              delta_fp,
               deltaIdxBuffer,
                this->bufferSize, 
                localBufIdx,
                 totalBufferIdx);

            // Check if temporary offset is reached
            if (processedKmers == tempOffsetList[offsetListIdx]) {
                if (processedKmers == 0) {
                    offsetList[processedOffset++] = {currentKmer, totalBufferIdx};
                    offsetListIdx++;
                } else {
                    check = true;
                    offsetListIdx++;
                }
            }
            processedKmers++;

            // If new temp. offset is reached and a new AA is seen, make a new offset point
            if (check && (lastKmer & mask_getAA) != (currentKmer & mask_getAA)) {
                offsetList[processedOffset++] = {currentKmer, totalBufferIdx};
                check = false;
            }
            lastKmer = currentKmer;
        }
        fwrite(offsetList, sizeof(ProtIdxSplit), par.splitNum, offsetList_fp);
    }
    IndexCreator::flushKmerBuf(deltaIdxBuffer, delta_fp, localBufIdx);
    if (!completed) {
        fwrite(&uniqKmerCnt, sizeof(size_t), 1, kmerNum_fp);
        fclose(kmerNum_fp);
    }
    free(deltaIdxBuffer);
    fclose(delta_fp);
    if (completed) {
        fclose(offsetList_fp);
    }
    delete[] tempOffsetList;
    delete[] offsetList;
}

void ProteinDbIndexer::mergeProteinIndexFiles() {
    string mergedFileName = dbDir + "/protein.mtbl";
    FILE * mergedFile = fopen(mergedFileName.c_str(), "wb");
    FILE * idxSplitFile = fopen(proteinIndexSplitFileName.c_str(), "wb");

    // Buffer
    uint64_t * kmerBuffer = (uint64_t *)malloc(sizeof(uint64_t) * bufferSize);
    size_t kmerBufferIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    uint64_t * lookingKmers = new uint64_t[flushCnt];
    auto * maxIdxOfEachFiles = new size_t[flushCnt];
    struct MmapedData<uint64_t> * protIdxFileList = new struct MmapedData<uint64_t>[flushCnt];
    size_t * currPositions = new size_t[flushCnt];
    memset(currPositions, 0, flushCnt * sizeof(size_t));

    for (int file = 0; file < flushCnt; file++) {
        protIdxFileList[file] = mmapData<uint64_t>((dbDir + "/" + to_string(file) + "_protein.mtbl").c_str());
        maxIdxOfEachFiles[file] = protIdxFileList[file].fileSize / sizeof(uint64_t);
        numOfKmerBeforeMerge += maxIdxOfEachFiles[file];
    }    
    
    // Temporary offsets
    size_t splitWidth = numOfKmerBeforeMerge / par.splitNum;
    size_t remainder = numOfKmerBeforeMerge % par.splitNum;
    size_t offsetList[par.splitNum + 1];
    offsetList[0] = 0;
    for (size_t i = 1; i < (size_t) par.splitNum; i++) {
        offsetList[i] = offsetList[i - 1] + splitWidth;
        if (remainder > 0) {
            offsetList[i]++;
            remainder--;
        }
    }
    offsetList[par.splitNum] = UINT64_MAX;

    ProtIdxSplit splitList[par.splitNum];
    memset(splitList, 0, sizeof(ProtIdxSplit) * par.splitNum);

    // Get the first k-mer of each ProtSplitFile
    for(size_t file = 0; file < flushCnt; file++){
        lookingKmers[file] =  protIdxFileList[file].data[0];
        currPositions[file] = 1;
    }

    size_t remainedFiles = flushCnt;
    size_t splitCnt = 0;
    uint64_t lastWrittenKmer = 0;
    size_t offsetListIdx = 0;
    bool check = false;
    size_t totalIdx = 0;
    while(true) {
        // Find the smallest k-mer
        size_t idxOfMin = getSmallestKmer(lookingKmers, flushCnt);
        uint64_t smallestKmer = lookingKmers[idxOfMin];        
        // Copy the smallest k-mer to the buffer
        // If the buffer is full, flush it 
        if (kmerBufferIdx == bufferSize) {
            fwrite(kmerBuffer, sizeof(uint64_t), kmerBufferIdx, mergedFile);
            kmerBufferIdx = 0;
        }
        kmerBuffer[kmerBufferIdx] = smallestKmer;
        
        // Check if offset is reached
        if (totalIdx == offsetList[offsetListIdx]) {
            if (totalIdx == 0) {
                splitList[splitCnt++] = {smallestKmer, totalIdx};
                offsetListIdx++;
            } else {
                check = true;
                offsetListIdx++;
            }
        }

        // If new offset is reached and new AA is seen, it is time to split
        if (check && (lastWrittenKmer & mask_getAA) != (smallestKmer & mask_getAA)) {
            splitList[splitCnt++] = {smallestKmer, totalIdx};
            check = false;
        }
        kmerBufferIdx++;
        totalIdx++;

        // Update looking k-mers
        if (currPositions[idxOfMin] == maxIdxOfEachFiles[idxOfMin]) {
            lookingKmers[idxOfMin] = UINT64_MAX;
            remainedFiles--;
            if (remainedFiles == 0) break;
        } else {
            lookingKmers[idxOfMin] = protIdxFileList[idxOfMin].data[currPositions[idxOfMin]++];
        }
    }

    fwrite(kmerBuffer, sizeof(uint64_t), kmerBufferIdx, mergedFile);
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, idxSplitFile);
    free(kmerBuffer);
    fclose(mergedFile);
    fclose(idxSplitFile);
 
    for(int file = 0; file < flushCnt; file++){
        munmap(protIdxFileList[file].data, protIdxFileList[file].fileSize + 1);
    }

    delete[] protIdxFileList;
    delete[] lookingKmers;
    delete[] currPositions;
    delete[] maxIdxOfEachFiles;
}

void ProteinDbIndexer::mergeProteinDeltaIndexFiles() {
    string mergedFileName = dbDir + "/protein.mtbl";
    FILE * mergedFile = fopen(mergedFileName.c_str(), "wb");
    FILE * idxSplitFile = fopen(proteinIndexSplitFileName.c_str(), "wb");

    // Buffer
    uint16_t * kmerBuffer = (uint16_t *)malloc(sizeof(uint16_t) * bufferSize);
    size_t kmerBufferIdx = 0;
    size_t totalBufferIdx = 0;

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    uint64_t * lookingKmers = new uint64_t[flushCnt];
    auto * maxIdxOfEachFiles = new size_t[flushCnt];
    struct MmapedData<uint16_t> * protIdxFileList = new struct MmapedData<uint16_t>[flushCnt];
    struct MmapedData<size_t> * kmerCntFileList = new struct MmapedData<size_t>[flushCnt];
    size_t * currPositions = new size_t[flushCnt];
    memset(currPositions, 0, flushCnt * sizeof(size_t));

    for (int file = 0; file < flushCnt; file++) {
        protIdxFileList[file] = mmapData<uint16_t>((dbDir + "/" + to_string(file) + "_protein.mtbl").c_str());
        kmerCntFileList[file] = mmapData<size_t>((dbDir + "/" + to_string(file) + "_protein_kmer_num.mtbl").c_str());
        maxIdxOfEachFiles[file] = protIdxFileList[file].fileSize / sizeof(uint16_t);
        numOfKmerBeforeMerge += kmerCntFileList[file].data[0];
        cout << "Num of kmers: " << kmerCntFileList[file].data[0] << endl;
    }    
    
    // Temporary offsets
    size_t splitWidth = numOfKmerBeforeMerge / par.splitNum;
    size_t remainder = numOfKmerBeforeMerge % par.splitNum;
    size_t offsetList[par.splitNum + 1];
    offsetList[0] = 0;
    for (size_t i = 1; i < (size_t) par.splitNum; i++) {
        offsetList[i] = offsetList[i - 1] + splitWidth;
        if (remainder > 0) {
            offsetList[i]++;
            remainder--;
        }
    }
    offsetList[par.splitNum] = UINT64_MAX;

    // Index splits to fill in
    ProtIdxSplit splitList[par.splitNum];
    memset(splitList, 0, sizeof(ProtIdxSplit) * par.splitNum);

    // Get the first k-mer of each ProtSplitFile
    for(int file = 0; file < flushCnt; file++){
        lookingKmers[file] =  FileMerger::getNextKmer(0, protIdxFileList[file], currPositions[file]);
    }

    size_t remainedFiles = flushCnt;
    size_t splitCnt = 0;
    uint64_t lastWrittenKmer = 0;
    size_t offsetListIdx = 0;
    bool check = false;
    size_t processedKmers = 0;
    while(true) {
        // Find the smallest k-mer
        size_t idxOfMin = getSmallestKmer(lookingKmers, flushCnt);
        uint64_t smallestKmer = lookingKmers[idxOfMin];
        IndexCreator::getDiffIdx(lastWrittenKmer,
                             smallestKmer,
                              mergedFile,
                               kmerBuffer,
                               bufferSize,
                                kmerBufferIdx,
                                 totalBufferIdx); 
        
        // Check if offset is reached
        if (processedKmers == offsetList[offsetListIdx]) {
            if (processedKmers == 0) {
                splitList[splitCnt++] = {smallestKmer, totalBufferIdx};
                offsetListIdx++;
            } else {
                check = true;
                offsetListIdx++;
            }
        }
        processedKmers++;

        // If new offset is reached and new AA is seen, it is time to split
        if (check && (lastWrittenKmer & mask_getAA) != (smallestKmer & mask_getAA)) {
            splitList[splitCnt++] = {smallestKmer, totalBufferIdx};
            check = false;
        }
        

        // Update looking k-mers
        if (currPositions[idxOfMin] == maxIdxOfEachFiles[idxOfMin]) {
            lookingKmers[idxOfMin] = UINT64_MAX;
            remainedFiles--;
            if(remainedFiles == 0) break;
        } else {
            lookingKmers[idxOfMin] =  FileMerger::getNextKmer(smallestKmer, protIdxFileList[idxOfMin], currPositions[idxOfMin]);
        }
        lastWrittenKmer = smallestKmer;
    }

    // flush buffer
    IndexCreator::flushKmerBuf(kmerBuffer, mergedFile, kmerBufferIdx);

    // write split
    fwrite(splitList, sizeof(ProtIdxSplit), par.splitNum, idxSplitFile);
    free(kmerBuffer);
    fclose(mergedFile);
    fclose(idxSplitFile);
 
    for(int file = 0; file < flushCnt; file++){
        munmap(protIdxFileList[file].data, protIdxFileList[file].fileSize + 1);
    }

    delete[] protIdxFileList;
    delete[] lookingKmers;
    delete[] currPositions;
    delete[] maxIdxOfEachFiles;
    delete[] kmerCntFileList;
}

void ProteinDbIndexer::splitFasta(const std::string &fastaFileName, std::vector<SequenceBlock> &blocks) {
  KSeqWrapper* kseq = KSeqFactory(fastaFileName.c_str());
  size_t idx = 0;
  while (kseq->ReadEntry()) {
    const KSeqWrapper::KSeqEntry & e = kseq->entry;
    this->proteinId2Index[e.name.s] = idx;
    proteinIndex2taxonomyId[idx++] = getTaxIdFromComment(kseq->entry.comment.s);
    blocks.emplace_back(e.headerOffset - 1,
                        e.headerOffset + e.sequence.l + e.newlineCount + e.sequenceOffset - e.headerOffset - 2,
                        e.sequence.l + e.newlineCount + e.sequenceOffset - e.headerOffset,
                        e.sequence.l);
  }
  delete kseq;
}

TaxID ProteinDbIndexer::getTaxIdFromComment(const std::string &header) {
// comment format: 
// peptidylprolyl isomerase n=1 Tax=Triplophysa tibetana TaxID=1572043 RepID=A0A5A9P0L4_9TELE
    size_t taxIDPos = header.find("TaxID=");
    size_t taxIDEndPos = header.find(" ", taxIDPos);
    return std::stoi(header.substr(taxIDPos + 6, taxIDEndPos - taxIDPos - 6));    
}

size_t ProteinDbIndexer::getSmallestKmer(const uint64_t lookingKmers[], size_t fileCnt) {
    size_t idxOfMin = 0;
    uint64_t minKmer = lookingKmers[0] & mask_getAA;
    uint64_t minId = lookingKmers[0] & mask_getId;
    for(size_t i = 0; i < fileCnt; i++) {
        size_t currentKmer = lookingKmers[i] & mask_getAA;
        size_t currentId = lookingKmers[i] & mask_getId;
        if (currentKmer < minKmer) {
            minKmer = currentKmer;
            minId = currentId;
            idxOfMin = i;
        } else if ((currentKmer == minKmer) && (currentId < minId)) {
            minId = currentId;
            idxOfMin = i;
        } 
    }
    return idxOfMin;
}

void ProteinDbIndexer::writePrtIdMap() {
    string prtIdMapFileName = dbDir + "/prtIdMap.mtbl";
    ofstream prtIdMapFile(prtIdMapFileName);
    for (auto & kv : proteinId2Index) {
        prtIdMapFile << kv.first << " " << kv.second << endl;
    }
    prtIdMapFile.close();
}

// void ProteinDbIndexer::inspectTaxonomy() {
//     size_t speciesCnt = 0;
//     size_t leafCnt = 0;
//     size_t genusOrHigherCnt = 0;
//     for (auto & kv : proteinIndex2taxonomyId) {
//         TaxID taxID = kv.second;
//         const TaxonNode * node = taxonomy->taxonNode(taxID);
//         taxonomy->getString(node->rankIdx)

//         // Check if it is a leaf


//     }
// }