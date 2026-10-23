#ifndef ADKMER4_INDEXCREATOR_H
#define ADKMER4_INDEXCREATOR_H
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <unordered_set>
#include <atomic>
#include <cstdint>
#include <iomanip>
#include <regex>
#ifdef OPENMP
    #include <omp.h>
#endif


#include "printBinary.h"
#include "Mmap.h"
#include "Kmer.h"
#include "SeqIterator.h"
#include "TaxonomyWrapper.h"
#include "BitManipulateMacros.h"
#include "common.h"
#include "FastSort.h"
#include "LocalParameters.h"
#include "NucleotideMatrix.h"
#include "SubstitutionMatrix.h"
#include "tantan.h"
#include "LocalUtil.h"
#include "TaxonomyWrapper.h"
#include "fasta_validate.h"
#include "GeneticCode.h"
#include "KmerExtractor.h"
#include "DeltaIdxReader.h"
#include "InfoIndex.h"
#include "UnirefTree.h"
#include "MetamerPattern.h"


enum class FilterMode { DB_CREATION, COMMON_KMER, UNIQ_KMER, UNIREF_LCA, DB_CREATION_POS};
struct Accession {
    Accession() = default;
    Accession(const string & accession, uint32_t whichFasta, uint32_t order, uint32_t length) 
        : accession(accession), whichFasta(whichFasta), order(order), length(length), speciesID(0), taxID(0) {}
    string accession;
    uint32_t whichFasta;
    uint32_t order;
    uint32_t length;
    TaxID speciesID;
    TaxID taxID;

    void print() {
        std::cout << accession << " " << whichFasta << " " << order << " " << length << " " << speciesID << " " << taxID << std::endl;
    }

    bool operator < (const Accession & a) const {
        if (speciesID != a.speciesID)
            return speciesID < a.speciesID;
        
        if (whichFasta != a.whichFasta)
            return whichFasta < a.whichFasta;

        if (order != a.order)
            return order < a.order;
        
        return false;
    }

    static bool compare(const Accession & a, const Accession & b) {
        if (a.speciesID != b.speciesID)
            return a.speciesID < b.speciesID;
        
        if (a.whichFasta != b.whichFasta)
            return a.whichFasta < b.whichFasta;

        if (a.order != b.order)
            return a.order < b.order;
            
        return false;
    }
};

struct AccessionBatch {
    uint32_t whichFasta;
    TaxID speciesID;
    uint32_t trainingSeqFasta;
    uint32_t trainingSeqIdx;
    uint64_t totalLength;
    vector<uint32_t> orders;
    vector<TaxID> taxIDs;
    vector<uint32_t> lengths;

    void print() const {
        std::cout << "whichFasta: " << whichFasta << " speciesID: " << speciesID << " trainingSeqFasta: " << trainingSeqFasta << " trainingSeqIdx: " << trainingSeqIdx << endl;
        for (size_t i = 0; i < orders.size(); ++i) {
            std::cout << "order: " << orders[i] << " taxID: " << taxIDs[i] << " length: " << lengths[i] << endl;
        }
    }

    AccessionBatch(uint32_t whichFasta, TaxID speciesID, uint32_t trainingSeqFasta, uint32_t trainingSeqIdx, uint64_t totalLength) 
        : whichFasta(whichFasta), speciesID(speciesID), trainingSeqFasta(trainingSeqFasta), trainingSeqIdx(trainingSeqIdx), totalLength(totalLength) {}
};

struct FastaBatch {
    TaxID taxId;
    uint32_t whichFasta;
    uint64_t length;
    uint64_t n50;
    uint64_t l50;
    vector<AccessionBatch> accessionBatches;

    FastaBatch(uint32_t whichFasta) : whichFasta(whichFasta), length(0), n50(0), l50(0) {}
};

struct SpeciesBatch {
    TaxID speciesID;
    uint64_t spTotalLength;
    uint64_t expectedKmerNum;
    uint32_t repGenomeFasta;
    uint64_t repGenomeSize;
    vector<FastaBatch> fastaBatches;

    SpeciesBatch(TaxID speciesID) : speciesID(speciesID), spTotalLength(0), expectedKmerNum(0), repGenomeFasta(0), repGenomeSize(0) {}
};

using namespace std;

class IndexCreator{
protected:
    // Parameters
    const LocalParameters & par;
    bool isUpdating;
    int kmerFormat;
    int kmerLen;
    int windowSize;

    uint64_t MARKER;
    BaseMatrix *subMat;
    bool removeRedundancyInfo;
    unordered_map<TaxID, TaxID> taxId2speciesId;
    std::vector<TaxID> taxaNotToMask;
    std::vector<TaxID> taxaToSkipProdigal;

    // Inputs
    TaxonomyWrapper * taxonomy = nullptr;
    UnirefTree * unirefTree = nullptr;
    MetamerPattern * metamerPattern = nullptr;
    GeneticCode * geneticCode = nullptr;
    KmerExtractor * kmerExtractor;

    bool externTaxonomy;
    string dbDir;
    string fnaListFileName;
    string acc2taxidFileName;

    // Outputs
    string taxidListFileName;
    string taxonomyBinaryFileName;
    string versionFileName;
    string paramterFileName;

    std::unordered_map<string, vector<CDSinfo>> cdsInfoMap;
    std::vector<AccessionBatch> accessionBatches;
    std::vector<SpeciesBatch> spBatches;

    // Per-flush chunk files (<n>_diffIdx/<n>_info/<n>_kmerpos) are merge inputs
    // only; once the merged diffIdx/info/split are written they are dead weight
    // (hundreds of GB on large builds), so mergeTargetFiles() removes them.
    // External merge inputs (e.g., updateDB's old database) live in other
    // directories under different names and are never touched.
    void removeFlushChunks() {
        for (size_t i = 0; i < numOfFlush; ++i) {
            const std::string prefix = dbDir + "/" + std::to_string(i) + "_";
            remove((prefix + "diffIdx").c_str());
            remove((prefix + "info").c_str());
            remove((prefix + "kmerpos").c_str()); // absent unless --kmer-pos; remove() just fails silently
        }
    }

    // --mask 3 (auto): tantan for eukaryotic species (catches the long-period
    // tandem repeats SDUST cannot see; its strand asymmetry is harmless for
    // single-representative euk genomes), SDUST for everything else (strand-
    // symmetric, so multi-strain species keep consistent k-mer sets across
    // genomes and sub-DBs).
    int effectiveMaskMode(TaxID speciesID) const {
        if (par.maskMode != 3) {
            return par.maskMode;
        }
        if (taxonomy->getEukaryotaTaxID() != 0
            && taxonomy->IsAncestor(taxonomy->getEukaryotaTaxID(), speciesID)) {
            return 1;
        }
        return 2;
    }
    // Per-batch multiplier on the k-mer-count reservation estimate. Grown (x2)
    // when a batch's actual k-mer count exceeds its reservation, so the batch is
    // re-processed with more room instead of overrunning the shared buffer.
    std::vector<float> batchEstimateScale;
    std::unordered_set<TaxID> taxIdSet;
    std::vector<std::string> fastaPaths;
    std::vector<std::string> unusedFastaPaths;
    size_t numOfFlush=0;
    size_t totalLength=0;

    // Database splits
    std::vector<std::string> deltaIdxFileNames;
    std::vector<std::string> infoFileNames;
    std::vector<std::string> posFileNames;
    std::string mergedDeltaIdxFileName;
    std::string mergedInfoFileName;
    std::string mergedPosFileName;
    std::string deltaIdxSplitFileName;
    InfoIndexMetadata finalInfoMetadata;
    struct Split{
        Split(size_t offset, size_t end) : offset(offset), end(end) {}
        size_t offset;
        size_t end;

        void print() {
            std::cout << "offset: " << offset << " end: " << end << std::endl;
        }
    };

    void writeTargetFiles(
        Buffer<Kmer> & kmerBuffer,
        const size_t * uniqeKmerIdx,
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges);
    
    void writeTargetFiles(
        Buffer<Kmer> & kmerBuffer);

    void writeTargetFilesAndSplits(
        Buffer<Kmer> & kmerBuffer,
        const size_t * uniqeKmerIdx,
        size_t & uniqKmerCnt,
        const vector<pair<size_t, size_t>> & uniqKmerIdxRanges,
        bool streamPackTaxIds = false);

    void writeTargetFilesAndSplits(
        Buffer<Kmer> & kmerBuffer);

    void writeDbParameters();
    void writeInfoMetadata();
    void finalizeInfoIndex(const std::string &infoFileName, uint32_t maxInfoId, size_t idCount);

    // Bit width for the final info index, decided up front so tax-ID databases
    // can stream packed straight to disk. Honors --pack-info 0 (stays uint32).
    uint8_t chooseInfoIdBits() const;
    // Record the info-index format in db.parameters after a streaming write
    // (no repack pass needed, unlike finalizeInfoIndex).
    void recordInfoMetadata(uint8_t idBits, size_t idCount);

    size_t fillTargetKmerBuffer(
        Buffer<Kmer> &kmerBuffer,                 
        std::vector<std::atomic<bool>> & batchChecker,
        size_t &processedSplitCnt,
        const LocalParameters &par);

    size_t fillTargetKmerBuffer2(
        Buffer<Kmer> &kmerBuffer,                 
        std::vector<std::atomic<bool>> & batchChecker,
        size_t &processedSpCnt,
        const LocalParameters &par);
    
    bool extractKmerFromSixFrames(
        Buffer<Kmer> & kmerBuffer,
        std::vector<std::atomic<bool>> & batchChecker,
        size_t &processedBatchCnt);

    void indexReferenceSequences(size_t bufferSize);

    void getAccessionBatches(std::vector<Accession> & observedAccessionsVec, size_t bufferSize);

    void getSpeciesBatches();

    void printSpeciesBatches();

    void getObservedAccessions(const string & fnaListFileName,
                               vector<Accession> & observedAccessionsVec,
                               unordered_map<string, size_t> & accession2index);

    void getTaxonomyOfAccessions(vector<Accession> & observedAccessionsVec,
                                 const unordered_map<string, size_t> & accession2index,
                                 const string & acc2taxidFileName);

    void unzipAndList(const string & folder, const string & fastaList_fname){
        system(("./../../util/unzip_and_list.sh " + folder + " " + fastaList_fname).c_str());
    }

    void load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid);

    // TaxID getMaxTaxID();

    // void editTaxonomyDumpFiles(const vector<pair<string, pair<TaxID, TaxID>>> & newAcc2taxid);

    template <FilterMode M>
    bool areKmersDuplicate(
        const Kmer & kmer1,
        const Kmer & kmer2
    );

    size_t AminoAcidPart(size_t kmer) {
        if (kmerFormat == 3 || kmerFormat == 4) {
            return kmer;
        }
        return (kmer) & MARKER;
    }

    void loadCdsInfo(const string & cdsInfoFileList);

    size_t calculateBufferSize(size_t maxRam, size_t perItem) {
        constexpr double GIGABYTE = 1024.0 * 1024.0 * 1024.0;
        constexpr double MEGABYTE = 1024.0 * 1024.0;
        constexpr double MEMORY_PER_THREAD_MB = 50.0;
        
        float c = 0.7f;
        if (maxRam < 16) {
            c = 0.5f;
        } else if (maxRam <= 32) {
            c = 0.6f;
        }
        double totalRamBytes = maxRam * GIGABYTE * c;
        double memoryForThreads = par.threads * MEMORY_PER_THREAD_MB * MEGABYTE;
        double availableMemory = totalRamBytes - memoryForThreads;

        if (availableMemory <= 0.0) {
            cerr << "Not enough memory to create index" << endl;
            cerr << "Please increase the RAM usage or decrease the number of threads" << endl;
            exit(EXIT_FAILURE);
        }

        return static_cast<size_t>(availableMemory / perItem);
    }

    size_t calculateBufferSizeForMerge(size_t maxRam, int fileCnt) {
        float c = 0.7;
        if (maxRam <= 32) {
            c = 0.6;
        } else if (maxRam < 16) {
            c = 0.5;
        }
        if ((maxRam * 1024.0 * 1024.0 * 1024.0 * c - (par.threads * 50.0 * 1024.0 * 1024.0)) <= 0.0) {
            cerr << "Not enough memory to create index" << endl;
            cerr << "Please increase the RAM usage or decrease the number of threads" << endl;
            exit(EXIT_FAILURE);
        }
        return static_cast<size_t>((maxRam * 1024.0 * 1024.0 * 1024.0 * c - (par.threads * 50.0 * 1024.0 * 1024.0))/ 
                                  (sizeof(Kmer) + sizeof(size_t)));
    }

    void loadMergedTaxIds(const std::string &mergedFile, unordered_map<TaxID, TaxID> & old2new);

    string addToLibrary(const std::string & dbDir,
                        const std::string & fileList,
                        const std::string & acc2taxIdFileName);

    void generateIntergenicKmerList(
        struct _gene *genes, 
        struct _node *nodes, 
        int numberOfGenes,
        vector<uint64_t> &intergenicKmerList, 
        const char *seq);

    void makeIntergenicKmerList(
        const std::string & fastaFileName,
        vector<uint64_t> & intergenicKmers,
        ProdigalWrapper * prodigal);

    bool compareMinHashList(
        const std::vector<uint64_t>& list1, 
        const std::vector<uint64_t>& list2, 
        size_t length1, 
        size_t length2);
    
    std::vector<uint64_t> getMinHashList(const char* seq);

    void devideToCdsAndNonCds(
        size_t seqLen,
        const vector<CDSinfo> &cdsInfo, 
        std::vector<SequenceBlock> & fragments);

    void recordSpeciesGenomeSize() {
        string fileName = dbDir + "/" + "species2genomeSize.tsv";
        ofstream out(fileName);
        if (!out.is_open()) {
            cerr << "Failed to open file for writing: " << fileName << endl;
            return;
        }
        for (const auto & spBatch : spBatches) {
            out << spBatch.speciesID << "\t" << taxonomy->getOriginalTaxID(spBatch.speciesID) << "\t" << spBatch.repGenomeSize << endl;
        }
        out.close();
    }
    



public:
    IndexCreator(const LocalParameters & par, TaxonomyWrapper * taxonomy, int kmerFormat);
    IndexCreator(const LocalParameters & par, UnirefTree * unirefTree, int kmerFormat);
    IndexCreator(const LocalParameters & par, int kmerFormat); // Used in create_unique_kmer_list.cpp
    ~IndexCreator();
    void createIndex();
    void createIndexWithPos();
    void createCommonKmerIndex();
    void createUniqueKmerIndex();
    void createLcaKmerIndex();

    template <FilterMode M>
    void mergeTargetFiles();

    template <FilterMode M>
    void filterKmers(
        Buffer<Kmer> & kmerBuffer,
        size_t * uniqeKmerIdx,
        size_t & uniqKmerCnt,
        vector<pair<size_t, size_t>> & uniqKmerIdxRanges);

    static void getDiffIdx(
        uint64_t & lastKmer,
        uint64_t entryToWrite,
        WriteBuffer<uint16_t> & diffBuffer);

    // Getters
    int getNumOfFlush() const { return numOfFlush; }
    TaxonomyWrapper* getTaxonomy() const { return taxonomy;}
    unordered_set<TaxID> getTaxIdSet() { return taxIdSet; }

    // Setters
    void setIsUpdating(bool isUpdating) { this->isUpdating = isUpdating; }
    void setIsNewFormat(int kmerFormat) { this->kmerFormat = kmerFormat; }
    void addFilesToMerge(string diffIdxFileName, string infoFileName, string posFileName = "");
    void updateTaxId2SpeciesTaxId(const string & taxIdListFileName);
    void setMergedFileNames(string diffFileName, string infoFileName, string splitFileName, string posFileName = "");

    void printFilesToMerge() {
        cout << "Files to merge :" << endl;
        for (size_t i = 0; i < deltaIdxFileNames.size(); i++) {
            cout << deltaIdxFileNames[i] << " " << infoFileNames[i];
            if (!posFileNames.empty()) {
                cout << " " << posFileNames[i];
            }
            cout << endl;
        }
    }
    
    static void printIndexSplitList(DiffIdxSplit * splitList) {
        for (int i = 0; i < 4096; i++) {
            std::cout << splitList[i].infoIdxOffset << " " << 
                    splitList[i].diffIdxOffset << " " << 
                    splitList[i].ADkmer << std::endl;
        }
    }
};

// --- Merge-write helpers (parallel delta/info encoding) ----------------------
// Fragment count of one 15-bit varint delta; matches getDiffIdx() exactly.
static inline size_t deltaFragmentCount(uint64_t delta) {
    if (delta == 0) {
        return 1;
    }
    const int bits = 64 - __builtin_clzll(delta);
    return static_cast<size_t>((bits + 14) / 15);
}

// Encode one delta as big-fragment-first 15-bit groups with the end flag on
// the last fragment; byte-identical to getDiffIdx(). Returns fragment count.
static inline size_t encodeDeltaFragments(uint64_t delta, uint16_t * out) {
    const size_t n = deltaFragmentCount(delta);
    for (size_t f = 0; f + 1 < n; ++f) {
        out[f] = static_cast<uint16_t>((delta >> (15 * (n - 1 - f))) & 0x7FFF);
    }
    out[n - 1] = static_cast<uint16_t>(0x8000u | (delta & 0x7FFF));
    return n;
}

template <FilterMode M>
void IndexCreator::mergeTargetFiles() {
    size_t bufferSize = 1024 * 1024 * 512;
    WriteBuffer<uint16_t> diffBuffer(mergedDeltaIdxFileName, bufferSize);
    // Tax-ID databases (DB_CREATION / DB_CREATION_POS) stream the info index in
    // its final packed encoding, chosen up front from the taxonomy's max ID, so
    // no unpacked file is written and repacked. Other modes write uint32 and let
    // finalizeInfoIndex() repack afterwards from the observed max ID.
    constexpr bool streamPackTaxIds =
        (M == FilterMode::DB_CREATION || M == FilterMode::DB_CREATION_POS);
    const uint8_t infoIdBits = streamPackTaxIds ? chooseInfoIdBits() : 32;
    InfoWriter infoBuffer(mergedInfoFileName, infoIdBits, bufferSize);
    WriteBuffer<uint16_t> posBuffer(mergedPosFileName, bufferSize);

    // Prepare files to merge
    size_t numOfKmerBeforeMerge = 0;
    size_t splitNum = deltaIdxFileNames.size();
    DeltaIdxReader ** deltaIdxReaders = new DeltaIdxReader*[splitNum];
    size_t valueBufferSize = 1024 * 1024 * 16;
    for (size_t file = 0; file < splitNum; file++) {
        if (posBuffer.buffer == nullptr)
        {
            deltaIdxReaders[file] = new DeltaIdxReader(deltaIdxFileNames[file],
                                                       infoFileNames[file],
                                                       valueBufferSize,
                                                       1024 * 1024 * 4);
        }
        else
        {
            deltaIdxReaders[file] = new DeltaIdxReader(deltaIdxFileNames[file],
                                                       infoFileNames[file],
                                                       posFileNames[file],
                                                       valueBufferSize,
                                                       1024 * 1024 * 4);
        }
        numOfKmerBeforeMerge += deltaIdxReaders[file]->getTotalValueNum();
    }

    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = numOfKmerBeforeMerge / (par.splitNum - 1);
    size_t offsetList[par.splitNum + 1];
    int offsetListIdx = 1;
    for(int os = 0; os < par.splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }
    offsetList[par.splitNum] = UINT64_MAX;
    DiffIdxSplit splitList[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    int splitListIdx = 1;

    // get the first k-mer to write
    unsigned int mask = ~((static_cast<unsigned int>(par.skipRedundancy == 0) << 31));
    int splitCheck = 0;

    cout << "Merging target files..." << endl;
    size_t kmerBufferSize = 1024 * 1024 * 1024;
    if (kmerBufferSize < valueBufferSize * splitNum * 2) {
        kmerBufferSize = valueBufferSize * splitNum * 2; 
    }
    Buffer<Kmer> kmerBuffer(kmerBufferSize); // 64GB
    std::atomic<int> hasOverflow{0};
    std::vector<std::atomic<bool>> completedSplits(splitNum);
    int remainingSplits = splitNum;
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    uint64_t lastKmer = 0;
    // DB_CREATION_POS streams k-mers straight through without the uniq index,
    // so skip its large allocation in that mode.
    size_t * uniqKmerIdx = nullptr;
    if (M != FilterMode::DB_CREATION_POS) {
        uniqKmerIdx = new size_t[kmerBuffer.bufferSize];
    }
    vector<size_t> splitToProcess;
    // Round-persistent scratch for the parallel write phase (grown once).
    vector<size_t> rangePrefix;
    vector<size_t> chunkFrag;
    vector<uint32_t> chunkMax;
    vector<uint16_t> roundDiff;
    vector<uint64_t> roundWords;
    vector<size_t> fragCkpt;
    while (remainingSplits > 0) {
        kmerBuffer.init();
        
        time_t start = time(nullptr);
        while (remainingSplits > 0
                && kmerBuffer.startIndexOfReserve + valueBufferSize * remainingSplits <= kmerBuffer.bufferSize) {
            size_t posToWrite = kmerBuffer.reserveMemory(valueBufferSize * remainingSplits);
            uint64_t max = UINT64_MAX;
            splitToProcess.clear();
            for (size_t i = 0; i < splitNum; i++) {
                if (completedSplits[i].load(std::memory_order_acquire))
                    continue; 
                if (deltaIdxReaders[i]->getLastValue() < max) {
                    max = deltaIdxReaders[i]->getLastValue();
                }
                splitToProcess.push_back(i);
            }
            #pragma omp parallel for default(none), shared(cout, kmerBuffer, deltaIdxReaders, splitToProcess, completedSplits, posToWrite, max, splitNum, valueBufferSize, taxId2speciesId, mask, remainingSplits)
            for (size_t i = 0; i < splitToProcess.size(); i ++) {
                size_t split = splitToProcess[i];
                size_t offset = posToWrite + i * valueBufferSize;
                size_t valueNum = deltaIdxReaders[split]->getValues(kmerBuffer.buffer + offset, max);
                if (deltaIdxReaders[split]->isCompleted()) {
                    completedSplits[split].store(true, std::memory_order_release);
                    __sync_fetch_and_sub(&remainingSplits, 1);
                }
                if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::COMMON_KMER) {            
                    for (size_t j = 0; j < valueNum; j++) {
                        kmerBuffer.buffer[offset + j].tInfo.speciesId 
                            = taxId2speciesId[kmerBuffer.buffer[offset + j].tInfo.taxId & mask];
                    }
                }
            }
        }
        time_t end = time(nullptr);
        cout << "K-mer loading       : " << (double) (end - start) << " s" << endl;

        time_t beforeSort = time(nullptr);

        if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::COMMON_KMER || M == FilterMode::DB_CREATION_POS) {
            SORT_PARALLEL(kmerBuffer.buffer, 
                          kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                          Kmer::compareTargetKmer);
        } else if constexpr (M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
            SORT_PARALLEL(kmerBuffer.buffer, 
                          kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                          Kmer::compareKmer);
        }
        
        time_t afterSort = time(nullptr);
        cout << "Sorting k-mer list  : " << afterSort - beforeSort << " s" << endl;

        // Reduce redundancy
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        time_t reduction = time(nullptr);
        if (M == FilterMode::DB_CREATION_POS) {
            for (size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++) {
                if (kmerBuffer.buffer[i].isEmpty()) {
                    continue;
                }

                if (kmerBuffer.buffer[i].value == UINT64_MAX) {
                    break;
                }

                uint16_t pos = static_cast<uint16_t>(kmerBuffer.buffer[i].tInfo.pos);
                posBuffer.write(&pos);
                infoBuffer.write(kmerBuffer.buffer[i].id);
                getDiffIdx(lastKmer, kmerBuffer.buffer[i].value, diffBuffer);
                if (AminoAcidPart(lastKmer) != AAofTempSplitOffset && splitCheck == 1) {
                    splitList[splitListIdx++] = {lastKmer, diffBuffer.writeCnt, infoBuffer.writeCnt};
                    splitCheck = 0;
                }
                if (infoBuffer.writeCnt == offsetList[offsetListIdx]) {
                    AAofTempSplitOffset = AminoAcidPart(lastKmer);
                    splitCheck = 1;
                    offsetListIdx++;
                }
            }
        } else {
            // filterKmers overwrites every slot of the index ranges it reports,
            // so uniqKmerIdx needs no zeroing (an 8 GB memset per round).
            filterKmers<M>(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
            reduction = time(nullptr);
            cout << "Filtering k-mers    : " << (double) (reduction - afterSort) << " s" << endl;
            cout << "Selected count      : " << uniqKmerCnt << endl;

            // ---- Write the selected k-mers (parallel chunked encode) ----
            // The output streams only look serial: a chunk's first delta needs
            // just the previous SELECTED k-mer's value (random-accessible), and
            // a packed-info bit position is a pure function of the global
            // ordinal. So chunks encode independently into one round buffer at
            // exact offsets from a counting pass, byte-identical to the old
            // one-at-a-time loop.
            const size_t nRanges = uniqKmerIdxRanges.size();
            rangePrefix.assign(nRanges + 1, 0);
            for (size_t r = 0; r < nRanges; ++r) {
                rangePrefix[r + 1] = rangePrefix[r] +
                    (uniqKmerIdxRanges[r].second - uniqKmerIdxRanges[r].first);
            }
            const size_t totalSel = rangePrefix[nRanges];
            if (totalSel > 0) {
                const Kmer * kb = kmerBuffer.buffer;
                // j-th selected k-mer (sparse accesses only).
                auto selAt = [&](size_t j) -> const Kmer & {
                    const size_t r = (size_t)(std::upper_bound(rangePrefix.begin(),
                        rangePrefix.end(), j) - rangePrefix.begin()) - 1;
                    return kb[uniqKmerIdx[uniqKmerIdxRanges[r].first + (j - rangePrefix[r])]];
                };
                // Sequential walker over the selected indices (ranges may be empty).
                struct SelWalk {
                    const size_t * idx;
                    const vector<pair<size_t, size_t>> * ranges;
                    size_t r = 0;
                    size_t pos = 0;
                    void seek(const vector<size_t> & prefix, size_t j) {
                        r = (size_t)(std::upper_bound(prefix.begin(), prefix.end(), j)
                                     - prefix.begin()) - 1;
                        pos = (*ranges)[r].first + (j - prefix[r]);
                    }
                    size_t next() {
                        while (pos >= (*ranges)[r].second) {
                            ++r;
                            pos = (*ranges)[r].first;
                        }
                        return idx[pos++];
                    }
                };

                const uint8_t vpw = InfoIndex::idsPerWord(infoIdBits);
                const size_t baseOrd = infoBuffer.writeCnt;
                const size_t diffBase = diffBuffer.writeCnt;
                const uint64_t roundStartKmer = lastKmer;
                const size_t firstWord = baseOrd / vpw;
                const size_t nWords = (baseOrd + totalSel + vpw - 1) / vpw - firstWord;

                size_t nChunks = (size_t) par.threads;
                const size_t maxChunks = totalSel / 4096 + 1;
                if (nChunks > maxChunks) {
                    nChunks = maxChunks;
                }

                roundWords.resize(nWords);
                #pragma omp parallel for schedule(static)
                for (size_t w = 0; w < nWords; ++w) {
                    roundWords[w] = 0;
                }

                // Pass 1: per-chunk fragment totals and ID maxima.
                chunkFrag.assign(nChunks + 1, 0);
                chunkMax.assign(nChunks, 0);
                #pragma omp parallel for schedule(static, 1)
                for (size_t c = 0; c < nChunks; ++c) {
                    const size_t j0 = totalSel * c / nChunks;
                    const size_t j1 = totalSel * (c + 1) / nChunks;
                    if (j0 >= j1) {
                        continue;
                    }
                    SelWalk walk{uniqKmerIdx, &uniqKmerIdxRanges};
                    walk.seek(rangePrefix, j0);
                    uint64_t prev = (j0 == 0) ? roundStartKmer : selAt(j0 - 1).value;
                    size_t frags = 0;
                    uint32_t mx = 0;
                    for (size_t j = j0; j < j1; ++j) {
                        const Kmer & k = kb[walk.next()];
                        frags += deltaFragmentCount(k.value - prev);
                        prev = k.value;
                        if (k.id > mx) {
                            mx = k.id;
                        }
                    }
                    chunkFrag[c + 1] = frags;
                    chunkMax[c] = mx;
                }
                uint32_t roundMax = 0;
                for (size_t c = 0; c < nChunks; ++c) {
                    chunkFrag[c + 1] += chunkFrag[c];
                    if (chunkMax[c] > roundMax) {
                        roundMax = chunkMax[c];
                    }
                }
                const size_t roundFrags = chunkFrag[nChunks];
                roundDiff.resize(roundFrags);
                fragCkpt.assign((totalSel >> 16) + 1, 0);

                // Pass 2: encode deltas and OR-pack IDs at their final offsets.
                #pragma omp parallel for schedule(static, 1)
                for (size_t c = 0; c < nChunks; ++c) {
                    const size_t j0 = totalSel * c / nChunks;
                    const size_t j1 = totalSel * (c + 1) / nChunks;
                    if (j0 >= j1) {
                        continue;
                    }
                    SelWalk walk{uniqKmerIdx, &uniqKmerIdxRanges};
                    walk.seek(rangePrefix, j0);
                    uint64_t prev = (j0 == 0) ? roundStartKmer : selAt(j0 - 1).value;
                    size_t fragOff = chunkFrag[c];
                    // Words at chunk boundaries are shared with a neighbor and
                    // need atomic merging; interior words are exclusive.
                    const size_t cFirstW = (baseOrd + j0) / vpw - firstWord;
                    const size_t cLastW = (baseOrd + j1 - 1) / vpw - firstWord;
                    uint64_t curWord = 0;
                    size_t curW = SIZE_MAX;
                    for (size_t j = j0; j < j1; ++j) {
                        if ((j & 0xFFFF) == 0) {
                            fragCkpt[j >> 16] = fragOff; // for sparse F() lookups
                        }
                        const Kmer & k = kb[walk.next()];
                        fragOff += encodeDeltaFragments(k.value - prev, roundDiff.data() + fragOff);
                        prev = k.value;
                        const size_t ord = baseOrd + j;
                        const size_t w = ord / vpw - firstWord;
                        if (w != curW) {
                            if (curW != SIZE_MAX) {
                                if (curW == cFirstW || curW == cLastW) {
                                    __atomic_or_fetch(&roundWords[curW], curWord, __ATOMIC_RELAXED);
                                } else {
                                    roundWords[curW] = curWord;
                                }
                            }
                            curW = w;
                            curWord = 0;
                        }
                        curWord |= static_cast<uint64_t>(k.id) << ((ord % vpw) * infoIdBits);
                    }
                    // The final word is always a boundary word.
                    __atomic_or_fetch(&roundWords[curW], curWord, __ATOMIC_RELAXED);
                }

                // Fragment-count prefix of the first k selected k-mers, from the
                // nearest checkpoint (used only around the few split anchors).
                auto fragPrefixAt = [&](size_t k) -> size_t { // k >= 1
                    size_t j = ((k - 1) >> 16) << 16;
                    size_t frags = fragCkpt[(k - 1) >> 16];
                    uint64_t prev = (j == 0) ? roundStartKmer : selAt(j - 1).value;
                    for (; j < k; ++j) {
                        const uint64_t v = selAt(j).value;
                        frags += deltaFragmentCount(v - prev);
                        prev = v;
                    }
                    return frags;
                };

                // Split anchors: replay the serial state machine, but touch only
                // the k-mers around each threshold (anchors are rare). The
                // serial loop checks the anchor BEFORE the threshold on the same
                // k-mer, so the anchor scan includes the threshold position.
                size_t scanJ = 0;
                while (true) {
                    size_t thrJ = totalSel; // next threshold inside this round
                    if (offsetList[offsetListIdx] >= baseOrd + 1 &&
                        offsetList[offsetListIdx] <= baseOrd + totalSel) {
                        thrJ = offsetList[offsetListIdx] - baseOrd - 1;
                    }
                    if (splitCheck == 1) {
                        bool found = false;
                        for (size_t s = scanJ; s < totalSel && s <= thrJ; ++s) {
                            const uint64_t v = selAt(s).value;
                            if (AminoAcidPart(v) != AAofTempSplitOffset) {
                                splitList[splitListIdx++] =
                                    {v, diffBase + fragPrefixAt(s + 1), baseOrd + s + 1};
                                splitCheck = 0;
                                found = true;
                                break;
                            }
                        }
                        if (!found && thrJ >= totalSel) {
                            break; // pending anchor carries over to the next round
                        }
                    }
                    if (thrJ >= totalSel) {
                        break;
                    }
                    AAofTempSplitOffset = AminoAcidPart(selAt(thrJ).value);
                    splitCheck = 1;
                    offsetListIdx++;
                    scanJ = thrJ + 1;
                }

                diffBuffer.writeBulk(roundDiff.data(), roundFrags);
                infoBuffer.writeWords(roundWords.data(), nWords, totalSel, roundMax);
                lastKmer = selAt(totalSel - 1).value;
            }
        }

  
        time_t writeTime = time(nullptr);
        cout << "Writing k-mers      : " << (double) (writeTime - reduction) << " s" << endl;
        cout << "Written k-mers      : " << infoBuffer.writeCnt << " " << std::fixed << std::setprecision(2) << ((float) infoBuffer.writeCnt / numOfKmerBeforeMerge) * 100 << "%" << endl;
        cout << "--------------------" << endl;
        cout << "remaining splits    : " << remainingSplits << endl;
        std::cout.unsetf(std::ios::fixed | std::ios::scientific);
        kmerBuffer.startIndexOfReserve = 0;
    }

    FILE * diffIdxSplitFile = fopen(deltaIdxSplitFileName.c_str(), "wb");
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, diffIdxSplitFile);
    fclose(diffIdxSplitFile);
    const size_t finalInfoCount = infoBuffer.writeCnt;
    infoBuffer.close();
    // Tax-ID databases already wrote the final packed encoding, so only the
    // metadata needs recording. Other modes still repack the uint32 stream from
    // the observed max ID.
    if constexpr (streamPackTaxIds) {
        recordInfoMetadata(infoIdBits, finalInfoCount);
    } else {
        finalizeInfoIndex(mergedInfoFileName, infoBuffer.maxId, finalInfoCount);
    }
    // for(int i = 0; i < par.splitNum; i++) {
    //     cout<<splitList[i].ADkmer<< " "<<splitList[i].diffIdxOffset<< " "<<splitList[i].infoIdxOffset<<endl;
    // }
    cout<<"DB creation completed"<<endl;
    cout<<"Total k-mer count   : " << finalInfoCount <<endl;

    cout<<"DB files you need   : " << endl;
    cout<<mergedDeltaIdxFileName<<endl;
    cout<<mergedInfoFileName<<endl;
    cout<<string(dbDir) + "/taxID_list"<<endl;
    cout<<deltaIdxSplitFileName<<endl;

    removeFlushChunks();
}


template <FilterMode M>
void IndexCreator::filterKmers(
    Buffer<Kmer> & kmerBuffer,
    size_t * selectedKmerIdx,
    size_t & selectedKmerCnt,
    vector<pair<size_t, size_t>> & selectedKmerIdxRanges) 
{
    // Find the first index of garbage k-mer (UINT64_MAX)
    for(size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
        if(kmerBuffer.buffer[checkN].value != UINT64_MAX){
            kmerBuffer.startIndexOfReserve = checkN + 1;
            break;
        }
    }

    // Find the first index of meaningful k-mer
    size_t startIdx = 0;
    for (size_t i = 0; i < kmerBuffer.startIndexOfReserve ; i++) {
        if(!kmerBuffer.buffer[i].isEmpty()){
            startIdx = i;
            break;
        }
    }
    cout << "Loaded k-mer count  : " << (kmerBuffer.startIndexOfReserve - startIdx) << endl;
    // Make splits
    vector<Split> splits;
    size_t splitWidth = (kmerBuffer.startIndexOfReserve - startIdx) / par.threads;
    for (int i = 0; i < par.threads - 1; i++) {
        for (size_t j = startIdx + splitWidth; j + 1 < kmerBuffer.startIndexOfReserve; j++) {
            if (AminoAcidPart(kmerBuffer.buffer[j].value) != AminoAcidPart(kmerBuffer.buffer[j + 1].value)) {
                splits.emplace_back(startIdx, j);
                selectedKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));
                startIdx = j + 1;
                break;
            }
        }
    }
    splits.emplace_back(startIdx, kmerBuffer.startIndexOfReserve - 1);
    selectedKmerIdxRanges.emplace_back(pair<size_t, size_t>(startIdx, 0));

    size_t * cntOfEachSplit = new size_t[splits.size()];
    for(size_t i = 0; i < splits.size(); i++) {
        cntOfEachSplit[i] = 0;
    }
#pragma omp parallel default(none), shared(kmerBuffer, cntOfEachSplit, splits, par, cout, selectedKmerCnt, selectedKmerIdx, selectedKmerIdxRanges)
    {
        Kmer * lookingKmer;
        size_t lookingIndex;
        vector<TaxID> taxIds;
        vector<uint32_t> ids;
        size_t * localSelectedIdx = new size_t[16 * 1024 * 1024];
        size_t tempSelectedKmerCnt = 0;
#pragma omp for schedule(static, 1)
        for(size_t split = 0; split < splits.size(); split ++) {
            size_t i = splits[split].offset;
            i++;
            while (i < splits[split].end + 1) {
                lookingKmer = & kmerBuffer.buffer[i - 1];
                lookingIndex = i - 1;
                bool selected = false;
                if constexpr (M == FilterMode::DB_CREATION) {
                    taxIds.clear();
                    taxIds.push_back(lookingKmer->tInfo.taxId);
                } else if constexpr (M == FilterMode::COMMON_KMER) {
                    taxIds.clear();
                    taxIds.push_back(lookingKmer->tInfo.speciesId);
                } else if constexpr (M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
                    ids.clear();
                    ids.push_back(lookingKmer->id);
                }

                while ((i < splits[split].end + 1) && areKmersDuplicate<M>(*lookingKmer, kmerBuffer.buffer[i])) {
                    if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::DB_CREATION_POS) {
                        taxIds.push_back(kmerBuffer.buffer[i].tInfo.taxId);
                    } else if constexpr (M == FilterMode::COMMON_KMER) {
                        taxIds.push_back(kmerBuffer.buffer[i].tInfo.speciesId);
                    } else if constexpr (M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
                        ids.push_back(kmerBuffer.buffer[i].id);
                    }
                    i++;
                }
                
                if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::UNIREF_LCA || M == FilterMode::DB_CREATION_POS) {
                    selected = true;
                } else if constexpr (M == FilterMode::COMMON_KMER) {
                    for (size_t i = 0; i < taxIds.size(); i++) {
                        if (taxIds[i] != lookingKmer->tInfo.speciesId) {
                            selected = true;
                            break;
                        }
                    }
                } else if constexpr (M == FilterMode::UNIQ_KMER) {
                    selected = true;
                    for (size_t i = 1; i < ids.size(); i++) {
                        if (ids[i] != lookingKmer->id) {
                            selected = false;
                            break;
                        }
                    }
                }

                if (selected) {
                    if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::COMMON_KMER || M == FilterMode::DB_CREATION_POS) {
                        lookingKmer->tInfo.taxId = taxonomy->LCA(taxIds)->taxId;
                    } else if constexpr (M == FilterMode::UNIREF_LCA) {
                        lookingKmer->id = unirefTree->getLCA(ids);
                    }

                    if (tempSelectedKmerCnt >= 16 * 1024 * 1024) {
                        memcpy(selectedKmerIdx + splits[split].offset, localSelectedIdx, tempSelectedKmerCnt * sizeof(size_t));
                        splits[split].offset += tempSelectedKmerCnt;
                        tempSelectedKmerCnt = 0;
                    }
                    localSelectedIdx[tempSelectedKmerCnt++] = lookingIndex;
                    cntOfEachSplit[split]++;
                }

                i++;                
            }

            // Check the last k-mer
            if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA || M == FilterMode::DB_CREATION_POS) {
                if(!areKmersDuplicate<M>(kmerBuffer.buffer[splits[split].end - 1], kmerBuffer.buffer[splits[split].end])){
                    if (tempSelectedKmerCnt >= 16 * 1024 * 1024) {
                        memcpy(selectedKmerIdx + splits[split].offset, localSelectedIdx, tempSelectedKmerCnt * sizeof(size_t));
                        splits[split].offset += tempSelectedKmerCnt;
                        tempSelectedKmerCnt = 0;
                    }
                    localSelectedIdx[tempSelectedKmerCnt++] = splits[split].end;
                    cntOfEachSplit[split]++;
                }
            }
            memcpy(selectedKmerIdx + splits[split].offset, localSelectedIdx, tempSelectedKmerCnt * sizeof(size_t));
            splits[split].offset += tempSelectedKmerCnt;
            selectedKmerIdxRanges[split].second = splits[split].offset;
            __sync_fetch_and_add(&selectedKmerCnt, cntOfEachSplit[split]); 
        }
        delete[] localSelectedIdx;
    }
    delete[] cntOfEachSplit;
}

template <FilterMode M>
bool IndexCreator::areKmersDuplicate(
    const Kmer & kmer1,
    const Kmer & kmer2) 
{
    if constexpr (M == FilterMode::DB_CREATION || M == FilterMode::DB_CREATION_POS) {
        return kmer1.tInfo.speciesId == kmer2.tInfo.speciesId &&
               kmer1.value == kmer2.value;
    } else if constexpr (M == FilterMode::COMMON_KMER || M == FilterMode::UNIQ_KMER || M == FilterMode::UNIREF_LCA) {
        return kmer1.value == kmer2.value;
    }
    return false; // Default case, should not be reached
}

#endif
