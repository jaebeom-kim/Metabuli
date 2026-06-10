#include "IndexCreator.h"
#include "FileUtil.h"
#include "LocalUtil.h"
#include "ProdigalWrapper.h"
#include <cstdint>
#include <cstdio>
#include <unordered_map>
#include <utility>
#include "NcbiTaxonomy.cpp"
#include "common.h"

extern const char *version;

IndexCreator::IndexCreator(
    const LocalParameters & par, 
    TaxonomyWrapper * taxonomy,
    int kmerFormat) 
    : par(par), taxonomy(taxonomy), kmerFormat(kmerFormat) 
{
    dbDir = par.filenames[0];
    fnaListFileName = par.filenames[1];
    if (par.filenames.size() >= 3) {
        acc2taxidFileName = par.filenames[2];
    }
    
    taxidListFileName = dbDir + "/taxID_list";
    taxonomyBinaryFileName = dbDir + "/taxonomyDB";
    versionFileName = dbDir + "/db.version";
    paramterFileName = dbDir + "/db.parameters";

    this->totalLength = par.dbTotalLength;

    if (kmerFormat == 1) { // Use the legacy metamer pattern
        metamerPattern = new LegacyPattern(std::make_unique<RegularGeneticCode>(), 8);
    } else if (!par.customMetamer.empty()) {
        // Use custom metamer pattern
        int codeNum = getCodeNum(par.customMetamer);
        int weightedPosNum = getWeightedPosNum(par.customMetamer);
        if (codeNum == 1) {
            if (par.spaceMask == "") {
                metamerPattern = new SingleCodePattern(par.customMetamer);
            } else {
                uint32_t mask = parseMask(par.spaceMask.c_str());
                int setBitNum = __builtin_popcount(mask);
                if (setBitNum != weightedPosNum) {
                    std::cerr << "Error: The number of set bits in space mask (" << setBitNum 
                              << ") does not match length (" 
                              << weightedPosNum << ") in custom metamer pattern." << std::endl;
                    exit(EXIT_FAILURE);
                }
                metamerPattern = new SpacedPattern(par.customMetamer, mask);
            }
        } else if (codeNum > 1) {
            if (par.spaceMask != "") {
                cout << "Error: Spaced k-mer isn't supported with a multi-code custom metamer pattern." << endl;
                exit(EXIT_FAILURE);
            }
            metamerPattern = new MultiCodePattern(par.customMetamer);
        }
    } else {
        if (par.spaceMask.empty()) {
            metamerPattern = new SingleCodePattern(std::make_unique<RegularGeneticCode>(), 8);
        } else {
            uint32_t mask = parseMask(par.spaceMask.c_str());
            metamerPattern = new SpacedPattern(std::make_unique<RegularGeneticCode>(), __builtin_popcount(mask), mask);
        }
    }
    kmerExtractor = new KmerExtractor(par, metamerPattern);
    MARKER = ~ metamerPattern->dnaMask;
    isUpdating = false;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);

    if (!par.noMaskTaxa.empty()) { 
        vector<string> taxaNotToMaskStr = Util::split(par.noMaskTaxa, ",");
        for (const string &taxIdStr : taxaNotToMaskStr) {
            TaxID taxId = taxonomy->getInternalTaxID(stoi(taxIdStr));
            taxaNotToMask.push_back(taxId);
        }
    }

}

IndexCreator::IndexCreator(
    const LocalParameters & par, 
    UnirefTree * unirefTree,
    int kmerFormat) 
    : par(par), kmerFormat(kmerFormat), unirefTree(unirefTree)
{
    dbDir = par.filenames[0];
    // if (par.reducedAA) {
    //     geneticCode = new ReducedGeneticCode();
    // } else {
        geneticCode = new RegularGeneticCode();
    // }
    kmerExtractor = new KmerExtractor(par, geneticCode, kmerFormat);
    isUpdating = false;
}

IndexCreator::IndexCreator(
    const LocalParameters & par, 
    int kmerFormat) : par(par), kmerFormat(kmerFormat) 
{
    dbDir = par.filenames[0];
    // if (par.reducedAA) {
    //     geneticCode = new ReducedGeneticCode();
    // } else {
        geneticCode = new RegularGeneticCode();
    // }
    kmerExtractor = new KmerExtractor(par, geneticCode, kmerFormat);
    isUpdating = false;
    subMat = new NucleotideMatrix(par.scoringMatrixFile.values.nucleotide().c_str(), 1.0, 0.0);
}


IndexCreator::~IndexCreator() {
    delete metamerPattern;
    if (geneticCode != nullptr) {
        delete geneticCode;
    }   
    delete kmerExtractor;
    delete subMat;
}

TaxID IndexCreator::getCollapseTaxId(TaxID taxId) const {
    TaxID collapseTaxId = taxonomy->getTaxIdAtRank(taxId, par.collapseRank);
    return collapseTaxId == 0 ? taxId : collapseTaxId;
}

void IndexCreator::createLcaKmerIndex() {
    Buffer<Kmer> kmerBuffer(Buffer<Kmer>::calculateBufferSize(par.ramUsage, par.threads, sizeof(Kmer) + sizeof(size_t)));
    Buffer<size_t> uniqKmerIdx(kmerBuffer.bufferSize);
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;

    string fileName = par.filenames[1];
    KSeqWrapper * kseq = KSeqFactory(fileName.c_str());
    std::unordered_map<string, uint32_t> name2id;

    std::cout << "Filling UniRef100 name to ID mapping ... " << std::endl;
    time_t start = time(nullptr);
    unirefTree->getName2Id(name2id);
    cout << "Mapping filled in " << time(nullptr) - start << " s." << endl;


    bool complete = false;
    SeqEntry savedSeq;
    uint32_t processedSeqCnt = 0;
    while (!complete) {
        // Extract k-mers
        time_t start = time(nullptr);
        cout << "K-mer extraction    : " << flush;
        bool moreData = kmerExtractor->extractUnirefKmers(kseq, kmerBuffer, name2id, processedSeqCnt, savedSeq);
        complete = !moreData;
        cout << double(time(nullptr) - start) << " s" << endl;
        cout << "Processed sequences : " << processedSeqCnt << endl;

        // Sort the k-mers
        start = time(nullptr);
        cout << "Sort k-mers         : " << flush;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
        cout << double(time(nullptr) - start) << " s" << endl;

        // Filter k-mers
        start = time(nullptr);
        size_t selectedKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<FilterMode::UNIREF_LCA>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        cout << "Filter k-mers       : " << time(nullptr) - start << " s" << endl; 
        cout << "Selected k-mers     : " << selectedKmerCnt << endl;

        // Write k-mers
        start = time(nullptr);
        if (complete && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
        }
        cout << "Write k-mers        : " << time(nullptr) - start << " s" << endl;

        if (!complete) {
            kmerBuffer.init();
            uniqKmerIdx.init();
        }
        cout << "--------" << endl;
    }
    
    // for (int i = 0; i < 47; i++) {
    //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
    //                     dbDir + "/" + to_string(i) + "_info");
    // }

    if (numOfFlush == 1) {
        cout << "Index creation completed." << endl;
        return;
    }
    cout << "Merge reference DB files ... " << endl;
    printFilesToMerge();
    setMergedFileNames(
        par.filenames[0] + "/diffIdx",  
        par.filenames[0] + "/info", 
        par.filenames[0] + "/split");
    mergeTargetFiles<FilterMode::UNIREF_LCA>();

}

void IndexCreator::createUniqueKmerIndex() {
    Buffer<Kmer> kmerBuffer(Buffer<Kmer>::calculateBufferSize(par.ramUsage, par.threads, sizeof(Kmer) + sizeof(size_t)));
    Buffer<size_t> uniqKmerIdx(kmerBuffer.bufferSize);
    bool complete = false;
    string fileName = par.filenames[1];
    KSeqWrapper * kseq = KSeqFactory(fileName.c_str());
    std::unordered_map<string, uint32_t> accession2index;
    uint32_t idOffset = 0;
    SeqEntry savedSeq;
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    while (!complete) {
        // Extract k-mers
        time_t start = time(nullptr);
        cout << "K-mer extraction    : " << flush;
        bool moreData = kmerExtractor->extractKmers(kseq, kmerBuffer, accession2index, idOffset, savedSeq);
        complete = !moreData;
        cout << double(time(nullptr) - start) << " s" << endl;
        cout << "Processed sequences : " << idOffset << endl;

        // Sort the k-mers
        start = time(nullptr);
        cout << "Sort k-mers         : " << flush;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareKmer);
        cout << double(time(nullptr) - start) << " s" << endl;

        // Filter k-mers
        
        start = time(nullptr);
        size_t selectedKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<FilterMode::UNIQ_KMER>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        cout << "Filter k-mers       : " << time(nullptr) - start << " s" << endl; 
        cout << "Selected k-mers     : " << selectedKmerCnt << endl;

        // Write k-mers
        
        start = time(nullptr);
        if (complete && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
        }
        cout << "Write k-mers        : " << time(nullptr) - start << " s" << endl;

        // Write accession to index mapping
        FILE * accIndexFile = fopen((dbDir + "/accession2index").c_str(), "a");
        for (const auto & entry : accession2index) {
            fprintf(accIndexFile, "%s\t%u\n", entry.first.c_str(), entry.second);
        }
        fclose(accIndexFile);

        if (!complete) {
            accession2index.clear();
            kmerBuffer.init();
            uniqKmerIdx.init();
        }
        cout << "--------" << endl;
    }
    

    if (numOfFlush == 1) {
        cout << "Index creation completed." << endl;
        return;
    }
    cout << "Merge reference DB files ... " << endl;

    // for (int i = 0; i < 10; i++) {
    //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
    //                     dbDir + "/" + to_string(i) + "_info");
    // }

    printFilesToMerge();
    setMergedFileNames(
        par.filenames[0] + "/diffIdx",  
        par.filenames[0] + "/info", 
        par.filenames[0] + "/split");
    mergeTargetFiles<FilterMode::UNIQ_KMER>();

}

void IndexCreator::createCommonKmerIndex() {
    size_t sizE = calculateBufferSize(par.ramUsage, sizeof(Kmer) + sizeof(size_t));
    indexReferenceSequences(sizE);
    Buffer<Kmer> kmerBuffer(sizE);
    Buffer<size_t> uniqKmerIdx(sizE);
    if (!par.cdsInfo.empty() && par.cdsInfo != "x") {
        cout << "Loading CDS info from: " << par.cdsInfo << endl;
        loadCdsInfo(par.cdsInfo);
    }

    // Write taxonomy ID list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    std::vector<std::atomic<bool>> batchChecker(accessionBatches.size());
    size_t processedBatchCnt = 0;
    
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
    while(processedBatchCnt < accessionBatches.size()) {
        // Extract target k-mers
        time_t start = time(nullptr);
        cout << "K-mer extraction : " << flush;
        if (par.cdsInfo == "x") {
            extractKmerFromSixFrames(kmerBuffer, batchChecker, processedBatchCnt);
        } else {
            fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);
        }
        cout << double(time(nullptr) - start) << " s" << endl;
        
        // Sort the k-mers
        start = time(nullptr);
        cout << "Sort k-mers      : " << flush;
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve, Kmer::compareTargetKmer);
        cout << time(nullptr) - start << " s" << endl;

        // Filter k-mers
        start = time(nullptr);
        size_t selectedKmerCnt = 0;
        filterKmers<FilterMode::DB_CREATION>(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        cout << "Filter k-mers    : " << time(nullptr) - start << " s" << endl; 
        cout << "Selected k-mers  : " << selectedKmerCnt << endl;

        // Write the target files
        start = time(nullptr);
        if(processedBatchCnt == accessionBatches.size() && numOfFlush == 0) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx.buffer, selectedKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx.buffer, uniqKmerIdxRanges);
        }
        cout << "Write k-mers     : " << time(nullptr) - start << " s" << endl;

        // Reset buffers
        if (processedBatchCnt < accessionBatches.size()) {
            kmerBuffer.init();
            uniqKmerIdx.init();
            uniqKmerIdxRanges.clear();
        }
    }

    taxonomy->writeTaxonomyDB(par.filenames[0] + "/taxonomyDB");
    writeDbParameters();
    
    if (numOfFlush == 1) {
        cout << "Index creation completed." << endl;
        return;
    }
    cout << "Merge reference DB files ... " << endl;

    // for (int i = 0; i < 243; i++) {
    //     addFilesToMerge(dbDir + "/" + to_string(i) + "_diffIdx",
    //                     dbDir + "/" + to_string(i) + "_info");
    // }
    // updateTaxId2SpeciesTaxId(dbDir + "/taxID_list");

    printFilesToMerge();
    setMergedFileNames(
        par.filenames[0] + "/diffIdx",  
        par.filenames[0] + "/info", 
        par.filenames[0] + "/split");
    mergeTargetFiles<FilterMode::COMMON_KMER>();
}


void IndexCreator::createIndexWithPos() {
    Buffer<Kmer> kmerBuffer(calculateBufferSize(par.ramUsage, sizeof(Kmer)));
    cout << "Target metamer buffer size: " << kmerBuffer.bufferSize << endl;

    indexReferenceSequences(kmerBuffer.bufferSize);
    getSpeciesBatches();
    recordSpeciesGenomeSize();
    printSpeciesBatches();
    cout << "Species batches prepared: " << spBatches.size() << " species." << endl;



    if (!par.cdsInfo.empty()) {
        cout << "Loading CDS info from: " << par.cdsInfo << endl;
        loadCdsInfo(par.cdsInfo);
    }

    // Write taxonomy id list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);


    // Process the splits until all are processed
    std::vector<std::atomic<bool>> batchChecker(spBatches.size());
    size_t processedSpCnt = 0;

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    while (processedSpCnt < spBatches.size()) {
        kmerBuffer.init();

        fillTargetKmerBuffer2(kmerBuffer, batchChecker, processedSpCnt, par);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      Kmer::compareKmer);
        time_t sort = time(nullptr);
        cout << "Reference k-mer sort : " << sort - start << " s" << endl;

        // Write database files
        if(processedSpCnt == spBatches.size() && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer);
        } else {
            writeTargetFiles(kmerBuffer);
        }
    }
    writeDbParameters();
}





















void IndexCreator::createIndex() {
    Buffer<Kmer> kmerBuffer(calculateBufferSize(par.ramUsage, sizeof(Kmer) + sizeof(size_t)));
    cout << "Target metamer buffer size: " << kmerBuffer.bufferSize << endl;
    
    indexReferenceSequences(kmerBuffer.bufferSize);
    cout << "Made blocks for each thread" << endl;

    if (!par.cdsInfo.empty()) {
        cout << "Loading CDS info from: " << par.cdsInfo << endl;
        loadCdsInfo(par.cdsInfo);
    }

    // Write taxonomy id list
    FILE * taxidListFile = fopen(taxidListFileName.c_str(), "w");
    for (auto & taxid: taxIdSet) {
        fprintf(taxidListFile, "%d\n", taxid);
    }
    fclose(taxidListFile);

    // Process the splits until all are processed
    std::vector<std::atomic<bool>> batchChecker(accessionBatches.size());
    size_t processedBatchCnt = 0;
    
    vector<pair<size_t, size_t>> uniqKmerIdxRanges;
#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif
    while(processedBatchCnt < accessionBatches.size()) {
        kmerBuffer.init();
        cout << "Buffer initialized" << endl;

        // Extract target k-mers
        fillTargetKmerBuffer(kmerBuffer, batchChecker, processedBatchCnt, par);

        // Sort the k-mers
        time_t start = time(nullptr);
        SORT_PARALLEL(kmerBuffer.buffer, kmerBuffer.buffer + kmerBuffer.startIndexOfReserve,
                      Kmer::compareTargetKmer);
        time_t sort = time(nullptr);
        cout << kmerBuffer.startIndexOfReserve << " k-mers extracted" << endl;
        cout << "Reference k-mer sort : " << sort - start << endl;

        // Reduce redundancy
        auto * uniqKmerIdx = new size_t[kmerBuffer.startIndexOfReserve + 1];
        size_t uniqKmerCnt = 0;
        uniqKmerIdxRanges.clear();
        filterKmers<FilterMode::DB_CREATION>(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        time_t reduction = time(nullptr);
        cout << "Unique k-mer count   : " << uniqKmerCnt << endl;
        cout << "Redundancy reduction : " << (double) (reduction - sort) << " s" << endl;

        // Write the target files
        if(processedBatchCnt == accessionBatches.size() && numOfFlush == 0 && !isUpdating) {
            writeTargetFilesAndSplits(kmerBuffer, uniqKmerIdx, uniqKmerCnt, uniqKmerIdxRanges);
        } else {
            writeTargetFiles(kmerBuffer, uniqKmerIdx, uniqKmerIdxRanges);
        }
        delete[] uniqKmerIdx;
    }

    writeDbParameters();
}

string IndexCreator::addToLibrary(
    const std::string & dbDir,
    const std::string & fnaListFileName,
    const std::string & acc2taxIdFileName) 
{
    // Make library directory
    time_t now = time(0);
    tm *ltm = localtime(&now);
    string timeStr = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday) + "-" + to_string(ltm->tm_hour) + "-" + to_string(ltm->tm_min);
    string libraryDir = dbDir + "/" + timeStr;
    if (!FileUtil::directoryExists(libraryDir.c_str())) {
        FileUtil::makeDir(libraryDir.c_str());
    }

    unordered_map<std::string, TaxID> accession2taxid;
    vector<std::string> fileNames;
    vector<string> libraryFiles;
    unordered_set<TaxID> observedSpecies;
    getObservedAccessionList(fnaListFileName, fileNames, accession2taxid);
    fillAcc2TaxIdMap(accession2taxid, acc2taxIdFileName);

    unordered_map<TaxID, TaxID> original2internalTaxId;
    if (taxonomy->hasInternalTaxID()) {
      taxonomy->getOriginal2InternalTaxId(original2internalTaxId);
    } 
    
    vector<std::string> unmapped;
    std::string accession;
    for (size_t i = 0; i < fileNames.size(); ++i) {
      KSeqWrapper* kseq = KSeqFactory(fileNames[i].c_str());
      while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry & e = kseq->entry;
        accession = string(e.name.s);
        size_t pos = accession.find('.');
        if (pos != std::string::npos) {
          accession = accession.substr(0, pos);
        }
        TaxID taxId = accession2taxid[accession];
        if (taxId == 0) {
          std::cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
               " is not found in the mapping file. It is skipped.\n";
          unmapped.push_back(e.name.s);
          continue;
        }
        if (taxonomy->hasInternalTaxID()) {
          if (original2internalTaxId.find(taxId) == original2internalTaxId.end()) {
            std::cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
            ", " << taxId << " is not included in the taxonomy. It is skipped.\n";
            unmapped.push_back(e.name.s);
            continue;
          }
          taxId = original2internalTaxId[taxId];
        }
        int speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
        if (speciesTaxID == 0) {
          cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
               " is not matched to any species. It is skipped.\n";
          unmapped.push_back(e.name.s);
          continue;
        }
        if (taxonomy->hasInternalTaxID()) {
            speciesTaxID = taxonomy->getOriginalTaxID(speciesTaxID);
        }
        string speciesFileName = libraryDir + "/" + to_string(speciesTaxID) + ".fna";
        if (observedSpecies.find(speciesTaxID) == observedSpecies.end()) {
            observedSpecies.insert(speciesTaxID);
            libraryFiles.push_back(speciesFileName);
        }
        FILE *file = fopen(speciesFileName.c_str(), "a");
        fprintf(file, ">%s %s\n", e.name.s, e.comment.s);
        fprintf(file, "%s\n", e.sequence.s);
        fclose(file);
      }
      delete kseq;
    }

    // Write unmapped accession to file
    if (unmapped.empty()) {
        cout << "All accessions are mapped to taxonomy" << endl;
    } else {
        FILE *file = fopen((libraryDir + "/unmapped.txt").c_str(), "w");
        for (const auto & i : unmapped) {
            fprintf(file, "%s\n", i.c_str());
        }
        fclose(file);
        cout << "Unmapped accessions are written to " << libraryDir + "/unmapped.txt" << endl;
    }  
    // Write the list of absolute paths of library files
    string libraryListFileName = libraryDir + "/library.list";
    FILE *libraryListFile = fopen(libraryListFileName.c_str(), "w");
    for (size_t i = 0; i < libraryFiles.size(); ++i) {
        fprintf(libraryListFile, "%s\n", libraryFiles[i].c_str());
    }
    fclose(libraryListFile);
    return libraryListFileName;
}

void IndexCreator::indexReferenceSequences(size_t bufferSize) {
    time_t start = time(nullptr);
    cout << "Make sequence batches for each thread : " << flush;
    
    vector<Accession> observedAccessionsVec;       // vector of observed accessions
    unordered_map<string, size_t> accession2index; // map accession to its index in the observedAccessionsVec 
    if (par.makeLibrary) {
        getObservedAccessions(addToLibrary(dbDir, fnaListFileName, acc2taxidFileName), observedAccessionsVec, accession2index);
    } else {
        getObservedAccessions(fnaListFileName, observedAccessionsVec, accession2index);
    }
    cout << "Number of observed accessions: " << observedAccessionsVec.size() << endl;
    getTaxonomyOfAccessions(observedAccessionsVec, accession2index, acc2taxidFileName);
    cout << "Taxonomy of accessions is obtained" << endl;
    vector<Accession> accessionsWithTaxonomy;
    getAccessionBatches(observedAccessionsVec, bufferSize);
    cout << "Number of accession batches: " << accessionBatches.size() << endl;

    if (par.storeKmerPos == 0) {
        sort(accessionBatches.begin(), accessionBatches.end(),
         [](const AccessionBatch &a, const AccessionBatch &b) {
             return a.totalLength < b.totalLength;
         });
    }


    time_t end = time(nullptr);
    cout << end - start << " s" << endl;
}


void IndexCreator::getObservedAccessions(
    const string & fnaListFileName,
    vector<Accession> & observedAccessionsVec,
    unordered_map<string, size_t> & accession2index) 
{
    ifstream fileListFile(fnaListFileName);
    if (fileListFile.is_open()) {
        for (string eachLine; getline(fileListFile, eachLine);) {
            fastaPaths.push_back(eachLine);
        }
    } else {
        cout << "Cannot open file for file list" << endl;
    } 

    // Iterate through the fasta files to get observed accessions
    size_t accCnt = 0;
    size_t copyCount = 0;
    #pragma omp parallel default(none), shared(accession2index, cout, observedAccessionsVec, accCnt, copyCount)
    {
        vector<Accession> localObservedAccessionsVec;
        localObservedAccessionsVec.reserve(4096 * 4);
        unordered_set<string> duplicateCheck;
        
        #pragma omp for schedule(static, 1)
        for (size_t i = 0; i < fastaPaths.size(); ++i) {
            KSeqWrapper* kseq = KSeqFactory(fastaPaths[i].c_str());
            uint32_t order = 0;
            while (kseq->ReadEntry()) {
                const KSeqWrapper::KSeqEntry & e = kseq->entry;
                // Get the accession ID without version
                char* pos = strchr(e.name.s, '.'); 
                if (pos != nullptr) {
                    *pos = '\0';
                }
                if (duplicateCheck.find(e.name.s) != duplicateCheck.end()) {
                    order++; 
                    continue;
                } else {
                    duplicateCheck.insert(e.name.s);
                } 
                localObservedAccessionsVec.emplace_back(string(e.name.s), i, order, e.sequence.l);
                order++; 
            }
            delete kseq;
        } 
        __sync_fetch_and_add(&accCnt, localObservedAccessionsVec.size()); 
        #pragma omp barrier
       
        #pragma omp critical
        {
            if (observedAccessionsVec.size() < accCnt) {
                observedAccessionsVec.resize(accCnt);
                accession2index.reserve(accCnt);
            }
        }   

        #pragma omp critical
        {
            size_t start = copyCount;
            for (size_t j = start ; j < start + localObservedAccessionsVec.size(); ++j) {
                observedAccessionsVec[j] = localObservedAccessionsVec[j - start];
                accession2index[observedAccessionsVec[j].accession] = j;
            }
            copyCount += localObservedAccessionsVec.size();
        }                     
    }
}

void IndexCreator::getTaxonomyOfAccessions(vector<Accession> & observedAccessionsVec,
                                           const unordered_map<string, size_t> & accession2index,
                                           const string & acc2taxidFileName) {
    unordered_map<TaxID, TaxID> old2merged;
    taxonomy->getMergedNodeMap(old2merged, true);
    
    vector<pair<string, pair<TaxID, TaxID>>> acc2accId;   
    int fd = open(acc2taxidFileName.c_str(), O_RDONLY);
    if (fd < 0) {
        cout << "Cannot open file for mapping from accession to tax ID" << endl;
        return;
    }

    // Get the size of the file
    struct stat sb;
    if (fstat(fd, &sb) == -1) {
        cout << "Cannot get the size of the file for mapping from accession to tax ID" << endl;
        close(fd);
        return;
    }

    size_t fileSize = sb.st_size;

    // Map the file to memory
    char* fileData = static_cast<char*>(mmap(nullptr, fileSize, PROT_READ, MAP_PRIVATE, fd, 0));
    if (fileData == MAP_FAILED) {
        cout << "mmap failed" << endl;
        close(fd);
        return;
    }
    close(fd);  // Close the file descriptor as it is no longer needed after mmap.

    // Parse the file
    char* current = fileData;
    char* end = fileData + fileSize;

    // Skip the header line
    while (current < end && *current != '\n') {
        ++current;
    }
    ++current;  // Move past the newline

    char accession[16384];
    TaxID taxID;
    std::unordered_set<TaxID> usedExternalTaxIDs;
    std::vector<NewTaxon> newTaxons;
    if (par.accessionLevel == 1) {
        taxonomy->getUsedExternalTaxIDs(usedExternalTaxIDs);
    }

    // First, label the accessions with external taxIDs
    while (current < end) {
        // Read a line
        char* lineStart = current;
        while (current < end && *current != '\n') {
            ++current;
        }
        std::string line(lineStart, current - lineStart);

        // Parse the line
        if (sscanf(line.c_str(), "%s\t%*s\t%d\t%*d", accession, &taxID) == 2) {
            // Get the accession ID without version
            // char* pos = strchr(accession, '.');
            // if (pos != nullptr) {
            //     *pos = '\0';
            // }
            auto it = accession2index.find(accession);
            if (it != accession2index.end()) {
                if (old2merged.count(taxID) > 0) {
                    taxID = old2merged[taxID];
                }
                if (par.accessionLevel == 1) {
                    TaxID accTaxId = taxonomy->getSmallestUnusedExternalTaxID(usedExternalTaxIDs);
                    acc2accId.emplace_back(accession, make_pair(taxID, accTaxId));
                    if (accTaxId == 0) {
                        cout << "accTaxId is 0 for accession " << accession << " " << taxID << endl;
                    }
                    observedAccessionsVec[it->second].taxID = accTaxId;
                    newTaxons.emplace_back(accTaxId, taxID, "accession", accession);
                } else {
                    observedAccessionsVec[it->second].taxID = taxID;
                }
            }
        }
        ++current;  // Move to the next line
    }

    if (munmap(fileData, fileSize) == -1) {
        cout << "munmap failed" << endl;
    }                 

    if (par.accessionLevel == 1) {
        TaxonomyWrapper * newTaxonomy = taxonomy->addNewTaxa(newTaxons);
        delete taxonomy;
        taxonomy = newTaxonomy;
    }

    // Second, convert external taxIDs to internal taxIDs
    cout << "Converting external taxIDs to internal taxIDs" << endl;
    vector<std::string> unmappedAccessions;
    std::unordered_map<TaxID, TaxID> external2internalTaxID;
    taxonomy->getExternal2internalTaxID(external2internalTaxID);
    std::unordered_set<TaxID> spIdSet;
    for (size_t i = 0; i < observedAccessionsVec.size(); ++i) {
        TaxID curId = 0;
        if (taxonomy->hasInternalTaxID()) {
            auto it = external2internalTaxID.find(observedAccessionsVec[i].taxID);
            if (it == external2internalTaxID.end() || observedAccessionsVec[i].taxID == 0) {
                cout << "TaxID is not found for accession " << observedAccessionsVec[i].accession << " " << observedAccessionsVec[i].taxID << endl;
                unmappedAccessions.push_back(observedAccessionsVec[i].accession);
                continue;
            }
            curId = it->second;
        } else {
            curId = observedAccessionsVec[i].taxID;
        }
        observedAccessionsVec[i].taxID = curId;
        observedAccessionsVec[i].speciesID = taxonomy->getTaxIdAtRank(curId, "species");
        spIdSet.insert(observedAccessionsVec[i].speciesID);
        taxIdSet.insert(curId);
        TaxID collapseTaxID = getCollapseTaxId(curId);
        taxId2speciesId[curId] = collapseTaxID;
        if (observedAccessionsVec[i].speciesID == 0) {
            cout << "Species ID is not found for accession " << observedAccessionsVec[i].accession << " " << curId << endl;
            unmappedAccessions.push_back(observedAccessionsVec[i].accession);
            exit(1);
        }
        taxId2speciesId[observedAccessionsVec[i].speciesID] = collapseTaxID;
        taxId2speciesId[collapseTaxID] = collapseTaxID;
    }

    cout << "Number of unique taxIDs: " << taxIdSet.size() << endl;
    cout << "Number of unique speciesIDs: " << spIdSet.size() << endl;
    
    string mappingFileName = dbDir + "/acc2taxid.map";
    FILE * acc2taxidFile = nullptr;
    if (isUpdating) {
        acc2taxidFile = fopen(mappingFileName.c_str(), "a");
    } else {
        acc2taxidFile = fopen(mappingFileName.c_str(), "w");
    }
    if (par.accessionLevel == 1) {
        for (auto & acc2taxid: acc2accId) {
            fprintf(acc2taxidFile, "%s\t%d\t%d\n", acc2taxid.first.c_str(), acc2taxid.second.first, acc2taxid.second.second);
        }
    } else {
        for (size_t i = 0; i < observedAccessionsVec.size(); ++i) {
            if (observedAccessionsVec[i].taxID == 0) {
                continue;
            }
            fprintf(acc2taxidFile, "%s\t%d\n", observedAccessionsVec[i].accession.c_str(), taxonomy->getOriginalTaxID(observedAccessionsVec[i].taxID));
        }        
    }   

    if (unmappedAccessions.empty()) {
        cout << "All accessions are mapped to taxonomy" << endl;
    } else {
        FILE *file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
        for (const auto & i : unmappedAccessions) {
            fprintf(file, "%s\n", i.c_str());
        }
        fclose(file);
        cout << "Unmapped accessions are written to " << dbDir + "/unmapped.txt" << endl;
    }

    fclose(acc2taxidFile);                   
}

void IndexCreator::getAccessionBatches(std::vector<Accession> & observedAccessionsVec, size_t bufferSize) {
    size_t accCnt = observedAccessionsVec.size();
    SORT_PARALLEL(observedAccessionsVec.begin(), observedAccessionsVec.end(), Accession::compare);
    vector<uint32_t> orders;
    vector<TaxID> taxIDs;
    vector<uint32_t> lengths;
    for (size_t i = 0; i < accCnt;) {
        if (observedAccessionsVec[i].speciesID == 0 || observedAccessionsVec[i].taxID == 0) {
            i++;
            continue;
        }
        TaxID currentSpeciesID = observedAccessionsVec[i].speciesID;
        uint32_t maxLength = 0, trainingFasta = 0, trainingSeq = 0;
        size_t firstBatchOfSpecies = accessionBatches.size();
        while (i < accCnt && currentSpeciesID == observedAccessionsVec[i].speciesID) {    
            uint32_t lengthSum = 0;
            size_t kmerCntSum = 0;
            orders.clear();
            lengths.clear();
            taxIDs.clear();
            uint32_t currentFasta = observedAccessionsVec[i].whichFasta;
            while (i < accCnt && currentSpeciesID == observedAccessionsVec[i].speciesID 
                   && currentFasta == observedAccessionsVec[i].whichFasta) {
                if (observedAccessionsVec[i].length > maxLength) {
                    maxLength = observedAccessionsVec[i].length;
                    trainingFasta = observedAccessionsVec[i].whichFasta;
                    trainingSeq = observedAccessionsVec[i].order;
                }
                lengthSum += observedAccessionsVec[i].length;
                totalLength += observedAccessionsVec[i].length;
                kmerCntSum += static_cast<size_t>(observedAccessionsVec[i].length * 0.4);
                orders.push_back(observedAccessionsVec[i].order);
                lengths.push_back(observedAccessionsVec[i].length);
                taxIDs.push_back(observedAccessionsVec[i].taxID);
                i++;
                if (kmerCntSum > bufferSize || lengthSum > 100'000'000 || orders.size() > 300 || (orders.size() > 100 && lengthSum > 50'000'000)) {
                    break;
                }
            }
            // Add the batch
            accessionBatches.emplace_back(currentFasta, currentSpeciesID, 0, 0, lengthSum);
            accessionBatches.back().orders = orders;
            accessionBatches.back().lengths = lengths;
            accessionBatches.back().taxIDs = taxIDs;
        }
        // Update training sequence information
        for (size_t j = firstBatchOfSpecies; j < accessionBatches.size(); ++j) {
            accessionBatches[j].trainingSeqFasta = trainingFasta;
            accessionBatches[j].trainingSeqIdx = trainingSeq;
        }
    }

    // for (size_t i = 0; i < accessionBatches.size(); ++i) {
    //     accessionBatches[i].print();
    // }

}

void IndexCreator::getSpeciesBatches() {

    std::unordered_map<TaxID, std::string> sp2assacc;
    regex regex1("(GC[AF]_[0-9]+\\.[0-9]+)");
    if (!par.repGenomeList.empty()) {
        ifstream repGenomeFile(par.repGenomeList);
        if (repGenomeFile.is_open()) {
            for (string eachLine; getline(repGenomeFile, eachLine);) {
                const size_t tabPos = eachLine.find('\t');
                if (tabPos == string::npos) {
                    cout << "Invalid line in representative genome list: " << eachLine << endl;
                    continue;
                }
                const TaxID spId = stoi(eachLine.substr(0, tabPos));
                sp2assacc[spId] = eachLine.substr(tabPos + 1);
            }
        } else {
            cout << "Cannot open file for representative genome list" << endl;
        }
    }

    std::vector<uint64_t> contigLengths;
    for (size_t i = 0; i < accessionBatches.size(); ) {
        TaxID currentSpeciesID = accessionBatches[i].speciesID;
        spBatches.emplace_back(currentSpeciesID);
        auto & currSpBatch = spBatches.back();
        uint64_t spTotalLength = 0;
        uint64_t spGenomeCnt = 0;

        // For each species
        while (i < accessionBatches.size() && currentSpeciesID == accessionBatches[i].speciesID) { 
            spGenomeCnt++;
            if (spGenomeCnt > 100) {
                unusedFastaPaths.push_back(fastaPaths[accessionBatches[i].whichFasta]);
                i++;
                continue;
            }
            uint32_t curFasta = accessionBatches[i].whichFasta;
            currSpBatch.fastaBatches.emplace_back(curFasta);
            auto & fastaBatch = currSpBatch.fastaBatches.back();
            fastaBatch.taxId = accessionBatches[i].taxIDs[0];
            
            // For each FASTA (Genome)
            contigLengths.clear();
            uint64_t fastaTotalLength = 0;
            while (i < accessionBatches.size() && currentSpeciesID == accessionBatches[i].speciesID
                   && curFasta == accessionBatches[i].whichFasta) {
                fastaBatch.accessionBatches.push_back(accessionBatches[i]);
                for (size_t j = 0; j < accessionBatches[i].lengths.size(); j++) {
                    contigLengths.push_back(accessionBatches[i].lengths[j]);
                }
                spTotalLength += accessionBatches[i].totalLength;
                fastaTotalLength += accessionBatches[i].totalLength;
                i++;
            }
            fastaBatch.length = fastaTotalLength;

            // Get contig N50 and L50
            sort(contigLengths.begin(), contigLengths.end(), greater<uint64_t>());
            uint64_t lengthSum = 0;
            for (size_t j = 0; j < contigLengths.size(); j++) {
                lengthSum += contigLengths[j];
                if (lengthSum >= fastaTotalLength / 2) {
                    fastaBatch.n50 = contigLengths[j];
                    fastaBatch.l50 = j + 1;
                    break;
                }
            }
        }
        currSpBatch.spTotalLength   = spTotalLength;
        currSpBatch.expectedKmerNum = spTotalLength * 0.4;

        // Select species representative genome
        // 1. Get midian FASTA length
        vector<uint64_t> fastaLengths;
        std::string repGenomeAssacc;
        if (sp2assacc.find(currentSpeciesID) != sp2assacc.end()){
            repGenomeAssacc = sp2assacc[currentSpeciesID];
        }
        for (size_t j = 0; j < currSpBatch.fastaBatches.size(); j++) {
            const std::string & fileName = fastaPaths[currSpBatch.fastaBatches[j].whichFasta];
            smatch assacc;
            if (regex_search(fileName, assacc, regex1)) {
                if (assacc[0] == repGenomeAssacc) {
                    currSpBatch.repGenomeFasta = currSpBatch.fastaBatches[j].whichFasta;
                    currSpBatch.repGenomeSize  = currSpBatch.fastaBatches[j].length;
                }
                
            }
            fastaLengths.push_back(currSpBatch.fastaBatches[j].length);
        }

        if (currSpBatch.repGenomeSize > 0) {
            continue;
        }

        // 2. Select representative genome using L50 and N50 among genomes with size between (0.9 * median, 1.1 * median)
        sort(fastaLengths.begin(), fastaLengths.end());
        size_t sz = fastaLengths.size();
        uint64_t medianFastaLength = (sz % 2 == 0) ? (fastaLengths[sz / 2 - 1] + fastaLengths[sz / 2]) / 2 : fastaLengths[sz / 2];
        uint64_t bestN50 = 0;
        uint64_t bestL50 = UINT64_MAX;
        for (size_t j = 0; j < currSpBatch.fastaBatches.size(); j++) {
            const auto & currFASTA = currSpBatch.fastaBatches[j];
            if (currFASTA.length >= 0.9 * medianFastaLength && currFASTA.length <= 1.1 * medianFastaLength) {
                if (currFASTA.n50 > bestN50 || (currFASTA.n50 == bestN50 && currFASTA.l50 < bestL50)) {
                    bestN50 = currFASTA.n50;
                    bestL50 = currFASTA.l50;
                    currSpBatch.repGenomeFasta = currFASTA.whichFasta;
                    currSpBatch.repGenomeSize  = currFASTA.length;
                }
            }
        }
    }

    if (unusedFastaPaths.size() > 0) {
        cout << "The following " << unusedFastaPaths.size() << " FASTA files are not used for index creation because their species have more than 100 genomes. They are written to " << dbDir + "/unused_fasta.list" << endl;
        FILE * unusedFastaFile = fopen((dbDir + "/unused_fasta.list").c_str(), "w");
        for (const auto & i : unusedFastaPaths) {
            fprintf(unusedFastaFile, "%s\n", i.c_str());
        }
        fclose(unusedFastaFile);
    }
}

// void IndexCreator::writeTargetFiles(
//     Buffer<Kmer> & kmerBuffer)
// {
//     string diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
//     string infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";
//     deltaIdxFileNames.push_back(diffIdxFileName);
//     infoFileNames.push_back(infoFileName);

//     numOfFlush++;

//     size_t bufferSize = 1024 * 1024 * 32;
//     WriteBuffer<uint16_t> diffBuffer(diffIdxFileName, bufferSize);
//     WriteBuffer<uint32_t> infoBuffer(infoFileName, bufferSize); 
//     uint64_t lastKmer = 0;

//     for (size_t i = 0; i < kmerBuffer.startIndexOfReserve; i ++) {
//         infoBuffer.write(&kmerBuffer.buffer[i].id);
//         getDiffIdx(lastKmer, kmerBuffer.buffer[i].value, diffBuffer);
//     }
//     infoBuffer.flush();
//     diffBuffer.flush();

//     kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch

// }

void IndexCreator::writeTargetFiles(
    Buffer<Kmer> & kmerBuffer, 
    const size_t * uniqKmerIdx,
    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges) 
{
    string diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    string infoFileName = dbDir + "/" + to_string(numOfFlush) + "_info";
    deltaIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);

    numOfFlush++;

    size_t bufferSize = 1024 * 1024 * 32;
    WriteBuffer<uint16_t> diffBuffer(diffIdxFileName, bufferSize);
    WriteBuffer<uint32_t> infoBuffer(infoFileName, bufferSize); 
    uint64_t lastKmer = 0;

    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
            infoBuffer.write(&kmerBuffer.buffer[uniqKmerIdx[j]].id);
            getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[j]].value, diffBuffer);
        }
    }
    infoBuffer.flush();
    diffBuffer.flush();

    kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch
}

void IndexCreator::writeTargetFiles(
    Buffer<Kmer> & kmerBuffer) 
{
    string diffIdxFileName = dbDir + "/" + to_string(numOfFlush) + "_diffIdx";
    string infoFileName    = dbDir + "/" + to_string(numOfFlush) + "_info";
    string posFileName     = dbDir + "/" + to_string(numOfFlush) + "_kmerpos";
    deltaIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);
    posFileNames.push_back(posFileName);

    numOfFlush++;

    size_t bufferSize = 1024 * 1024 * 32;
    WriteBuffer<uint16_t> diffBuffer(diffIdxFileName, bufferSize);
    WriteBuffer<uint32_t> infoBuffer(infoFileName, bufferSize); 
    WriteBuffer<uint16_t> posBuffer(posFileName, bufferSize);
    uint64_t lastKmer = 0;

    // Find the first index of garbage k-mer (UINT64_MAX)
    for (size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
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

    for (size_t i = startIdx; i < kmerBuffer.startIndexOfReserve ; i++) {
        uint16_t pos = static_cast<uint16_t>(kmerBuffer.buffer[i].tInfo.pos);
        posBuffer.write(&pos);
        infoBuffer.write(&kmerBuffer.buffer[i].id);
        getDiffIdx(lastKmer, kmerBuffer.buffer[i].value, diffBuffer);
    }
    kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch
}

void IndexCreator::writeTargetFilesAndSplits(
    Buffer<Kmer> & kmerBuffer,
    const size_t * uniqKmerIdx,
    size_t & uniqKmerCnt,
    const vector<pair<size_t, size_t>> & uniqKmerIdxRanges)
{    
    // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t sizeOfSplit = uniqKmerCnt / (par.splitNum - 1);
    size_t offsetList[par.splitNum + 1];
    int offsetListIdx = 1;
    for(int os = 0; os < par.splitNum; os++){
        offsetList[os] = os * sizeOfSplit;
    }
    offsetList[par.splitNum] = UINT64_MAX;
    DiffIdxSplit * splitList = new DiffIdxSplit[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    int splitListIdx = 1;
    int splitCheck = 0;

    numOfFlush++;
    size_t bufferSize = 1024 * 1024 * 32;
    uint64_t lastKmer = 0;
    WriteBuffer<uint16_t> diffBuffer(dbDir + "/diffIdx", bufferSize);
    
    WriteBuffer<uint32_t> infoBuffer(dbDir + "/info", bufferSize); 
    for (size_t i = 0; i < uniqKmerIdxRanges.size(); i ++) {
        for (size_t j = uniqKmerIdxRanges[i].first; j < uniqKmerIdxRanges[i].second; j ++) {
            infoBuffer.write(&kmerBuffer.buffer[uniqKmerIdx[j]].id);
            getDiffIdx(lastKmer, kmerBuffer.buffer[uniqKmerIdx[j]].value, diffBuffer);
            // Write split info
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
    }
    cout << "Written k-mer count : " << infoBuffer.writeCnt << endl;

    FILE * deltaIdxSplitFile = fopen((dbDir + "/split").c_str(), "wb");
    if (deltaIdxSplitFile == nullptr) {
        cout << "Cannot open the file for writing target DB" << endl;
        return;
    }
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, deltaIdxSplitFile);
    delete[] splitList;
    fclose(deltaIdxSplitFile);
    
    kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch
}

void IndexCreator::writeTargetFilesAndSplits(
    Buffer<Kmer> & kmerBuffer)
{    

    DiffIdxSplit * splitList = new DiffIdxSplit[par.splitNum];
    memset(splitList, 0, sizeof(DiffIdxSplit) * par.splitNum);
    int splitListIdx = 1;
    int splitCheck = 0;

    numOfFlush++;
    size_t bufferSize = 1024 * 1024 * 32;
    uint64_t lastKmer = 0;
    WriteBuffer<uint16_t> diffBuffer(dbDir + "/diffIdx", bufferSize);
    WriteBuffer<uint16_t> posBuffer(dbDir + "/kmerpos", bufferSize);
    WriteBuffer<uint32_t> infoBuffer(dbDir + "/info", bufferSize); 

    // Find the first index of garbage k-mer (UINT64_MAX)
    for (size_t checkN = kmerBuffer.startIndexOfReserve - 1; checkN != 0; checkN--){
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



        // To make differential index splits
    uint64_t AAofTempSplitOffset = UINT64_MAX;
    size_t splitSize = (kmerBuffer.startIndexOfReserve - startIdx) / (par.splitNum - 1);
    vector<size_t> splitOffsets;
    for (size_t i = 1; i < par.splitNum; i++) {
        splitOffsets.push_back(i * splitSize);
    }
    splitOffsets.push_back(UINT64_MAX);
    int offsetListIdx = 1;
    for (size_t i = startIdx; i < kmerBuffer.startIndexOfReserve ; i++) {
        uint16_t pos = static_cast<uint16_t>(kmerBuffer.buffer[i].tInfo.pos);
        posBuffer.write(&pos);
        infoBuffer.write(&kmerBuffer.buffer[i].id);
        getDiffIdx(lastKmer, kmerBuffer.buffer[i].value, diffBuffer);

        // Write split info
        if (AminoAcidPart(lastKmer) != AAofTempSplitOffset && splitCheck == 1) {
            splitList[splitListIdx++] = {lastKmer, diffBuffer.writeCnt, infoBuffer.writeCnt};
            splitCheck = 0;
        }
        if (infoBuffer.writeCnt == splitOffsets[offsetListIdx]) {
            AAofTempSplitOffset = AminoAcidPart(lastKmer);
            splitCheck = 1;
            offsetListIdx++;
        }
    }
    cout << "Written k-mer count : " << infoBuffer.writeCnt << endl;
    FILE * deltaIdxSplitFile = fopen((dbDir + "/split").c_str(), "wb");
    if (deltaIdxSplitFile == nullptr) {
        cout << "Cannot open the file for writing target DB" << endl;
        return;
    }
    fwrite(splitList, sizeof(DiffIdxSplit), par.splitNum, deltaIdxSplitFile);
    delete[] splitList;
    fclose(deltaIdxSplitFile);
    
    kmerBuffer.startIndexOfReserve = 0; // Reset the buffer for the next batch
}



void IndexCreator::getDiffIdx(
    uint64_t & lastKmer,
    uint64_t entryToWrite,
    WriteBuffer<uint16_t> & diffBuffer) 
{
    uint64_t kmerdiff = entryToWrite - lastKmer;
    uint16_t buffer[5];
    int idx = 3;
    buffer[4] = SET_END_FLAG(GET_15_BITS(kmerdiff));
    kmerdiff >>= 15U;
    while (kmerdiff) {
        uint16_t toWrite = GET_15_BITS(kmerdiff);
        kmerdiff >>= 15U;
        buffer[idx] = toWrite;
        idx--;
    }
    diffBuffer.write((buffer + idx + 1), (4 - idx));
    lastKmer = entryToWrite;
}


void IndexCreator::load_assacc2taxid(const string & mappingFile, unordered_map<string, int> & assacc2taxid){
    string key, value;
    ifstream map;
    map.open(mappingFile);
    if(map.is_open()){
        while(getline(map,key,'\t')){
            getline(map, value, '\n');
            assacc2taxid[key] = stoi(value);
            if(key[2] == 'F'){
                key[2] = 'A';
                assacc2taxid[key] = stoi(value);
            }
        }
    } else{
        cout<<"Cannot open file for mappig from assemlby accession to tax ID"<<endl;
    }
    map.close();
}


bool IndexCreator::extractKmerFromSixFrames(
    Buffer<Kmer> & kmerBuffer,
    std::vector<std::atomic<bool>> & batchChecker,
    size_t &processedBatchCnt
) {
    std::atomic<int> hasOverflow{0};
    #pragma omp parallel default(none), shared(kmerBuffer, batchChecker, processedBatchCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        size_t estimatedKmerCnt = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t batchIdx = 0; batchIdx < accessionBatches.size(); batchIdx ++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (batchChecker[batchIdx].exchange(true, std::memory_order_acq_rel))
                continue; 
            
            // Estimate the number of k-mers to be extracted from current split
            size_t totalLength = 0;
            for (size_t p = 0; p < accessionBatches[batchIdx].lengths.size(); p++) {
                totalLength += accessionBatches[batchIdx].lengths[p];
            }

            if (par.syncmer) {
                estimatedKmerCnt = static_cast<size_t>(
                    (totalLength * 2.5) / ((12 - par.smerLen + 1) / 2.0)
                );
            } else {
                estimatedKmerCnt = static_cast<size_t>(
                    totalLength * 2.5
                );
            }
                
            // Process current split if buffer has enough space.
            size_t posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
            if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                KSeqWrapper* kseq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].whichFasta].c_str());
                size_t seqCnt = 0;
                size_t idx = 0;
                while (kseq->ReadEntry()) {
                    if (seqCnt == accessionBatches[batchIdx].orders[idx]) {
                        if (accessionBatches[batchIdx].taxIDs[idx] == 0) {
                            #pragma omp critical
                            {
                                accessionBatches[batchIdx].print();
                                exit(1);
                            }
                        }
                        const KSeqWrapper::KSeqEntry & e = kseq->entry;

                        // Mask low complexity regions
                        char *maskedSeq = nullptr;
                        if (par.maskMode) {
                            maskedSeq = new char[e.sequence.l + 1]; 
                            SeqIterator::maskLowComplexityRegions((unsigned char *) e.sequence.s, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                            maskedSeq[e.sequence.l] = '\0';
                        } else {
                            maskedSeq = e.sequence.s;
                        }

                        kmerExtractor->extractKmer_dna2aa(
                            maskedSeq,
                            e.sequence.l,
                            kmerBuffer,
                            posToWrite,
                            accessionBatches[batchIdx].taxIDs[idx],
                            accessionBatches[batchIdx].speciesID);
                            
                        idx++;
                        if (par.maskMode) {
                            delete[] maskedSeq;
                        }
                        if (idx == accessionBatches[batchIdx].lengths.size()) {
                            break;
                        }
                    }
                    seqCnt++;
                }
                delete kseq;
                __sync_fetch_and_add(&processedBatchCnt, 1);
                #pragma omp critical
                {
                    cout << processedBatchCnt << " batches processed out of " << accessionBatches.size() << endl;
                        // cout << fastaPaths[accessionBatches[batchIdx].whichFasta] << " processed\n";
                }
            } else {
                batchChecker[batchIdx].store(false, std::memory_order_release);
                hasOverflow.fetch_add(1, std::memory_order_relaxed);
                __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
            }
        }
    }

    // cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
    
}


















void IndexCreator::makeIntergenicKmerList(
    const std::string & fastaFileName,
    vector<uint64_t> & intergenicKmers,
    ProdigalWrapper * prodigal)
{
    KSeqWrapper* kseq = KSeqFactory(fastaFileName.c_str());
    while (kseq->ReadEntry()) {
        const KSeqWrapper::KSeqEntry & e = kseq->entry;
        prodigal->getPredictedGenes((unsigned char *)(e.sequence.s), e.sequence.l);
        generateIntergenicKmerList(
            prodigal->genes, 
            prodigal->nodes,
            prodigal->getNumberOfPredictedGenes(),
            intergenicKmers,
            e.sequence.s);
    }
    delete kseq;
}

size_t IndexCreator::fillTargetKmerBuffer2(
    Buffer<Kmer> &kmerBuffer,
    std::vector<std::atomic<bool>> & batchChecker,
    size_t &processedSpCnt,
    const LocalParameters &par) 
{
    std::atomic<int> hasOverflow{0};
#pragma omp parallel default(none), shared(kmerBuffer, batchChecker, processedSpCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        size_t posToWrite;
        size_t orfNum;
        vector<SequenceBlock> fragments;
        vector<uint64_t> intergenicKmers;
        size_t maxSeqLen = 1000;
        char *maskedSeq = nullptr;
        if (par.maskMode) {
            maskedSeq = new char[maxSeqLen + 1];
        }
        int tempCheck = 1;

#pragma omp for schedule(dynamic, 1)
        for (size_t spIdx = 0; spIdx < spBatches.size(); spIdx ++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (batchChecker[spIdx].exchange(true, std::memory_order_acq_rel))
                continue; 
            
            intergenicKmers.clear();

            // Estimate the number of k-mers to be extracted from current split
            size_t totalLength = spBatches[spIdx].spTotalLength;
            size_t estimatedKmerCnt = static_cast<size_t>((totalLength * 1.3) / 3.0);
            if (par.syncmer) {
                estimatedKmerCnt = static_cast<size_t>(
                    (totalLength * 1.3 / 3.0) / ((metamerPattern->windowSize - par.smerLen + 1) / 2.0)
                );
            } 

            if (estimatedKmerCnt > kmerBuffer.bufferSize) {
                cout << "Estimated k-mer count for species " << spBatches[spIdx].speciesID << " exceeds buffer size. Stop processing." << endl;
                exit(1);
            }
                
            
            // Process current species if buffer has enough space.
            posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
            size_t startPosToWrite = posToWrite;
            if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                // Train Prodigal
                ProdigalWrapper * prodigal = new ProdigalWrapper();
                if (spBatches[spIdx].repGenomeSize < 100'000 || 
                    ((taxonomy->getEukaryotaTaxID() != 0) && 
                     (taxonomy->IsAncestor(spBatches[spIdx].speciesID, taxonomy->getEukaryotaTaxID()))
                    )) {
                    prodigal->is_meta = 1;
                    prodigal->trainMeta(fastaPaths[spBatches[spIdx].repGenomeFasta]);
                } else {
                    prodigal->trainASpecies(fastaPaths[spBatches[spIdx].repGenomeFasta]);
                }

                // Make intergenic k-mer list to guide ORF extension
                makeIntergenicKmerList(
                    fastaPaths[spBatches[spIdx].repGenomeFasta],
                    intergenicKmers,
                    prodigal);
                
                // Extract k-mers from each fasta batch
                const auto & fastaBatches = spBatches[spIdx].fastaBatches;
                for (size_t i = 0; i < fastaBatches.size(); i++) {
                    KSeqWrapper* kseq = KSeqFactory(fastaPaths[fastaBatches[i].whichFasta].c_str());
                    const auto & accessions = fastaBatches[i].accessionBatches;
                    TaxID taxId = fastaBatches[i].taxId;
                    
                    // Process each sequence
                    uint64_t genomicPos  = 0;
                    uint64_t scaleFactor = 0;
                    if (fastaBatches[i].whichFasta == spBatches[spIdx].repGenomeFasta) {
                        if (spBatches[spIdx].repGenomeSize >= 65535) {
                            scaleFactor = (65535ULL << 32) / spBatches[spIdx].repGenomeSize; 
                        } else {
                            scaleFactor = UINT64_MAX;
                        }
                    }

                    while (kseq->ReadEntry()) {
                        const KSeqWrapper::KSeqEntry & e = kseq->entry;

                        if (par.maskMode) {
                            if (e.sequence.l > maxSeqLen) {
                                delete[] maskedSeq;
                                maxSeqLen = e.sequence.l;
                                maskedSeq = new char[maxSeqLen + 1];
                            }
                            SeqIterator::maskLowComplexityRegions((unsigned char *) e.sequence.s, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                            maskedSeq[e.sequence.l] = '\0';
                        } else {
                            maskedSeq = e.sequence.s;
                        }

                        fragments.clear();
                        if (cdsInfoMap.find(string(e.name.s)) != cdsInfoMap.end()) {
                            // USE PROVIDED CDS ANNOTATION
                            devideToCdsAndNonCds(e.sequence.l, cdsInfoMap[string(e.name.s)], fragments);

                            for (size_t f = 0; f < fragments.size(); f++) {
                                tempCheck = kmerExtractor->extractTargetKmers(
                                                maskedSeq,
                                                kmerBuffer,
                                                posToWrite,
                                                genomicPos,
                                                taxId,
                                                fragments[f],
                                                scaleFactor);
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                }
                            }
                        } else {
                            // PREDICT GENES USING PRODIGAL
                            orfNum = 0;                            
                            prodigal->getPredictedGenes((unsigned char *) e.sequence.s, e.sequence.l);
                            prodigal->removeCompletelyOverlappingGenes();
                            prodigal->getExtendedORFs(prodigal->finalGenes, prodigal->nodes, fragments,
                                                      prodigal->fng, e.sequence.l, orfNum, intergenicKmers, e.sequence.s);
                            // for (size_t f = 0; f < fragments.size(); f++) {
                            //     fragments[f].printSequenceBlock();
                            // }
                            for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                tempCheck = kmerExtractor->extractTargetKmers(
                                                maskedSeq,
                                                kmerBuffer,
                                                posToWrite,
                                                genomicPos,
                                                taxId,
                                                fragments[orfCnt],
                                                scaleFactor);
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                }
                            }
                        }
                        genomicPos += e.sequence.l;
                    } 
                    delete kseq; // End of processing each genome
                }                

                // Sort the k-mers extracted from the current species
                std::sort(kmerBuffer.buffer + startPosToWrite, kmerBuffer.buffer + posToWrite, Kmer::compareTargetKmerPerSpecies);

                size_t i = startPosToWrite;
                uint32_t distinctAAgroupCnt = 0;
                uint32_t distinctDNAgroupCnt = 0;
                while (i < posToWrite) {
                    size_t currAA = AminoAcidPart(kmerBuffer.buffer[i].value);
                    distinctAAgroupCnt++;

                    // ==========================================
                    // PHASE 1: PRE-SCAN THE AA GROUP FOR A REP. K-MER POSITION
                    // ==========================================
                    uint32_t sharedPos = kmerBuffer.buffer[i].tInfo.pos;
                    size_t j = i;
                    while ((sharedPos == 0) 
                            && j < posToWrite 
                            && AminoAcidPart(kmerBuffer.buffer[j].value) == currAA) {
                        sharedPos = std::max(sharedPos, kmerBuffer.buffer[j].tInfo.pos);
                        j++;
                    }

                    // ==========================================
                    // PHASE 2: DEDUPLICATE AND APPLY POSITION
                    // ==========================================
                    // Process all DNA variants within this AA group (from 'i' up to 'endOfAA')
                    while (i < posToWrite && AminoAcidPart(kmerBuffer.buffer[i].value) == currAA) {                        
                        Kmer * firstOfCurrDNA = &kmerBuffer.buffer[i];

                        // All synonymous DNA k-mers inherit the shared position (Physical or Virtual)
                        firstOfCurrDNA->tInfo.pos = sharedPos;  
                        bool duplicated = false;
                        bool isFirst = true;
                        distinctDNAgroupCnt++;

                        // Deduplicate exact identical DNA k-mers
                        while (i < posToWrite && kmerBuffer.buffer[i].value == firstOfCurrDNA->value) {
                            if (kmerBuffer.buffer[i].tInfo.taxId != firstOfCurrDNA->tInfo.taxId) {
                                duplicated = true;
                            }
                            if (isFirst) {
                                isFirst = false;
                            } else {
                                kmerBuffer.buffer[i] = {0, 0, 0}; // Tombstone duplicates
                            }
                            i++;
                        }
                    
                        // Promote shared DNA k-mers to the Species level LCA
                        if (duplicated) {
                            firstOfCurrDNA->tInfo.taxId = spBatches[spIdx].speciesID;
                        }
                    }
                }


                // // ==========================================
                // // PHASE 3: CONVERT TO 65536 FIXED BINS
                // // ==========================================
                // // The total pan-genome size is now known. 
                // uint64_t totalPanGenomeSize = virtualPos;
                // if (totalPanGenomeSize > 0) {
                //     uint64_t scaleFactor = (65535ULL << 32) / totalPanGenomeSize;
                //     for (size_t j = startPosToWrite; j < posToWrite; ++j) {
                //         Kmer * currKmer = &kmerBuffer.buffer[j];
                //         if (currKmer->tInfo.taxId != 0) { 
                //             uint64_t binID = (static_cast<uint64_t>(currKmer->tInfo.pos) * scaleFactor) >> 32;
                //             if (binID > 65535) binID = 65535;
                //             currKmer->tInfo.pos = static_cast<uint32_t>(binID);
                //         }
                //     }
                // }

                // cout << "totalPanGenomeSize: " << totalPanGenomeSize << endl;
                // cout << "distinctAAgroupCnt: " << distinctAAgroupCnt << endl;
                // cout << "distinctDNAgroupCnt: " << distinctDNAgroupCnt << endl;

                __sync_fetch_and_add(&processedSpCnt, 1);
                #pragma omp critical
                {
                    cout << processedSpCnt << " batches processed out of " << spBatches.size() << endl;
                }
                delete prodigal;
                // --- End of processing current species   
            } else {
                batchChecker[spIdx].store(false, std::memory_order_release);
                hasOverflow.fetch_add(1, std::memory_order_relaxed);
                __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
            }
        }
        if (par.maskMode) {
            delete[] maskedSeq;
        }
    } // End of parallel region

    return 0;
}












void IndexCreator::devideToCdsAndNonCds(
    size_t seqLen,
    const vector<CDSinfo> &cdsInfo, 
    std::vector<SequenceBlock> & fragments) 
{
    for (size_t i = 0; i < cdsInfo.size(); i++) {
        size_t locNum = cdsInfo[i].loc.size();
        int currStartCodonPos = 0;
        for (size_t j = 0; j < locNum; j++) {
            // Extend 21 bases to both sides for k-mer from CDS boudaries
            size_t begin = cdsInfo[i].loc[j].first - 1; // convert to 0-based (inclusive)
            size_t end = cdsInfo[i].loc[j].second - 1;  // convert to 0-based (inclusive)
            if (j == 0) {
                int k = 0;
                while (k < (this->kmerLen - 1) && begin >= 3) {
                    begin -= 3;
                    currStartCodonPos += 3;
                    k++;
                }
            }
            if (j == locNum - 1) {
                int k = 0;
                while (k < (this->kmerLen - 1) && end + 3 < seqLen) {
                    end += 3;
                    k++;
                }
            }    
            fragments.emplace_back(begin, end, cdsInfo[i].isComplement ? -1 : 1);
        }
    }

    // Get non-CDS that are not in the CDS list
    bool * checker = new bool[seqLen];
    memset(checker, true, seqLen * sizeof(bool));
    for (size_t i = 0; i < cdsInfo.size(); i++) {
        for (size_t j = 0; j < cdsInfo[i].loc.size(); j++) {
            for (size_t k = cdsInfo[i].loc[j].first - 1; k < cdsInfo[i].loc[j].second; k++) {
                checker[k] = false;
            }
        }
    }

    size_t len;
    size_t i = 0;
    while (i < seqLen) {
        len = 0;
        while (i < seqLen && checker[i]) {
            i++;
            len++;
        }
        if (len > 32) {
            fragments.emplace_back(i - len, i - 1, 1);
        }
        i ++;
    }
}


size_t IndexCreator::fillTargetKmerBuffer(Buffer<Kmer> &kmerBuffer,
                                          std::vector<std::atomic<bool>> & batchChecker,
                                          size_t &processedBatchCnt,
                                          const LocalParameters &par) {
    std::atomic<int> hasOverflow{0};
#pragma omp parallel default(none), shared(kmerBuffer, batchChecker, processedBatchCnt, hasOverflow, par, cout)
    {
        ProbabilityMatrix probMatrix(*subMat);
        SeqIterator seqIterator(par);
        size_t posToWrite;
        size_t orfNum;
        vector<SequenceBlock> extendedORFs;
        std::vector<uint64_t> standardList;
        std::vector<uint64_t> currentList;
        size_t lengthOfTrainingSeq = 0;
        char *rcomp;
        vector<uint64_t> intergenicKmers;
        vector<string> cds;
        vector<string> nonCds;
        bool trained = false;
        size_t estimatedKmerCnt = 0;
#pragma omp for schedule(dynamic, 1)
        for (size_t batchIdx = 0; batchIdx < accessionBatches.size(); batchIdx ++) {
            if (hasOverflow.load(std::memory_order_acquire))
                continue;
            
            if (batchChecker[batchIdx].exchange(true, std::memory_order_acq_rel))
                continue; 
            
            intergenicKmers.clear();
            standardList.clear();
            currentList.clear();

            // Estimate the number of k-mers to be extracted from current split
            size_t totalLength = 0;
            for (size_t p = 0; p < accessionBatches[batchIdx].lengths.size(); p++) {
                totalLength += accessionBatches[batchIdx].lengths[p];
            }

            if (par.syncmer) {
                estimatedKmerCnt = static_cast<size_t>(
                    (totalLength * 1.3 / 3.0) / ((metamerPattern->kmerLen - par.smerLen + 1) / 2.0)
                );
            } else {
                estimatedKmerCnt = static_cast<size_t>(
                    (totalLength * 1.3) / 3.0
                );
            }
                
            ProdigalWrapper * prodigal = new ProdigalWrapper();
            trained = false;

            // Process current split if buffer has enough space.
            posToWrite = kmerBuffer.reserveMemory(estimatedKmerCnt);
            if (posToWrite + estimatedKmerCnt < kmerBuffer.bufferSize) {
                KSeqWrapper* kseq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].whichFasta].c_str());
                size_t seqCnt = 0;
                size_t idx = 0;
                while (kseq->ReadEntry()) {
                    if (seqCnt == accessionBatches[batchIdx].orders[idx]) {
                        bool doMasking = par.maskMode;
                        bool isNotToMask = taxonomy->isAunderB(accessionBatches[batchIdx].taxIDs[idx], taxaNotToMask);
                        if (isNotToMask) {
                            #pragma omp critical
                            {
                                cout << "Masking is not applied to " << kseq->entry.name.s << " because it belongs to a taxon that is set to be not masked.\n";
                            }
                        }
                        doMasking = doMasking && !isNotToMask;

                        if (accessionBatches[batchIdx].taxIDs[idx] == 0) {
                            #pragma omp critical
                            {
                            accessionBatches[batchIdx].print();
                            exit(1);
                            }
                        }
                        TaxID collapseTaxID = getCollapseTaxId(accessionBatches[batchIdx].taxIDs[idx]);
                        const KSeqWrapper::KSeqEntry & e = kseq->entry;
                        // Mask low complexity regions
                        char *maskedSeq = nullptr;
                        if (doMasking) {
                            maskedSeq = new char[e.sequence.l + 1]; // TODO: reuse the buffer
                            SeqIterator::maskLowComplexityRegions((unsigned char *) e.sequence.s, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                            maskedSeq[e.sequence.l] = '\0';
                        } else {
                            maskedSeq = e.sequence.s;
                        }

                        orfNum = 0;
                        extendedORFs.clear();
                        int tempCheck = 0;
                        if (par.readingFrame != 0) {
                            // USE PROVIDED READING FRAME
                
                            tempCheck = kmerExtractor->extractTargetKmers(
                                            maskedSeq,
                                            kmerBuffer,
                                            posToWrite,
                                            accessionBatches[batchIdx].taxIDs[idx],
                                            collapseTaxID,
                                            {par.readingFrame-1, (int) e.sequence.l - 1, par.readingFrame < 3 ? 1 : -1});
                            if (tempCheck == -1) {
                                cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                            }   
                        }
                        else if (cdsInfoMap.find(string(e.name.s)) != cdsInfoMap.end()) {
                            // USE PROVIDED CDS ANNOTATION

                            cds.clear();
                            nonCds.clear();
                            seqIterator.devideToCdsAndNonCds(maskedSeq,
                                                             e.sequence.l,
                                                             cdsInfoMap[string(e.name.s)],
                                                             cds,
                                                             nonCds);

                            for (size_t cdsCnt = 0; cdsCnt < cds.size(); cdsCnt ++) {
                                tempCheck = kmerExtractor->extractTargetKmers(
                                                cds[cdsCnt].c_str(),
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                collapseTaxID,
                                                {0, (int) cds[cdsCnt].length() - 1, 1});
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                }
                            }
                            for (size_t nonCdsCnt = 0; nonCdsCnt < nonCds.size(); nonCdsCnt ++) {
                                tempCheck = kmerExtractor->extractTargetKmers(
                                                nonCds[nonCdsCnt].c_str(),
                                                kmerBuffer,
                                                posToWrite,
                                                accessionBatches[batchIdx].taxIDs[idx],
                                                collapseTaxID,
                                                {0, (int) cds[nonCdsCnt].length() - 1, 1});
                                if (tempCheck == -1) {
                                    cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                }
                            }
                        } else {
                            // USE PRODIGAL
                            if (!trained) {
                                KSeqWrapper* training_seq = KSeqFactory(fastaPaths[accessionBatches[batchIdx].trainingSeqFasta].c_str()); 
                                size_t seqCnt = 0;
                                while (training_seq->ReadEntry()) {
                                    if (seqCnt == accessionBatches[batchIdx].trainingSeqIdx) {
                                        break;
                                    }
                                    seqCnt++;
                                }
                                lengthOfTrainingSeq = training_seq->entry.sequence.l;
                                prodigal->is_meta = 0;
                                if (lengthOfTrainingSeq < 100'000 || 
                                    ((taxonomy->getEukaryotaTaxID() != 0) && 
                                     (taxonomy->IsAncestor(accessionBatches[batchIdx].speciesID, taxonomy->getEukaryotaTaxID()))
                                    )) {
                                    prodigal->is_meta = 1;
                                    prodigal->trainMeta((unsigned char *) training_seq->entry.sequence.s, 
                                                        training_seq->entry.sequence.l);
                                } else {
                                    prodigal->trainASpecies((unsigned char *) training_seq->entry.sequence.s,
                                                            training_seq->entry.sequence.l);
                                }

                                // Generate intergenic 23-mer list. It is used to determine extension direction of intergenic sequences.
                                prodigal->getPredictedGenes((unsigned char *) training_seq->entry.sequence.s,
                                                            training_seq->entry.sequence.l);
                                generateIntergenicKmerList(prodigal->genes, prodigal->nodes,
                                                                       prodigal->getNumberOfPredictedGenes(),
                                                                       intergenicKmers, training_seq->entry.sequence.s);
                                // Get min k-mer hash list for determining strandness
                                standardList = getMinHashList(training_seq->entry.sequence.s);
                                delete training_seq;
                                trained = true;
                            }
                            currentList.clear();
                            currentList = getMinHashList(e.sequence.s);
                            if (compareMinHashList(standardList, currentList, lengthOfTrainingSeq, e.sequence.l)) {
                                // Get extended ORFs
                                prodigal->getPredictedGenes((unsigned char *) e.sequence.s, e.sequence.l);
                                prodigal->removeCompletelyOverlappingGenes();
                                prodigal->getExtendedORFs_fixed(
                                    prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                    prodigal->fng, e.sequence.l,
                                    orfNum, intergenicKmers, e.sequence.s, metamerPattern->windowSize * 3);
                                // Get k-mers from extended ORFs
                                for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                    tempCheck = kmerExtractor->extractTargetKmers(
                                                    maskedSeq,
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    collapseTaxID,
                                                    extendedORFs[orfCnt]);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }
                            } else { // Reverse complement
                                rcomp = reverseComplement(e.sequence.s, e.sequence.l);
                                // Get extended ORFs
                                prodigal->getPredictedGenes((unsigned char *) rcomp, e.sequence.l);
                                prodigal->removeCompletelyOverlappingGenes();
                                prodigal->getExtendedORFs_fixed(
                                    prodigal->finalGenes, prodigal->nodes, extendedORFs,
                                    prodigal->fng, e.sequence.l,
                                    orfNum, intergenicKmers, rcomp, metamerPattern->windowSize * 3);

                                // Get reverse masked sequence
                                if (doMasking) {
                                    delete[] maskedSeq;
                                    maskedSeq = new char[e.sequence.l + 1];
                                    SeqIterator::maskLowComplexityRegions((unsigned char *) rcomp, (unsigned char *) maskedSeq, probMatrix, par.maskProb, subMat);
                                    maskedSeq[e.sequence.l] = '\0';
                                } else {
                                    maskedSeq = rcomp;
                                }

                                for (size_t orfCnt = 0; orfCnt < orfNum; orfCnt++) {
                                    tempCheck = kmerExtractor->extractTargetKmers(
                                                    maskedSeq,
                                                    kmerBuffer,
                                                    posToWrite,
                                                    accessionBatches[batchIdx].taxIDs[idx],
                                                    collapseTaxID,
                                                    extendedORFs[orfCnt]);
                                    if (tempCheck == -1) {
                                        cout << "ERROR: Buffer overflow " << e.name.s << e.sequence.l << endl;
                                    }
                                }
                                free(rcomp);  
                            }                            
                        }
                        idx++;
                        if (doMasking) {
                            delete[] maskedSeq;
                        }
                        if (idx == accessionBatches[batchIdx].lengths.size()) {
                            break;
                        }
                    }
                    seqCnt++;
                }
                delete kseq;
                __sync_fetch_and_add(&processedBatchCnt, 1);
                // #pragma omp critical
                // {
                //     cout << processedBatchCnt << " batches processed out of " << accessionBatches.size() << endl;
                //         // cout << fastaPaths[accessionBatches[batchIdx].whichFasta] << " processed\n";
                // }
            } else {
                batchChecker[batchIdx].store(false, std::memory_order_release);
                hasOverflow.fetch_add(1, std::memory_order_relaxed);
                __sync_fetch_and_sub(&kmerBuffer.startIndexOfReserve, estimatedKmerCnt);
            }
            delete prodigal;   
        }
    }

    // cout << "Before return: " << kmerBuffer.startIndexOfReserve << endl;
    return 0;
}


void IndexCreator::writeDbParameters() {
    FILE *handle = fopen(paramterFileName.c_str(), "w");
    if (handle == NULL) {
        Debug(Debug::ERROR) << "Could not open " << paramterFileName << " for writing\n";
        EXIT(EXIT_FAILURE);
    }
    fprintf(handle, "DB_name\t%s\n", par.dbName.c_str());
    fprintf(handle, "Creation_date\t%s\n", par.dbDate.c_str());
    fprintf(handle, "Metabuli commit used to create the DB\t%s\n", version);
    fprintf(handle, "Spaced_kmer_mask\t%s\n", par.spaceMask.c_str());
    fprintf(handle, "Accession_level\t%d\n", par.accessionLevel);
    fprintf(handle, "Mask_mode\t%d\n", par.maskMode);
    fprintf(handle, "Mask_prob\t%f\n", par.maskProb);
    fprintf(handle, "Skip_redundancy\t1\n");
    fprintf(handle, "Syncmer\t%d\n", par.syncmer);
    if (par.syncmer == 1) {
        fprintf(handle, "Syncmer_len\t%d\n", par.smerLen);
    }
    fprintf(handle, "Kmer_format\t%d\n", kmerFormat);
    fprintf(handle, "Total_seq_length\t%lu\n", totalLength);
    fprintf(handle, "Kmer_position\t%d\n", par.storeKmerPos);
    fprintf(handle, "Collapse_rank\t%s\n", par.collapseRank.c_str());

    if (!par.customMetamer.empty()) {
        // Read the custom metamer file and write to parameter file
        ifstream metamerFile(par.customMetamer);
        if (!metamerFile.is_open()) {
            Debug(Debug::ERROR) << "Could not open " << par.customMetamer << " for reading\n";
            EXIT(EXIT_FAILURE);
        }
        string line;
        bool inSection = false;
        while (getline(metamerFile, line)) {
            if (line.empty()) {
                continue;
            }

            if (line.find("===BEGIN_CUSTOM_METAMER===") != string::npos) {
                inSection = true;
            }

            if (inSection) {            
               fprintf(handle, "%s\n", line.c_str());
            }

            // Stop printing AFTER we have printed the END tag
            if (line.find("===END_CUSTOM_METAMER===") != string::npos) {
                inSection = false;
            }
        }
        metamerFile.close();
    }
    fclose(handle);
}


void IndexCreator::loadCdsInfo(const string & cdsInfoFileList) {
    uint32_t prtId = 1;
    ifstream cdsInfoList(cdsInfoFileList);
    if (!cdsInfoList.is_open()) {
        Debug(Debug::ERROR) << "Could not open " << cdsInfoFileList << " for reading\n";
        EXIT(EXIT_FAILURE);
    }

    string cdsInfoFile;
    while (getline(cdsInfoList, cdsInfoFile)) {
        if (!FileUtil::fileExists(cdsInfoFile.c_str())) {
            Debug(Debug::ERROR) << "Could not open " << cdsInfoFile << " for reading\n";
            EXIT(EXIT_FAILURE);
        }
        KSeqWrapper* kseq = KSeqFactory(cdsInfoFile.c_str());
        while (kseq->ReadEntry()) {
            const KSeqWrapper::KSeqEntry & e = kseq->entry;
            string ename = e.name.s;
            size_t start = ename.find('|') + 1;
            size_t end = ename.find('.', start);
            string accession = ename.substr(start, end - start + 2);
            int frame = 1;
            string comment = e.comment.s;
            while (true) {
                start = comment.find('[', end) + 1;
                end = comment.find(']', start);
                if (start == string::npos) { break;}
                size_t equalPos = comment.find('=', start);
                string feature = comment.substr(start, equalPos - start);
                string value = comment.substr(equalPos + 1, end - equalPos - 1);
                if (feature == "pseudo") {
                    break;
                } else if (feature == "protein" && value == "hypothetical protein") {
                    break;
                } else if (feature == "frame") {
                    frame = stoi(value);
                } else if (feature == "protein_id") {
                    cdsInfoMap[accession].emplace_back(CDSinfo(prtId++, frame));
                } else if (feature == "location") {
                    size_t complementPos = value.find('c');
                    bool isComplement = (complementPos != string::npos);
                    if (isComplement) {
                        cdsInfoMap[accession].back().isComplement = true;
                        value = value.substr(complementPos + 11, value.size() - complementPos - 12);
                    } else {
                        cdsInfoMap[accession].back().isComplement = false;
                    }
                    size_t joinPos = value.find('j');
                    if (joinPos != string::npos) {
                        value = value.substr(joinPos + 5, value.size() - joinPos - 6);
                    }                
                    size_t commaPos = value.find(',');
                    size_t dotPos;
                    string locationBegin, locationEnd;
                    while (commaPos != string::npos) {
                        dotPos = value.find('.');
                        if (dotPos > commaPos) {
                            locationBegin = value.substr(0, commaPos);
                            locationEnd = value.substr(0, commaPos);
                        } else {
                            locationBegin = value.substr(0, dotPos);
                            locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                        }
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
                    if (dotPos == string::npos) {
                        locationBegin = value;
                        locationEnd = value;
                    } else {
                        locationBegin = value.substr(0, dotPos);
                        locationEnd = value.substr(dotPos + 2, commaPos - dotPos - 2);
                    }
                    // cout << value << endl;
                    // cout << locationBegin << endl;
                    // cout << locationEnd << endl;
                    if (locationBegin[0] == '<') {
                        locationBegin = locationBegin.substr(1, locationBegin.size() - 1);
                    }
                    if (locationEnd[0] == '>') {
                        locationEnd = locationEnd.substr(1, locationEnd.size() - 1);
                    }
                    cdsInfoMap[accession].back().loc.emplace_back(stoi(locationBegin), stoi(locationEnd));

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
        delete kseq;
    }
    for (auto & cdsInfo : cdsInfoMap) {
        cout << "CDS info for " << cdsInfo.first << ": ";   
    }
    cdsInfoList.close();
}

void IndexCreator::loadMergedTaxIds(const std::string &mergedFile, unordered_map<TaxID, TaxID> & old2new) {
    std::ifstream ss(mergedFile);
    if (ss.fail()) {
        cout << "File " << mergedFile << " not found!\n";
        EXIT(EXIT_FAILURE);
    }

    std::string line;
    while (std::getline(ss, line)) {
        std::vector<std::string> result = splitByDelimiter(line, "\t|\t", 2);
        if (result.size() != 2) {
            Debug(Debug::ERROR) << "Invalid name entry!\n";
            EXIT(EXIT_FAILURE);
        }
        TaxID oldId = (TaxID) atoi(result[0].c_str());
        TaxID newId = (TaxID) atoi(result[1].c_str());
        old2new[oldId] = newId;
    }
}


void IndexCreator::addFilesToMerge(string diffIdxFileName, string infoFileName, string posFileName) {
    deltaIdxFileNames.push_back(diffIdxFileName);
    infoFileNames.push_back(infoFileName);
    if (!posFileName.empty()) {
        posFileNames.push_back(posFileName);
    }
}

void IndexCreator::setMergedFileNames(string diffFileName, string infoFileName, string splitFileName, string posFileName) {
    mergedDeltaIdxFileName = diffFileName;
    mergedInfoFileName = infoFileName;
    deltaIdxSplitFileName = splitFileName;
    mergedPosFileName = posFileName;
}

void IndexCreator::updateTaxId2SpeciesTaxId(const string & taxIdListFileName) {
    FILE * taxIdFile;
    if((taxIdFile = fopen(taxIdListFileName.c_str(),"r")) == NULL){
        cout << "Cannot open the taxID list file: " << taxIdListFileName << endl;
        return;
    }

    char taxID[100];
    while(fscanf(taxIdFile,"%s",taxID) == 1) {
        TaxID taxId = atol(taxID);
        TaxonNode const * taxon = taxonomy->taxonNode(taxId);
        TaxID collapseTaxID = getCollapseTaxId(taxId);
        if (taxId == taxon->taxId){
            while (taxon->taxId != collapseTaxID) {
                taxId2speciesId[taxon->taxId] = collapseTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[collapseTaxID] = collapseTaxID;
        } else { // merged
            while (taxon->taxId != collapseTaxID) {
                taxId2speciesId[taxon->taxId] = collapseTaxID;
                taxon = taxonomy->taxonNode(taxon->parentTaxId);
            }
            taxId2speciesId[collapseTaxID] = collapseTaxID;
            taxId2speciesId[taxId] = collapseTaxID;
        }
    }
    fclose(taxIdFile);
    Debug(Debug::INFO) << "Collapse-rank taxonomy IDs are prepared.\n";
}

void IndexCreator::generateIntergenicKmerList(
    _gene *genes, 
    _node *nodes, 
    int numberOfGenes,
    vector <uint64_t> &intergenicKmerList,
    const char *seq) 
{
    if (numberOfGenes == 0) return;

    int k = 12;
    char *kmer = (char *) malloc(sizeof(char) * (k + 1));
    char *reverseKmer = (char *) malloc(sizeof(char) * (k + 1));

    // Use the frame of the first gene for the first intergenic region
    int beginOfFisrtGene = genes[0].begin - 1;
    if (beginOfFisrtGene > k - 1) {
        strncpy(kmer, seq + beginOfFisrtGene - k, k);
        if (nodes[genes[0].start_ndx].strand == 1) {
            intergenicKmerList.push_back(XXH64(kmer, k, 0));
        } else {
            for (int j = k - 1; j >= 0; j--) {
                reverseKmer[k - j - 1] = iRCT[kmer[j]];
            }
            intergenicKmerList.push_back(XXH64(reverseKmer, k, 0));
        }
    }

    //
    for (int i = 0; i < numberOfGenes; i++) {
        strncpy(kmer, seq + genes[i].end, k);
        if (nodes[genes[i].start_ndx].strand == 1) {
            intergenicKmerList.push_back(XXH64(kmer, k, 0));
        } else {
            for (int j = k - 1; j >= 0; j--) {
                reverseKmer[k - j - 1] = iRCT[kmer[j]];
            }
            intergenicKmerList.push_back(XXH64(reverseKmer, k, 0));
        }
    }

    free(reverseKmer);
    free(kmer);
}


bool IndexCreator::compareMinHashList(
    const std::vector<uint64_t>& list1, 
    const std::vector<uint64_t>& list2, 
    size_t length1, 
    size_t length2) 
{      
    if (list1.empty() || list2.empty()) return false;

    // std::set_intersection requires sorted input
    // Since getMinHashList produces a sorted vector, we can use a simple loop
    // or std::set_intersection.
    
    size_t identicalCount = 0;
    size_t i = 0, j = 0;
    
    while (i < list1.size() && j < list2.size()) {
        if (list1[i] < list2[j]) {
            i++;
        } else if (list1[i] > list2[j]) {
            j++;
        } else {
            identicalCount++;
            i++;
            j++;
        }
    }
    float lengthRatio = (float)length2 / (float)length1;
    return identicalCount > (list1.size() * lengthRatio * 0.5f);
}

std::vector<uint64_t> IndexCreator::getMinHashList(const char* seq) {
    std::vector<uint64_t> hashes;
    if (seq == nullptr) return hashes;
    size_t seqLength = strlen(seq);
    size_t kmerLength = 24;
    size_t maxLength = 3000;
    if (seqLength < kmerLength) return hashes;
    // Use a priority queue locally just to maintain the bottom-k
    std::priority_queue<uint64_t> maxHeap; 
    size_t lastStart = seqLength - kmerLength;
    for (size_t i = 0; i <= lastStart; ++i) {
        uint64_t currHash = XXH64(seq + i, kmerLength, 0);
        if (maxHeap.size() < maxLength) {
            maxHeap.push(currHash);
        } else if (currHash < maxHeap.top()) {
            maxHeap.pop();
            maxHeap.push(currHash);
        }
    }
    // Convert heap to sorted vector
    hashes.reserve(maxHeap.size());
    while (!maxHeap.empty()) {
        hashes.push_back(maxHeap.top());
        maxHeap.pop();
    }
    // Heap comes out largest-first, so sort or reverse. 
    // Standard sort is safest for set_intersection logic.
    std::sort(hashes.begin(), hashes.end());
    
    return hashes;
}

void IndexCreator::printSpeciesBatches(){
    for (size_t i = 0; i < spBatches.size(); i++) {
        cout << "Species batch " << i << ": SpeciesID=" << spBatches[i].speciesID << endl;
        cout << "  Genome Num: " << spBatches[i].fastaBatches.size() << "\tTotal Length: " << spBatches[i].spTotalLength << endl;
        cout << "  Rep Genome: " << spBatches[i].repGenomeFasta << "\t" << fastaPaths[spBatches[i].repGenomeFasta] << endl;
        cout << "  FASTA list: " << endl;
        for (size_t j = 0; j < spBatches[i].fastaBatches.size(); j++) {
            cout << "    " << fastaPaths[spBatches[i].fastaBatches[j].whichFasta] << endl;
        }
    }
}
