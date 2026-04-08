#include "LocalParameters.h"
#include <Command.h>
#include <cstddef>
#include <string>
#include <iostream>
#include "KmerMatcher.h"
#include "common.h"
#include "SeqIterator.h"
#include <cstdint>
#include <chrono>
#include <fcntl.h>
#include "DeltaIdxReader.h"
using namespace std;

int inspect_db(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string dbDir = par.filenames[0];
    loadDbParameters(par, dbDir);
    int kmerFormat = par.kmerFormat;

    MetamerPattern * metamerPattern = nullptr;

    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, "");

    if (kmerFormat == 1) {
        metamerPattern = new LegacyPattern(std::make_unique<RegularGeneticCode>(), 8);
        cout << "Using legacy k-mer format." << endl;
    } else if (!par.customMetamer.empty()) {
        int codeNum = getCodeNum(par.customMetamer);
        if (codeNum == 1) {
            metamerPattern = new SingleCodePattern(par.customMetamer);
        } else if (codeNum > 1) {
            metamerPattern = new MultiCodePattern(par.customMetamer);
        }
    } else if (kmerFormat == 2) {
        bool multicode = false;
        if (multicode) {
            vector<std::unique_ptr<GeneticCode>> geneticCodes;
            geneticCodes.push_back(std::make_unique<RegularGeneticCode>());
            vector<int> codePatterns{0, 0, 0, 0, 0, 0, 0, 0};
            metamerPattern = new MultiCodePattern(std::move(geneticCodes), codePatterns);
        } else {
            std::cout << "Using standard k-mer format." << std::endl;
            metamerPattern = new SingleCodePattern(std::make_unique<RegularGeneticCode>(), 8);
        }
    }
  
    size_t startIdx = par.kmerBegin;
    size_t endIdx = par.kmerEnd;

    DeltaIdxReader * deltaIdxReader = new DeltaIdxReader(dbDir + "/diffIdx", dbDir + "/info", 1024 * 1024 * 4, 1024 * 1024 * 4);
    deltaIdxReader->setReadPosition(DiffIdxSplit{0, 0, 0});

    size_t idx = 0;
    while (!deltaIdxReader->isCompleted()) {
        Kmer kmer = deltaIdxReader->next();
        if (taxonomy->getOriginalTaxID(kmer.tInfo.taxId) == 2567792) {
            metamerPattern->printAA(kmer.value); cout << "\t"; metamerPattern->printDNA(kmer.value); cout << "\t";
            cout <<kmer.tInfo.taxId << "\n";
        }
        idx++;        
    }
                   
    delete deltaIdxReader;
    return 0;
}