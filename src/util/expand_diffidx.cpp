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

int expand_diffidx(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string diffIdxFileName = par.filenames[0];
    size_t bufferSize = 16'777'216; // 1000'000'000;
   
    // RegularGeneticCode geneticCode();
    GeneticCode * geneticCode = new RegularGeneticCode();

    size_t startIdx = par.kmerBegin;
    size_t endIdx = par.kmerEnd;
    size_t idx = 0;

    KmerDbReader * kmerDbReader = new KmerDbReader(diffIdxFileName, 1024 * 1024, 1024 * 1024);
    while (!kmerDbReader->isCompleted()) {
        uint64_t kmer = kmerDbReader->next();
        if (idx >= startIdx && idx < endIdx) {
            Kmer kmerObj(kmer, 0);
            kmerObj.printAA(geneticCode, 12); cout << endl;
            // seqIterator->printKmerInDNAsequence(kmer); cout << "\t";
            // print_binary64(64, kmer); cout << "\n";  
        }
        idx++;        
    }
                   
    delete kmerDbReader;
    delete geneticCode;
    return 0;
}