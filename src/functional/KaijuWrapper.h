#ifndef METABULI_KAIJUWRAPPER_H
#define METABULI_KAIJUWRAPPER_H

#include <stdint.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <unordered_map>
#include <algorithm>
#include <string>
#include <deque>
#include <stdexcept>
#include <vector>
#include <cstring>
#include "LocalParameters.h"
// #include "ProducerConsumerQueue/src/ProducerConsumerQueue.hpp"

#include "ReadItem.hpp"
#include "Config.hpp"
#include "util.hpp"
#include "ConsumerThreadp.hpp"

extern "C" {
#include "bwt/bwt.h"
}

class KaijuWrapper //: public ConsumerThreadp
{
private:
    const LocalParameters & par;
    Config * config;
    std::map<char, std::vector<char>> blosum_subst;
    uint8_t nuc2int[256];
    uint8_t compnuc2int[256];
    uint8_t aa2int[256];
    int8_t blosum62diag[20];
    int8_t b62[20][20];
    char codon2aa[256];

    uint8_t codon_to_int(const char* codon);

    unsigned int calcScore(const std::string &, size_t, size_t, int);
    unsigned int calcScore(const std::string &, int);
    unsigned int calcScore(const std::string &);
    void eval_match_scores(SI *si, Fragment * frag, std::vector<SI *> & best_matches_SI, unsigned int best_match_score, std::vector<std::string> & best_matches);

    void ids_from_SI(SI *si, std::set<char *> & match_ids);
    void ids_from_SI_recursive(SI *si, std::set<char *> & match_ids);
    
    Fragment * getNextFragment(
        std::multimap<unsigned int,Fragment *,std::greater<unsigned int>> & fragments,
        unsigned int);
    
    
    
public:
    KaijuWrapper(const LocalParameters & par);
    ~KaijuWrapper();
    void translate(const std::string & cds, std::string & aa);
    void classify_length(std::string sequence, std::set<char *> & match_ids, int & maxLength);
    void classify_greedyblosum(std::string sequence, std::set<char *> & match_ids);
    void addAllMismatchVariantsAtPosSI(const Fragment *,
        std::multimap<unsigned int,Fragment *,std::greater<unsigned int>> & fragments,
        unsigned int,
        unsigned int best_match_score,
        size_t,
        SI *);
    

};

inline uint8_t KaijuWrapper::codon_to_int(const char* codon)  {
 return (uint8_t)(nuc2int[(uint8_t)codon[0]] << 4 | nuc2int[(uint8_t)codon[1]]  << 2 | nuc2int[(uint8_t)codon[2]]);
}

#endif // METABULI_KAIJUWRAPPER_H