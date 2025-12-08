#ifndef METABULI_GENETIC_CODE_H
#define METABULI_GENETIC_CODE_H

#include <iostream>

#define nuc2int(x) (x & 14u)>>1u

class GeneticCode {
    public:
        const std::string atcg = "................................................................"
                                 ".AGCG..GT..G.CN...ACTG.A.T.......agcg..gt..g.cn...actg.a.t......"
                                 "................................................................"
                                 "................................................................";
        const std::string iRCT = "................................................................"
                                 ".TVGH..CD..M.KN...YSAABW.R.......tvgh..cd..m.kn...ysaabw.r......"
                                 "................................................................"
                                 "................................................................";

        int codon2AA[8][8][8];
        int codon2codonIdx[8][8][8];
        std::string aminoacids;
        std::vector<std::vector<std::string>> aa2codon;
        int maxCodonPerAA;
        int alphabetSize;
        int bitPerAA;
        int bitPerCodon;
        uint8_t * hammingLookup;
        explicit GeneticCode() {

        }
        virtual void initialize() {

        }
        virtual ~GeneticCode() {
            delete[] hammingLookup;
        }

        int getAA(const char nuc1, const char nuc2, const char nuc3) const {
            return codon2AA[nuc2int(nuc1)][nuc2int(nuc2)][nuc2int(nuc3)];
        }

        int getCodon(const char nuc1, const char nuc2, const char nuc3) const {
            return codon2codonIdx[nuc2int(nuc1)][nuc2int(nuc2)][nuc2int(nuc3)];
        }
        // virtual int getAA(const char nuc1, const char nuc2, const char nuc3) {

        // }
        // virtual int getCodon(const char nuc1, const char nuc2, const char nuc3) const = 0;
        virtual uint8_t getHammingDist(int aaIdx, int codon1Idx, int codon2Idx) const = 0;
        // virtual uint8_t getHammingDist(uint8_t codon1, uint8_t codon2) const = 0;

};

class RegularGeneticCode final : public GeneticCode {
    public:
        RegularGeneticCode() : GeneticCode() {
            maxCodonPerAA = 6;
            alphabetSize = 21;
            bitPerAA = 5;
            bitPerCodon = 3;
            hammingLookup = new uint8_t[8 * 8] {
                0, 1, 1, 1, 2, 1, 3, 3,
                1, 0, 1, 1, 2, 2, 3, 2,
                1, 1, 0, 1, 2, 2, 2, 3,
                1, 1, 1, 0, 1, 2, 3, 3,
                2, 2, 2, 1, 0, 1, 4, 4,
                1, 2, 2, 2, 1, 0, 4, 4,
                3, 3, 2, 3, 4, 4, 0, 1,
                3, 2, 3, 3, 4, 4, 1, 0                
            };
            
            aminoacids = "ARNDCQEGHILKMFPSTWYVX";

            // A
            codon2AA[3][1][0] = 0;
            codon2AA[3][1][1] = 0;
            codon2AA[3][1][2] = 0; 
            codon2AA[3][1][3] = 0;
            aa2codon.push_back({"GCA", "GCC", "GCT", "GCG"});
            
            // R
            codon2AA[1][3][0] = 1;
            codon2AA[1][3][1] = 1;
            codon2AA[1][3][2] = 1;
            codon2AA[1][3][3] = 1; 
            codon2AA[0][3][0] = 1; 
            codon2AA[0][3][3] = 1;
            aa2codon.push_back({"CGA", "CGC", "CGT", "CGG", "AGG", "AGA"});
            // N
            codon2AA[0][0][2] = 2;
            codon2AA[0][0][1] = 2;
            aa2codon.push_back({"ERR", "AAC", "AAT"});

            // D
            codon2AA[3][0][2] = 3;
            codon2AA[3][0][1] = 3;
            aa2codon.push_back({"ERR", "GAC", "GAT"});
                
            // C
            codon2AA[2][3][2] = 4;
            codon2AA[2][3][1] = 4;
            aa2codon.push_back({"ERR", "TGC", "TGT"});

            // Q
            codon2AA[1][0][0] = 5;
            codon2AA[1][0][3] = 5;
            aa2codon.push_back({"CAA", "ERR", "ERR", "CAG"});

            // E
            codon2AA[3][0][0] = 6;
            codon2AA[3][0][3] = 6;
            aa2codon.push_back({"GAA", "ERR", "ERR", "GAG"});

            // G
            codon2AA[3][3][0] = 7;
            codon2AA[3][3][1] = 7;
            codon2AA[3][3][2] = 7;
            codon2AA[3][3][3] = 7;
            aa2codon.push_back({"GGA", "GGC", "GGT", "GGG"});

            // H
            codon2AA[1][0][2] = 8;
            codon2AA[1][0][1] = 8;
            aa2codon.push_back({"ERR", "CAC", "CAT"});
            
            // I
            codon2AA[0][2][2] = 9;
            codon2AA[0][2][1] = 9;
            codon2AA[0][2][0] = 9;
            aa2codon.push_back({"ATA", "ATC", "ATT"});
            
            // L
            codon2AA[2][2][0] = 10;
            codon2AA[2][2][3] = 10;
            codon2AA[1][2][0] = 10;
            codon2AA[1][2][1] = 10;
            codon2AA[1][2][2] = 10;
            codon2AA[1][2][3] = 10;
            aa2codon.push_back({"CTA", "CTC", "CTT", "CTG", "TTG", "TTA"});
            
            // K
            codon2AA[0][0][0] = 11; 
            codon2AA[0][0][3] = 11;
            aa2codon.push_back({"AAA", "ERR", "ERR", "AAG"});
            
            // M
            codon2AA[0][2][3] = 12;
            aa2codon.push_back({"ERR", "ERR", "ERR", "ATG"});
            
            // F
            codon2AA[2][2][2] = 13;
            codon2AA[2][2][1] = 13;
            aa2codon.push_back({"ERR", "TTC", "TTT"});
            
            // P
            codon2AA[1][1][0] = 14; 
            codon2AA[1][1][1] = 14;
            codon2AA[1][1][2] = 14;
            codon2AA[1][1][3] = 14;
            aa2codon.push_back({"CCA", "CCC", "CCT", "CCG"});
            
            // S
            codon2AA[2][1][0] = 15;
            codon2AA[2][1][1] = 15;
            codon2AA[2][1][2] = 15;
            codon2AA[2][1][3] = 15;
            codon2AA[0][3][2] = 15;
            codon2AA[0][3][1] = 15;
            aa2codon.push_back({"TCA", "TCC", "TCT", "TCG", "ERR", "ERR", "AGT", "AGC"});
            
            // T
            codon2AA[0][1][0] = 16;
            codon2AA[0][1][1] = 16;
            codon2AA[0][1][2] = 16;
            codon2AA[0][1][3] = 16;
            aa2codon.push_back({"ACA", "ACC", "ACT", "ACG"});
            // W
            codon2AA[2][3][3] = 17;
            aa2codon.push_back({"ERR", "ERR", "ERR", "TGG"});
            
            // Y
            codon2AA[2][0][2] = 18;
            codon2AA[2][0][1] = 18;
            aa2codon.push_back({"ERR", "TAC", "TAT"});
            // V
            codon2AA[3][2][0] = 19;
            codon2AA[3][2][1] = 19;
            codon2AA[3][2][2] = 19;
            codon2AA[3][2][3] = 19;
            aa2codon.push_back({"GTA", "GTC", "GTT", "GTG"});
            
            // Stop
            codon2AA[2][0][0] = 20; 
            codon2AA[2][3][0] = 20;
            codon2AA[2][0][3] = 20;
            aa2codon.push_back({"TAA", "ERR", "ERR", "TAG", "ERR", "TGA"});
            
            // triplet code with N's
            for(int i = 0; i < 8; i++){
                for(int j = 0; j < 8; j++){
                    codon2AA[7][i][j] = -1;
                    codon2AA[i][7][j] = -1;
                    codon2AA[i][j][7] = -1;
                    codon2codonIdx[7][i][j] = -1;
                    codon2codonIdx[i][7][j] = -1;
                    codon2codonIdx[i][j][7] = -1;
                }
            }
            
            // For encoding DNA information in k-mer
            for(int i =0; i < 4 ; i++){
                for(int i2 = 0; i2 < 4 ; i2++){
                    codon2codonIdx[i][i2][0] = 0;
                    codon2codonIdx[i][i2][1] = 1;
                    codon2codonIdx[i][i2][2] = 2;
                    codon2codonIdx[i][i2][3] = 3;
                }
            }
            // for Arg
            codon2codonIdx[0][3][3] = 4; //AGG
            codon2codonIdx[0][3][0] = 5;

            // for Leu
            codon2codonIdx[2][2][3] = 4;
            codon2codonIdx[2][2][0] = 5;

            // for Ser
            codon2codonIdx[0][3][2] = 6;
            codon2codonIdx[0][3][1] = 7;

            // for stop codon
            codon2codonIdx[2][3][0] = 5;

        }

        inline uint8_t getHammingDist(int aaIdx, int codon1Idx, int codon2Idx) const override {
            return hammingLookup[codon1Idx * 8 + codon2Idx];
        }

        // uint8_t getHammingDist(uint8_t codon1, uint8_t codon2) const override {
        //     return hammingLookup[codon1 * 8 + codon2];
        // }
};


class ReducedGeneticCode final : public GeneticCode {
    public:
        std::vector<std::vector<uint8_t>> hammingMatrix;
        inline uint8_t getHammingDist(int aaIdx, int codon1Idx, int codon2Idx) const override {
            return hammingLookup[aaIdx * maxCodonPerAA * maxCodonPerAA + codon1Idx * maxCodonPerAA + codon2Idx];
        }

        ReducedGeneticCode() : GeneticCode() {
            bitPerAA = 4;
            bitPerCodon = 3;
            maxCodonPerAA = 7;
            alphabetSize = 16; // 15 A.A. + 1 stop
            hammingLookup = new uint8_t[alphabetSize * maxCodonPerAA * maxCodonPerAA] {
                // 4 means not used
                // A 
                0, 1, 1, 1, 4, 4, 4,
                1, 0, 1, 1, 4, 4, 4,
                1, 1, 0, 1, 4, 4, 4,
                1, 1, 1, 0, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // R
                0, 1, 1, 1, 1, 2, 4,
                1, 0, 1, 1, 2, 2, 4,
                1, 1, 0, 1, 2, 2, 4,
                1, 1, 1, 0, 2, 1, 4,
                1, 2, 2, 2, 0, 1, 4,
                2, 2, 2, 1, 1, 0, 4,
                4, 4, 4, 4, 4, 4, 4,
                // N
                0, 1, 4, 4, 4, 4, 4,
                1, 0, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // D
                0, 1, 4, 4, 4, 4, 4,
                1, 0, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // C
                0, 1, 4, 4, 4, 4, 4,
                1, 0, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // QE
                0, 1, 1, 2, 4, 4, 4,
                1, 0, 2, 1, 4, 4, 4,
                1, 2, 0, 1, 4, 4, 4,
                2, 1, 1, 0, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // G
                0, 1, 1, 1, 4, 4, 4,
                1, 0, 1, 1, 4, 4, 4,
                1, 1, 0, 1, 4, 4, 4,
                1, 1, 1, 0, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // H
                0, 1, 4, 4, 4, 4, 4,
                1, 0, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // IV
                0, 1, 1, 1, 1, 2, 2,
                1, 0, 1, 1, 2, 1, 2,
                1, 1, 0, 1, 2, 2, 1,
                1, 1, 1, 0, 2, 2, 2,
                1, 2, 2, 2, 0, 1, 1,
                2, 1, 2, 2, 1, 0, 1,
                2, 2, 1, 2, 1, 1, 0,
                // ML
                0, 1, 1, 1, 1, 2, 2,
                1, 0, 1, 1, 2, 2, 2,
                1, 1, 0, 1, 2, 2, 2,
                1, 1, 1, 0, 2, 1, 1,
                1, 2, 2, 2, 0, 1, 2,
                2, 2, 2, 1, 1, 0, 1,
                2, 2, 2, 1, 2, 1, 0,
                // K
                0, 1, 4, 4, 4, 4, 4,
                1, 0, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // FYW
                0, 1, 1, 2, 2, 4, 4,
                1, 0, 2, 1, 2, 4, 4,
                1, 2, 0, 1, 2, 4, 4,
                2, 1, 1, 0, 2, 4, 4,
                2, 2, 2, 2, 0, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // P
                0, 1, 1, 1, 4, 4, 4,
                1, 0, 1, 1, 4, 4, 4,
                1, 1, 0, 1, 4, 4, 4,
                1, 1, 1, 0, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // S
                0, 1, 1, 1, 3, 3, 4,
                1, 0, 1, 1, 2, 3, 4,
                1, 1, 0, 1, 3, 2, 4,
                1, 1, 1, 0, 3, 3, 4,
                3, 2, 3, 3, 0, 1, 4,
                3, 3, 2, 3, 1, 0, 4,
                4, 4, 4, 4, 4, 4, 4,
                // T
                0, 1, 1, 1, 4, 4, 4,
                1, 0, 1, 1, 4, 4, 4,
                1, 1, 0, 1, 4, 4, 4,
                1, 1, 1, 0, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                // Stop
                0, 1, 1, 4, 4, 4, 4,
                1, 0, 2, 4, 4, 4, 4,
                1, 2, 0, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4,
                4, 4, 4, 4, 4, 4, 4
            };
            aminoacids = "ARNDCQGHILKFPSTX";
            // A
            codon2AA[3][1][0] = 0; codon2codonIdx[3][1][0] = 0;
            codon2AA[3][1][1] = 0; codon2codonIdx[3][1][1] = 1;
            codon2AA[3][1][2] = 0; codon2codonIdx[3][1][2] = 2;
            codon2AA[3][1][3] = 0; codon2codonIdx[3][1][3] = 3;
            aa2codon.push_back({"GCA", "GCC", "GCT", "GCG"});
            // R
            codon2AA[1][3][0] = 1; codon2codonIdx[1][3][0] = 0;
            codon2AA[1][3][1] = 1; codon2codonIdx[1][3][1] = 1;
            codon2AA[1][3][2] = 1; codon2codonIdx[1][3][2] = 2;
            codon2AA[1][3][3] = 1; codon2codonIdx[1][3][3] = 3;
            codon2AA[0][3][0] = 1; codon2codonIdx[0][3][0] = 4;
            codon2AA[0][3][3] = 1; codon2codonIdx[0][3][3] = 5;
            aa2codon.push_back({"CGA", "CGC", "CGT", "CGG", "AGA", "AGG"});
            // N
            codon2AA[0][0][1] = 2; codon2codonIdx[0][0][1] = 0;
            codon2AA[0][0][2] = 2; codon2codonIdx[0][0][2] = 1;
            aa2codon.push_back({"AAC", "AAT"});
            // D
            codon2AA[3][0][1] = 3; codon2codonIdx[3][0][1] = 0;
            codon2AA[3][0][2] = 3; codon2codonIdx[3][0][2] = 1;
            aa2codon.push_back({"GAC", "GAT"});
            // C
            codon2AA[2][3][2] = 4; codon2codonIdx[2][3][2] = 0;
            codon2AA[2][3][1] = 4; codon2codonIdx[2][3][1] = 1;
            aa2codon.push_back({"TGC", "TGT"});
            // QE
            codon2AA[1][0][0] = 5; codon2codonIdx[1][0][0] = 0;
            codon2AA[1][0][3] = 5; codon2codonIdx[1][0][3] = 1;
            codon2AA[3][0][0] = 5; codon2codonIdx[3][0][0] = 2;
            codon2AA[3][0][3] = 5; codon2codonIdx[3][0][3] = 3;
            aa2codon.push_back({"CAA", "CAG", "GAA", "GAG"});
            // G
            codon2AA[3][3][0] = 6; codon2codonIdx[3][3][0] = 0;
            codon2AA[3][3][1] = 6; codon2codonIdx[3][3][1] = 1;
            codon2AA[3][3][2] = 6; codon2codonIdx[3][3][2] = 2;
            codon2AA[3][3][3] = 6; codon2codonIdx[3][3][3] = 3;
            aa2codon.push_back({"GGA", "GGC", "GGT", "GGG"});
            // H
            codon2AA[1][0][1] = 7; codon2codonIdx[1][0][1] = 0;
            codon2AA[1][0][2] = 7; codon2codonIdx[1][0][2] = 1;
            aa2codon.push_back({"CAC", "CAT"});
            // IV
            codon2AA[3][2][0] = 8; codon2codonIdx[3][2][0] = 0;
            codon2AA[3][2][1] = 8; codon2codonIdx[3][2][1] = 1;
            codon2AA[3][2][2] = 8; codon2codonIdx[3][2][2] = 2;
            codon2AA[3][2][3] = 8; codon2codonIdx[3][2][3] = 3;
            codon2AA[0][2][0] = 8; codon2codonIdx[0][2][0] = 4;
            codon2AA[0][2][1] = 8; codon2codonIdx[0][2][1] = 5;
            codon2AA[0][2][2] = 8; codon2codonIdx[0][2][2] = 6;
            aa2codon.push_back({"GTA", "GTC", "GTT", "GTG", "ATA", "ATC", "ATT"});
            // ML
            codon2AA[1][2][0] = 9; codon2codonIdx[1][2][0] = 0;
            codon2AA[1][2][1] = 9; codon2codonIdx[1][2][1] = 1;
            codon2AA[1][2][2] = 9; codon2codonIdx[1][2][2] = 2;
            codon2AA[1][2][3] = 9; codon2codonIdx[1][2][3] = 3;
            codon2AA[2][2][0] = 9; codon2codonIdx[2][2][0] = 4;
            codon2AA[2][2][3] = 9; codon2codonIdx[2][2][3] = 5;
            codon2AA[0][2][3] = 9; codon2codonIdx[0][2][3] = 6;
            aa2codon.push_back({"CTA", "CTC", "CTT", "CTG", "TTA", "TTG", "ATG"});
            // K
            codon2AA[0][0][0] = 10; codon2codonIdx[0][0][0] = 0;
            codon2AA[0][0][3] = 10; codon2codonIdx[0][0][3] = 1;
            aa2codon.push_back({"AAA", "AAG"});
            // FYW
            codon2AA[2][0][1] = 11; codon2codonIdx[2][0][1] = 0;
            codon2AA[2][0][2] = 11; codon2codonIdx[2][0][2] = 1;
            codon2AA[2][2][1] = 11; codon2codonIdx[2][2][1] = 2;
            codon2AA[2][2][2] = 11; codon2codonIdx[2][2][2] = 3;
            codon2AA[2][3][3] = 11; codon2codonIdx[2][3][3] = 4;
            aa2codon.push_back({"TAC", "TAT", "TTC", "TTT", "TGG"});
            // P
            codon2AA[1][1][0] = 12; codon2codonIdx[1][1][0] = 0;
            codon2AA[1][1][1] = 12; codon2codonIdx[1][1][1] = 1;
            codon2AA[1][1][2] = 12; codon2codonIdx[1][1][2] = 2;
            codon2AA[1][1][3] = 12; codon2codonIdx[1][1][3] = 3;
            aa2codon.push_back({"CCA", "CCC", "CCT", "CCG"});
            // S
            codon2AA[2][1][0] = 13; codon2codonIdx[2][1][0] = 0;
            codon2AA[2][1][1] = 13; codon2codonIdx[2][1][1] = 1;
            codon2AA[2][1][2] = 13; codon2codonIdx[2][1][2] = 2;
            codon2AA[2][1][3] = 13; codon2codonIdx[2][1][3] = 3;
            codon2AA[0][3][1] = 13; codon2codonIdx[0][3][1] = 4;
            codon2AA[0][3][2] = 13; codon2codonIdx[0][3][2] = 5;
            aa2codon.push_back({"TCA", "TCC", "TCT", "TCG", "AGC", "AGT"});
            // T
            codon2AA[0][1][0] = 14; codon2codonIdx[0][1][0] = 0;
            codon2AA[0][1][1] = 14; codon2codonIdx[0][1][1] = 1;
            codon2AA[0][1][2] = 14; codon2codonIdx[0][1][2] = 2;
            codon2AA[0][1][3] = 14; codon2codonIdx[0][1][3] = 3;
            aa2codon.push_back({"ACA", "ACC", "ACT", "ACG"});
            // Stop
            codon2AA[2][0][0] = 15; codon2codonIdx[2][0][0] = 0;
            codon2AA[2][0][3] = 15; codon2codonIdx[2][0][3] = 1;
            codon2AA[2][3][0] = 15; codon2codonIdx[2][3][0] = 2;
            aa2codon.push_back({"TAA", "TAG", "TGA"});
            // Triplet code with N's
            for (int i = 0; i < 8; i++) {
                for (int j = 0; j < 8; j++) {
                    codon2AA[7][i][j] = -1;
                    codon2AA[i][7][j] = -1;
                    codon2AA[i][j][7] = -1;
                    codon2codonIdx[7][i][j] = -1;
                    codon2codonIdx[i][7][j] = -1;
                    codon2codonIdx[i][j][7] = -1;
                }
            }
        }
};

#endif //METABULI_GENETIC_CODE_H