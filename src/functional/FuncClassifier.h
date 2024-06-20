#include "Classifier.h"
#include "FuncTaxonomer.h"
#include "common.h"

class QueryF: public Query{
    vector<uint32_t> targetRegionId;
    vector<uint32_t> targetUniRefId;
    QueryF(int queryId, int classification, float score, float coverage, int hammingDist, int queryLength,
          int queryLength2, int kmerCnt, bool isClassified, bool newSpecies, std::string name)
            : Query(queryId, classification, score, coverage, hammingDist, queryLength, queryLength2, kmerCnt, isClassified, newSpecies, name) {}
             
    QueryF() : Query() {}
};

class FuncClassifier : public Classifier {
    protected:
        FuncTaxonomer * taxonomer;
        void sortMatches(Buffer<MatchF> *matchBuffer);

    public:
        FuncClassifier(LocalParameters & par) : Classifier(par, true) {
            taxonomer = new FuncTaxonomer(par, taxonomy);
        }
        virtual ~FuncClassifier() {}
        void startClassify();
};