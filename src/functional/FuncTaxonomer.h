#include "Taxonomer.h"

struct Result {
    int queryId;
    vector<uint32_t> targetRegionId;
    vector<uint32_t> targetUniRefId;
};

class FuncTaxonomer : public Taxonomer {
    protected:
       
    public:
        FuncTaxonomer(LocalParameters & par, NcbiTaxonomy * taxonomy) : Taxonomer(par, taxonomy) {}
        void assignTaxonomyAndFunction(const Buffer<MatchF> & matchList,
                                       std::vector<Query> & queryList,
                                       std::vector<Result> & resultList);
        virtual ~FuncTaxonomer() {}
};