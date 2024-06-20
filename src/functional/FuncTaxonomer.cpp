# include "FuncTaxonomer.h"

void FuncTaxonomer::assignTaxonomyAndFunction(const Buffer<MatchF> & matchList,
                                              std::vector<Query> & queryList,
                                              std::vector<Result> & resultList) {
    time_t beforeAnalyze = time(nullptr);
    cout << "Analyzing matches ..." << endl;

    // Divide matches into blocks for multi threading
    size_t seqNum = queryList.size();
    MatchBlock *matchBlocks = new MatchBlock[seqNum];
    size_t matchIdx = 0;
    size_t blockIdx = 0;
    uint32_t currentQuery;
    while (matchIdx < matchList.startIndexOfReserve) {
        currentQuery = matchList.buffer[matchIdx].qInfo.sequenceID;
        matchBlocks[blockIdx].id = currentQuery;
        matchBlocks[blockIdx].start = matchIdx;
        while ((currentQuery == matchList.buffer[matchIdx].qInfo.sequenceID) 
            && (matchIdx < matchList.startIndexOfReserve)) {
            ++matchIdx;
        } 
        matchBlocks[blockIdx].end = matchIdx - 1;
        blockIdx++;
    }

    // Process each block
#pragma omp parallel default(none), shared(cout, matchBlocks, matchList, seqNum, queryList, blockIdx)
    {
#pragma omp for schedule(dynamic, 1)
        for (size_t i = 0; i < blockIdx; ++i) {
            chooseBestTaxon(matchBlocks[i].id,
                            matchBlocks[i].start,
                            matchBlocks[i].end,
                            matchList,
                            queryList,
                            par);
        }
    }

    for (size_t i = 0; i < seqNum; i++) {
        ++taxCounts[queryList[i].classification];
    }
    delete[] matchBlocks;
    cout << "Time spent for analyzing: " << double(time(nullptr) - beforeAnalyze) << endl;

    
}