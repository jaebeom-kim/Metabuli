#include "Taxonomer.h"
#include "BitManipulateMacros.h"
#include "Match.h"
#include "NcbiTaxonomy.h"
#include "common.h"
#include "printBinary.h"
#include <cstdint>
#include <sys/types.h>
#include <unordered_map>

template <typename MatchType>
Taxonomer<MatchType>::Taxonomer(const LocalParameters &par, TaxonomyWrapper *taxonomy, const MetamerPattern *metamerPattern) :  
    par(par), 
    taxonomy(taxonomy), 
    metamerPattern(metamerPattern), 
    kmerLen(metamerPattern->kmerLen), 
    windowSize(metamerPattern->windowSize),
    windowMask((1U << windowSize) - 1)
{
    if (par.substitutionMatrices.empty()) {
        substitutionMatrix = nullptr;
    } else {
        substitutionMatrix = new SubstitutionMatrix(par.scoringMatrixFile.values.aminoacid().c_str(), 2.0, par.scoreBias);
    }
    // Parameters
    accessionLevel = par.accessionLevel;
    minConsCnt = (size_t) par.minConsCnt;
    minConsCntEuk = (size_t) par.minConsCntEuk;
    eukaryotaTaxId = taxonomy->getEukaryotaTaxID();
    tieRatio = par.tieRatio;


    maxCodonShift = par.maxShift;
    dnaShift = 3;
    if (par.syncmer) {
        dnaShift = (kmerLen - par.smerLen) * 3;
    }
    // if (par.syncmer) {
    //     dnaShift = (8 - par.smerLen) * 3;
    //     maxCodonShift = 8 - par.smerLen;
    // }

    if (par.seqMode == 1 || par.seqMode == 2) {
        denominator = 100;
    } else {
        denominator = 1000;
    }

    // chooseBestTaxon
    taxCnt.reserve(4096);

    // getBestSpeciesMatches
    matchPaths.reserve(4096);
    combinedMatchPaths.reserve(4096);
    maxSpecies.reserve(4096);

    // lowerRankClassification
    cladeCnt.reserve(4096);

    // filterRedundantMatches
    arraySize_filterRedundantMatches = 4096;
    // bestMatchForQuotient = new const Match*[arraySize_filterRedundantMatches]();
    bestMatchTaxIdForQuotient = new TaxID[arraySize_filterRedundantMatches]();
    minHammingForQuotient = new uint8_t[arraySize_filterRedundantMatches]();

    // Output
    taxCounts.reserve(4096);

    if (par.maxEValue > 0) {
        logMaxEValue = std::log(par.maxEValue);
        useEvalueFilter = true;
    } else {
        useEvalueFilter = false;
    }
    dbSize = par.dbTotalLength == 0 ? readDbSize(par.filenames[1 + (par.seqMode == 2)]) : par.dbTotalLength;
    // size_t tempSize = 562762599 * 2; // for nr database
}

template <typename MatchType>
Taxonomer<MatchType>::~Taxonomer() {
    delete substitutionMatrix;
    // delete[] bestMatchForQuotient;
    delete[] bestMatchTaxIdForQuotient;
    delete[] minHammingForQuotient;
    delete evaluer;
}

template <typename MatchType>
void Taxonomer<MatchType>::chooseBestTaxon(
    uint32_t currentQuery,
    size_t offset,
    size_t end,
    const MatchType *matchList,
    vector<Query> & queryList) 
{
    // cout << "Current query: " << currentQuery << endl;
    // for (size_t i = offset; i < end+1; i ++) {
    //     matchList[i].printMatch();
    // }
    TaxonScore speciesScore;
    std::pair<size_t, size_t> bestSpeciesRange;
    speciesScore = getBestSpeciesMatches(bestSpeciesRange,
                                         matchList,
                                         end,
                                         offset,                        
                                         queryList[currentQuery]);
    
    // If there is no proper species for current query, it is un-classified.
    if (speciesScore.score.idScore == 0 || speciesScore.score.idScore < par.minScore) {
        queryList[currentQuery].classification = 0;
        queryList[currentQuery].idScore = speciesScore.score.idScore;
        queryList[currentQuery].subScore = speciesScore.score.subScore;
        queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        return;
    }

    // If there are two or more good species level candidates, find the LCA.
    if (speciesScore.LCA) {
        queryList[currentQuery].classification = speciesScore.taxId;
        queryList[currentQuery].idScore = speciesScore.score.idScore;
        queryList[currentQuery].subScore = speciesScore.score.subScore;
        queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        return;
    }

    // Filter redundant matches
    taxCnt.clear();
    filterRedundantMatches(matchList,
                           bestSpeciesRange,
                           taxCnt,
                           queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    for (auto & tax : taxCnt) {
      queryList[currentQuery].taxCnt[tax.first] = tax.second;    
    }    


    if constexpr (std::is_same_v<MatchType, MatchWithPos>) {

        auto & speciesBins = sp2coverage[speciesScore.taxId];
        if (speciesBins.empty()) { speciesBins.resize(65536, 0); }
        uint8_t * bins = speciesBins.data();

        for (size_t i = speciesScore.matchPathRange.first; i < speciesScore.matchPathRange.second; i ++) {
            const vector<const MatchType*> & currentChain = combinedMatchPaths[i].chain;
            for (size_t j = 0; j < currentChain.size(); j ++) {
                const uint16_t binId = currentChain[j]->posId; 
                if (bins[binId] < 255) {
                    bins[binId]++;
                }
            }
        }
    }

    // If score is not enough, classify to the parent of the selected species
    if (speciesScore.score.idScore < par.minSpScore) {
      queryList[currentQuery].classification = taxonomy->taxonNode(
              taxonomy->getTaxIdAtRank(speciesScore.taxId, "species"))->parentTaxId;
      queryList[currentQuery].idScore = speciesScore.score.idScore;
      queryList[currentQuery].subScore = speciesScore.score.subScore;
      queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
      queryList[currentQuery].hammingDist = speciesScore.hammingDist;
      return;
    }
    // Store classification results
    queryList[currentQuery].idScore = speciesScore.score.idScore;
    queryList[currentQuery].subScore = speciesScore.score.subScore;
    queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
    queryList[currentQuery].hammingDist = speciesScore.hammingDist;

    if (!par.em) {
        queryList[currentQuery].classification
            = lowerRankClassification(
                taxCnt,
                speciesScore.taxId,
                queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2);
    } else {
        queryList[currentQuery].classification = speciesScore.taxId;
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::filterRedundantMatches(
    const MatchType *matchList,
    const std::pair<size_t, size_t> & bestSpeciesRange,
    unordered_map<TaxID, unsigned int> & taxCnt,
    int queryLength) 
{    
    // Determine the maximum quotient we need to handle
    size_t maxQuotient = (queryLength + 3) / dnaShift;
    
    ensureArraySize(maxQuotient + 1);

    std::fill_n(bestMatchTaxIdForQuotient, maxQuotient + 1, 0);

    for (size_t i = bestSpeciesRange.first; i < bestSpeciesRange.second; i ++) {
        size_t currQuotient = matchList[i].qKmer.qInfo.pos / dnaShift;
        uint8_t hamming = metamerPattern->hammingDistSum(matchList[i].qKmer.value, matchList[i].tKmer.value, kmerLen, true);

        if (bestMatchTaxIdForQuotient[currQuotient] == 0) {
            // bestMatchForQuotient[currQuotient] = matchList + i;
            bestMatchTaxIdForQuotient[currQuotient] = matchList[i].tKmer.tInfo.taxId;
            minHammingForQuotient[currQuotient] = hamming;
        } else {
            if (hamming < minHammingForQuotient[currQuotient]) {
                // bestMatchForQuotient[currQuotient] = matchList + i;
                bestMatchTaxIdForQuotient[currQuotient] = matchList[i].tKmer.tInfo.taxId;
                minHammingForQuotient[currQuotient] = hamming;
            } else if (hamming == minHammingForQuotient[currQuotient]) {
                bestMatchTaxIdForQuotient[currQuotient] = taxonomy->LCA(
                    bestMatchTaxIdForQuotient[currQuotient], matchList[i].tKmer.tInfo.taxId);
            }
        }
    }

    for (size_t i = 0; i <= maxQuotient; ++i) {
        if (bestMatchTaxIdForQuotient[i] != 0) {
            taxCnt[bestMatchTaxIdForQuotient[i]]++;
        }
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::printSpeciesMatches(
    const MatchType *matchList,
    const std::pair<size_t, size_t> & bestSpeciesRange) 
{
    for (size_t i = bestSpeciesRange.first; i < bestSpeciesRange.second; i ++) {
        matchList[i].printMatch();
    }
}

template <typename MatchType>
TaxID Taxonomer<MatchType>::lowerRankClassification(const unordered_map<TaxID, unsigned int> & taxCnt, TaxID spTaxId, int queryLength) {
    unsigned int minSubSpeciesMatch = ((queryLength - 1)/denominator);
    cladeCnt.clear();
    getSpeciesCladeCounts(taxCnt, cladeCnt, spTaxId);
    if (accessionLevel == 2) { // Don't do accession-level classification
        // Remove leaf nodes
        for (auto it = cladeCnt.begin(); it != cladeCnt.end(); it++) {
            TaxonNode const * taxon = taxonomy->taxonNode(it->first);
            if (strcmp(taxonomy->getString(taxon->rankIdx), "") == 0 || strcmp(taxonomy->getString(taxon->rankIdx), "accession") == 0) {
                // Remove current node from its parent's children list
                cladeCnt[taxon->parentTaxId].children.erase(find(cladeCnt[taxon->parentTaxId].children.begin(),
                                                                 cladeCnt[taxon->parentTaxId].children.end(),
                                                                 it->first));
            } 
        }
        return BFS(cladeCnt, spTaxId, minSubSpeciesMatch);
    } else {
        return BFS(cladeCnt, spTaxId, minSubSpeciesMatch);
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::getSpeciesCladeCounts(const unordered_map<TaxID, unsigned int> &taxCnt,
                                                 unordered_map<TaxID, TaxonCounts> & cladeCount,
                                                 TaxID speciesTaxID) {
    for (auto it = taxCnt.begin(); it != taxCnt.end(); ++it) {
        TaxonNode const * taxon = taxonomy->taxonNode(it->first);
        cladeCount[taxon->taxId].taxCount = it->second;
        cladeCount[taxon->taxId].cladeCount += it->second;
        while (taxon->taxId != speciesTaxID) {
            if (find(cladeCount[taxon->parentTaxId].children.begin(),
                     cladeCount[taxon->parentTaxId].children.end(),
                     taxon->taxId) == cladeCount[taxon->parentTaxId].children.end()) {
                cladeCount[taxon->parentTaxId].children.push_back(taxon->taxId);
            }
            cladeCount[taxon->parentTaxId].cladeCount += it->second;
            taxon = taxonomy->taxonNode(taxon->parentTaxId);
        }
    }
}

template <typename MatchType>
TaxID Taxonomer<MatchType>::BFS(const unordered_map<TaxID, TaxonCounts> & cladeCnt, TaxID root, unsigned int maxCnt) {
    unsigned int maxCnt2 = maxCnt;
    if (cladeCnt.at(root).children.empty()) { // root is a leaf
        return root;
    }
    unsigned int currentCnt;
    vector<TaxID> bestChildren;
    for (auto it = cladeCnt.at(root).children.begin(); it != cladeCnt.at(root).children.end(); it++) {
        currentCnt = cladeCnt.at(*it).cladeCount;
        if (currentCnt > maxCnt) {
            bestChildren.clear();
            bestChildren.push_back(*it);
            maxCnt = currentCnt;
        } else if (currentCnt == maxCnt) {
            bestChildren.push_back(*it);
        }
    }
    if (bestChildren.size() == 1) {
        return BFS(cladeCnt, bestChildren[0], maxCnt2);
    } else {
        return root;
    }
}

template <typename MatchType>
TaxonScore Taxonomer<MatchType>::getBestSpeciesMatches(std::pair<size_t, size_t> & bestSpeciesRange,
                                                       const MatchType *matchList,
                                                       size_t end,
                                                       size_t offset,
                                                       Query & query) {
    matchPaths.clear();
    combinedMatchPaths.clear();
    pair<size_t, size_t> bestMatchPathRange;
    vector<pair<TaxID, MatchScore>> sp2score;
    int queryLength = query.queryLength + query.queryLength2;
    if (par.printLog==1) {
        cout << "## " << query.name << " ##" << endl;
    }
    TaxonScore bestScore;
    MatchScore bestSpScore;
    size_t i = offset;
    size_t meaningfulSpecies = 0;
    if (par.printLog==1) {
        cout << "## " << query.name << " ##" << endl;
        for (size_t j = offset; j < end + 1; j++) {
            metamerPattern->printDNA(matchList[j].qKmer.value);
            cout << " " ;
            metamerPattern->printAA(matchList[j].qKmer.value);
            cout << " | " ;
            metamerPattern->printDNA(matchList[j].tKmer.value);
            cout << " " ;
            metamerPattern->printAA(matchList[j].tKmer.value);
            cout << " | " << taxonomy->getOriginalTaxID(matchList[j].tKmer.tInfo.taxId) << " " << taxonomy->getOriginalTaxID(matchList[j].tKmer.tInfo.speciesId) << " ";
            cout << (int) matchList[j].qKmer.qInfo.frame << " " << (int) matchList[j].qKmer.qInfo.pos << endl;
        }
    }
    while (i  < end + 1) {
        TaxID currentSpecies = matchList[i].tKmer.tInfo.speciesId;
        size_t start = i;
        size_t previousPathSize = matchPaths.size();
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].tKmer.tInfo.speciesId) {
            uint8_t curFrame = matchList[i].qKmer.qInfo.frame;
            size_t start = i;
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].tKmer.tInfo.speciesId && curFrame == matchList[i].qKmer.qInfo.frame) {
                i ++;
            }
            if (i - start > 1) {
                if (windowSize == kmerLen) {
                    getMatchPaths_lookbackDP(matchList + start, i - start, matchPaths, currentSpecies);
                } else {
                    getSpacedMatchPaths_lookbackDP(matchList + start, i - start, matchPaths, currentSpecies);
                }
            }
        }
        size_t pathSize = matchPaths.size();
        // Combine MatchPaths
        if (par.printLog==1) {
            if (pathSize > previousPathSize) {
                cout << "Current species: " << taxonomy->getOriginalTaxID(currentSpecies) << " " << currentSpecies << endl;
                for (size_t kk = previousPathSize; kk < matchPaths.size(); kk++) {
                    matchPaths[kk].printMatchPath();
                }
            }
        }
        if (pathSize > previousPathSize) {
            size_t spStart = combinedMatchPaths.size();
            MatchScore score = combineMatchPaths(matchPaths, previousPathSize, combinedMatchPaths, combinedMatchPaths.size(), queryLength);
            size_t spEnd = combinedMatchPaths.size();
            if (par.printLog==1) {   
                cout << "Combined score: " << score.idScore << " " << score.subScore << endl;
            }
            score.idScore = score.idScore / queryLength;
            score.idScore = min(score.idScore, 1.0f);
            if ((score.idScore < par.minScore) || (score.logE > logMaxEValue)) {
                continue;
            }

            sp2score.emplace_back(currentSpecies, score);
            if (score.idScore > 0.f) {
                meaningfulSpecies++;
            }
            if (score.isLargerThan(bestSpScore, par.scoreMode)) {
                bestSpScore = score;
                bestSpeciesRange = make_pair(start, i);
                bestMatchPathRange = make_pair(spStart, spEnd);
            }
        }
    }
    
    // If there are no meaningful species
    if (meaningfulSpecies == 0) {
        bestScore.score = {0.0f, 0.0f};
        return bestScore;
    }

    // if (par.em && !sp2score.empty()) {
    //     sort(sp2score.begin(), sp2score.end(),
    //          [](const pair<TaxID, float> &a, const pair<TaxID, float> &b) {
    //              return a.second > b.second;
    //          });
    //     query.topSpeciesId = sp2score[0].first;
    //     for (size_t i = 0; i < 10 && i < sp2score.size(); i++) {
    //         query.species2Score.emplace_back(sp2score[i].first, sp2score[i].second * sp2score[i].second);
    //     }
    // }

    maxSpecies.clear();
    for (size_t i = 0; i < sp2score.size(); i++) {
        if (sp2score[i].second.isLargerThan(bestSpScore * tieRatio, par.scoreMode)) {
            maxSpecies.push_back(sp2score[i].first);
            bestScore.score += sp2score[i].second;   
        }
    }
    
    // More than one species --> LCA    
    if (maxSpecies.size() > 1) {
        bestScore.LCA = true;
        bestScore.taxId = taxonomy->LCA(maxSpecies)->taxId;
        bestScore.score.idScore /= maxSpecies.size();
        bestScore.score.subScore /= maxSpecies.size();
        bestScore.score.logE = bestSpScore.logE;
        return bestScore;
    }
    
    // One species
    bestScore.taxId = maxSpecies[0];
    bestScore.score.logE = bestSpScore.logE;
    bestScore.matchPathRange = bestMatchPathRange;
    
    return bestScore;                                  
}

template <typename MatchType>
void Taxonomer<MatchType>::sortMatchPath(
    std::vector<MatchPath<MatchType>> & matchPaths, 
    size_t i) 
{
    if (par.scoreMode == 0) {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath<MatchType> &a, const MatchPath<MatchType> &b) {
           if (a.score.idScore != b.score.idScore) {
             return a.score.idScore > b.score.idScore;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });
    } else if (par.scoreMode == 1) {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath<MatchType> &a, const MatchPath<MatchType> &b) {
           if (a.score.subScore != b.score.subScore) {
             return a.score.subScore > b.score.subScore;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });
    } else {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath<MatchType> &a, const MatchPath<MatchType> &b) {
           float aTotal = a.score.idScore + a.score.subScore;
           float bTotal = b.score.idScore + b.score.subScore;
           if (aTotal != bTotal) {
             return aTotal > bTotal;
           }
           if (a.hammingDist != b.hammingDist) {
             return a.hammingDist < b.hammingDist;
           }
           return a.start > b.start;
         });
    }
}

template <typename MatchType>
MatchScore Taxonomer<MatchType>::combineMatchPaths(
    vector<MatchPath<MatchType>> & matchPaths,
    size_t matchPathStart,
    vector<MatchPath<MatchType>> & combinedMatchPaths,
    size_t combMatchPathStart,
    int queryLength) 
{   
    // Combine matchPaths
    // 1. Add the matchPath with the highest score to combinedMatchPaths
    // 2. Add the matchPath with the highest score that is not overlapped with the matchPath in combinedMatchPaths
    // 3. Repeat 2 until no matchPath can be added
    sortMatchPath(matchPaths, matchPathStart);
    MatchScore score;
    int spanLength = 0;
    for (size_t i = matchPathStart; i < matchPaths.size(); i++) {  
        if (combMatchPathStart == combinedMatchPaths.size()) {
            combinedMatchPaths.push_back(matchPaths[i]);
            score += matchPaths[i].score;
            spanLength += matchPaths[i].end - matchPaths[i].start + 1;
        } else {
            bool isOverlapped = false;
            for (size_t j = combMatchPathStart; j < combinedMatchPaths.size(); j++) {
                if (isMatchPathOverlapped(matchPaths[i], combinedMatchPaths[j])) { // overlap!
                    int overlappedLength = min(matchPaths[i].end, combinedMatchPaths[j].end) 
                                            - max(matchPaths[i].start, combinedMatchPaths[j].start) + 1;
                    if (overlappedLength == matchPaths[i].end - matchPaths[i].start + 1) { // Current path is completely overlapped
                        isOverlapped = true;
                        break;
                    }
                    
                    if (overlappedLength < 24) {
                        trimMatchPath(matchPaths[i], combinedMatchPaths[j], overlappedLength);

                        if (matchPaths[i].startMatch == nullptr || matchPaths[i].endMatch == nullptr) { // Current path is completely trimmed
                            isOverlapped = true;
                            break;
                        }
                        continue;
                    } else {
                        isOverlapped = true;
                        break;
                    }
                    if (matchPaths[i].end - matchPaths[i].start + 1 < windowSize * 3) { // Current path trimmed too much 
                        isOverlapped = true;
                        break;
                    }  
                } 
            }
            if (!isOverlapped) {
                score += matchPaths[i].score;     
                spanLength += matchPaths[i].end - matchPaths[i].start + 1;      
                combinedMatchPaths.push_back(std::move(matchPaths[i]));
            }
        }
    }

    score.logE = computeLEM_logE(
        queryLength,
        spanLength,
        dbSize,
        score.logP);

    return score;
}

template <typename MatchType>
bool Taxonomer<MatchType>::isMatchPathOverlapped(
    const MatchPath<MatchType> & matchPath1,
    const MatchPath<MatchType> & matchPath2) {
    return !((matchPath1.end < matchPath2.start) || (matchPath2.end < matchPath1.start));                                       
}

template <typename MatchType>
void Taxonomer<MatchType>::trimMatchPath(
    MatchPath<MatchType> & newPath, 
    const MatchPath<MatchType> & path2, 
    int overlapLength) 
{
    const int overlapCodons = overlapLength / 3;
    const int overlapRemainder = overlapLength % 3;

    bool isTrimEnd = (newPath.start < path2.start);
    const auto* targetMatch = isTrimEnd ? newPath.endMatch : newPath.startMatch;
    
    if (isTrimEnd) {
        newPath.end = path2.start - 1;
        while (!newPath.chain.empty() && newPath.chain.back()->qKmer.qInfo.pos > newPath.end) {
            newPath.chain.pop_back();
        }
    } else {
        newPath.start = path2.end + 1;
        auto it = newPath.chain.begin();
        while (it != newPath.chain.end() && (*it)->qKmer.qInfo.pos < newPath.start) {
            ++it;
        }
        newPath.chain.erase(newPath.chain.begin(), it);
    }

    if (!newPath.chain.empty()) {
        newPath.startMatch = newPath.chain.front();
        newPath.endMatch = newPath.chain.back();
    } else {
        // Edge Case: The overlap was so massive it completely consumed the path
        newPath.startMatch = nullptr;
        newPath.endMatch = nullptr;
        return;
    }

    bool isForwardFrame = (targetMatch->qKmer.qInfo.frame < 3);
    bool strandFlag = (isTrimEnd == isForwardFrame);

    newPath.hammingDist = std::max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
        targetMatch->qKmer.value,
        targetMatch->tKmer.value,
        overlapCodons,
        strandFlag));

    newPath.score -= metamerPattern->calMatchScore(
        targetMatch->qKmer.value,
        targetMatch->tKmer.value,
        overlapCodons,
        *substitutionMatrix,
        strandFlag);

    newPath.score.idScore -= overlapRemainder;
}



// void Taxonomer::trimMatchPath2(
//     MatchPath & newPath, 
//     const MatchPath & path2, 
//     int overlapLength) 
// {
//     const int shift = overlapLength / 3;
//     if (newPath.start < path2.start) { 
//         newPath.end = path2.start - 1;
//         if (newPath.endMatch->qKmer.qInfo.frame < 3) {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value,
//                 overlapLength/3,
//                 true));
            
//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value,
//                 newPath.historyMask & ((1U << shift) - 1),
//                 *substitutionMatrix);    

//             newPath.score.idScore -= overlapLength % 3;
//         } else {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value, 
//                 overlapLength/3,
//                 false));
            
//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.endMatch->qKmer.value,
//                 newPath.endMatch->tKmer.value,
//                 newPath.historyMask & (~0U << (32 - shift)),
//                 *substitutionMatrix);    

//             newPath.score.idScore -= overlapLength % 3;
//         }
//     } else {
//         newPath.start = path2.end + 1;
//         if (newPath.startMatch->qKmer.qInfo.frame < 3) {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value, 
//                 overlapLength/3,
//                 false));
            
//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value,
//                 newPath.firstHistoryMask & (~0U << (32 - shift)),
//                 *substitutionMatrix); 
//             newPath.score.idScore -= overlapLength % 3;
//         } else {
//             newPath.hammingDist = max(0, newPath.hammingDist - metamerPattern->hammingDistSum(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value, 
//                 overlapLength/3,
//                 true));

//             newPath.score -= metamerPattern->calMatchScore(
//                 newPath.startMatch->qKmer.value,
//                 newPath.startMatch->tKmer.value,
//                 newPath.firstHistoryMask & ((1U << shift) - 1),
//                 *substitutionMatrix);
//             newPath.score.idScore -= overlapLength % 3;
//         }
//     }
// }

template <typename MatchType>
MatchPath<MatchType> Taxonomer<MatchType>::makeMatchPath(
    const MatchType * match) 
{
    return MatchPath<MatchType>(
        match, 
        metamerPattern->calMatchScore(
            match->qKmer.value, 
            match->tKmer.value, 
            metamerPattern->windowSize, 
            *substitutionMatrix,
            true),
        metamerPattern->hammingDistSum(
            match->qKmer.value, 
            match->tKmer.value),
        metamerPattern->windowSize * 3);
}


template <typename MatchType>
void Taxonomer<MatchType>::getMatchPaths_lookbackDP(
    const MatchType * matchList,
    size_t matchNum,
    vector<MatchPath<MatchType>> & filteredMatchPaths,
    TaxID speciesId) 
{
    bool isForward = matchList[0].qKmer.qInfo.frame < 3;
    int MIN_DEPTH = (int) minConsCnt;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk;
    }

    connectedToNext.resize(matchNum);
    localMatchPaths.resize(matchNum);

    for (size_t i = 0; i < matchNum; i++) {
        localMatchPaths[i] = makeMatchPath(matchList + i);
        connectedToNext[i] = false;
    }
    
    for (size_t i = 1; i < matchNum; ++i) {
        uint32_t currentPos = matchList[i].qKmer.qInfo.pos;

        int bestParentIdx = -1;
        MatchScore bestScore;

        for (int j = i - 1; j >= 0; --j) {
            uint32_t prevPos = matchList[j].qKmer.qInfo.pos;
            int shift = (currentPos - prevPos) / 3;
            if (shift > maxCodonShift) {
                break;
            }
            if (shift <= 0) {
                continue;
            }

            auto tKmerPrev = isForward ? matchList[j].tKmer.value : matchList[i].tKmer.value;
            auto tKmerCurr = isForward ? matchList[i].tKmer.value : matchList[j].tKmer.value;
            if (metamerPattern->checkOverlap(tKmerPrev, tKmerCurr, shift)) {
                connectedToNext[j] = true;
                MatchScore newScore = localMatchPaths[j].score +
                    metamerPattern->calMatchScore(
                        matchList[i].qKmer.value, 
                        matchList[i].tKmer.value, 
                        shift,
                        *substitutionMatrix,
                        isForward);

                if (newScore.isLargerThan(bestScore, par.scoreMode)) {
                    bestScore = newScore;
                    bestParentIdx = j;
                }
            }
        }

        if (bestParentIdx != -1) {
            int shift = (currentPos - matchList[bestParentIdx].qKmer.qInfo.pos) / 3;

            // Update the current node's path data using the best parent
            localMatchPaths[i].start = localMatchPaths[bestParentIdx].start;
            localMatchPaths[i].score = bestScore;
            localMatchPaths[i].hammingDist = localMatchPaths[bestParentIdx].hammingDist + 
                metamerPattern->hammingDistSum(
                    matchList[i].qKmer.value, 
                    matchList[i].tKmer.value, 
                    shift, 
                    isForward);
                
            localMatchPaths[i].depth = localMatchPaths[bestParentIdx].depth + shift;
            localMatchPaths[i].startMatch = localMatchPaths[bestParentIdx].startMatch;
            localMatchPaths[i].prevMatchIdx = bestParentIdx;   
        }
    }

    // 4. Harvest Terminal Paths
    for (size_t i = 0; i < matchNum; ++i) {
        if (!connectedToNext[i] && localMatchPaths[i].depth >= MIN_DEPTH) {
            int currIdx = i;
            while (currIdx != -1) {
                localMatchPaths[i].chain.push_back(matchList + currIdx);
                currIdx = localMatchPaths[currIdx].prevMatchIdx;
            }
            std::reverse(localMatchPaths[i].chain.begin(), localMatchPaths[i].chain.end());

            filteredMatchPaths.push_back(std::move(localMatchPaths[i]));
        }
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::makeMatchPath(
    const MatchType * match,
    size_t index)
{
    localMatchPaths[index] = MatchPath(match, metamerPattern->windowSize * 3);
    localMatchPaths[index].score = metamerPattern->calMatchScore(
                                        match->qKmer.value, 
                                        match->tKmer.value,
                                        *substitutionMatrix);

    localMatchPaths[index].hammingDist = metamerPattern->hammingDistSum(
                                            match->qKmer.value, 
                                            match->tKmer.value);
    localMatchPaths[index].historyMask = metamerPattern->spaceMask;
}


template <typename MatchType>
void Taxonomer<MatchType>::getSpacedMatchPaths_lookbackDP(
    const MatchType * matchList,
    size_t matchNum,
    vector<MatchPath<MatchType>> & filteredMatchPaths,
    TaxID speciesId) 
{
    if (matchNum == 0) return;

    uint64_t frame = matchList[0].qKmer.qInfo.frame;
    bool isForward = (frame < 3); 
    
    int MIN_DEPTH = (int) minConsCnt;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_DEPTH = (int) minConsCntEuk;
    }

    connectedToNext.assign(matchNum, false);
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    // 1. Initialize all match paths sequentially
    for (size_t i = 0; i < matchNum; ++i) {
        makeMatchPath(matchList + i, i); 
    }

    // 2. Lookback Dynamic Programming
    for (size_t i = 1; i < matchNum; ++i) {
        uint32_t currentPos = matchList[i].qKmer.qInfo.pos;

        int bestParentIdx = -1;
        MatchScore bestScore;
        uint32_t bestHistoryMask = 0;

        for (int j = i - 1; j >= 0; --j) {
            uint32_t prevPos = matchList[j].qKmer.qInfo.pos;
            int shift = (currentPos - prevPos) / 3;
            
            if (shift > maxCodonShift) {
                break;
            }
            if (shift <= 0) {
                continue;
            }

            // The overlap check arguments swap order depending on the strand
            bool isOverlap = isForward ? 
                metamerPattern->checkOverlap(matchList[j].tKmer.value, matchList[i].tKmer.value, shift) :
                metamerPattern->checkOverlap(matchList[i].tKmer.value, matchList[j].tKmer.value, shift);

            if (isOverlap) {
                connectedToNext[j] = true;

                uint32_t shiftedHistoryMask = isForward ? 
                    ((localMatchPaths[j].historyMask << shift) & windowMask) : 
                    (localMatchPaths[j].historyMask >> shift);
                
                uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                
                MatchScore totalScore = localMatchPaths[j].score + 
                    metamerPattern->calMatchScore(
                        matchList[i].qKmer.value, 
                        matchList[i].tKmer.value, 
                        validPosMask,
                        *substitutionMatrix);

                if (totalScore.isLargerThan(bestScore, par.scoreMode)) {
                    bestScore = totalScore;
                    bestParentIdx = j;
                    bestHistoryMask = shiftedHistoryMask | metamerPattern->spaceMask;
                }
            }
        }

        // 3. Update the best path found for node i
        if (bestParentIdx != -1) {
            int shift = (currentPos - matchList[bestParentIdx].qKmer.qInfo.pos) / 3;

            localMatchPaths[i].start = localMatchPaths[bestParentIdx].start;                        
            localMatchPaths[i].score = bestScore;
            localMatchPaths[i].hammingDist = localMatchPaths[bestParentIdx].hammingDist + 
                metamerPattern->hammingDistSum(
                    matchList[i].qKmer.value, 
                    matchList[i].tKmer.value, 
                    shift,
                    isForward); 

            localMatchPaths[i].depth = localMatchPaths[bestParentIdx].depth + shift;
            localMatchPaths[i].startMatch = localMatchPaths[bestParentIdx].startMatch;
            localMatchPaths[i].historyMask = bestHistoryMask;
            localMatchPaths[i].prevMatchIdx = bestParentIdx;   
        }
    }

    // 4. Harvest Terminal Paths (Safely outside the main loop)
    for (size_t i = 0; i < matchNum; ++i) {
        if (!connectedToNext[i] && localMatchPaths[i].depth >= MIN_DEPTH) {
            int currIdx = i;
            while (currIdx != -1) {
                localMatchPaths[i].chain.push_back(matchList + currIdx);
                currIdx = localMatchPaths[currIdx].prevMatchIdx;
            }
            std::reverse(localMatchPaths[i].chain.begin(), localMatchPaths[i].chain.end());

            filteredMatchPaths.push_back(localMatchPaths[i]);
        }
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::ensureArraySize(size_t newSize) {
    if (newSize > arraySize_filterRedundantMatches) {
        delete[] bestMatchTaxIdForQuotient;
        delete[] minHammingForQuotient;
        bestMatchTaxIdForQuotient = new TaxID[newSize]();
        minHammingForQuotient = new uint8_t[newSize];
        arraySize_filterRedundantMatches = newSize;
    }
    std::memset(minHammingForQuotient, std::numeric_limits<uint8_t>::max(), newSize * sizeof(uint8_t));
}


template class Taxonomer<Match>;
template class Taxonomer<MatchWithPos>;