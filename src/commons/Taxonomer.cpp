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
    // Parameters
    accessionLevel = par.accessionLevel;
    eukaryotaTaxId = taxonomy->getEukaryotaTaxID();


    dnaShift = 3;
    maxCodonShift = 1; // 1 by default when NO sycmer AND NO spacing
    if (par.syncmer && (kmerLen != windowSize)) {
        // Both syncmer AND spacing
        maxCodonShift = getFirstOneAfterFirstZero(metamerPattern->spaceMask) * 2;
        dnaShift = (kmerLen - par.smerLen) * 3;
        bitPerCodon = (static_cast<const SpacedPattern*>(metamerPattern))->bitPerCodon;
        bitPerAA = (static_cast<const SpacedPattern*>(metamerPattern))->bitPerAA;
    } else if (par.syncmer && (kmerLen == windowSize)) {
        // Only syncmer
        maxCodonShift = (kmerLen - par.smerLen);
        dnaShift = maxCodonShift * 3;
    } else if (kmerLen != windowSize) {
        // Only spacing
        maxCodonShift = getFirstOneAfterFirstZero(metamerPattern->spaceMask);
        bitPerCodon = (static_cast<const SpacedPattern*>(metamerPattern))->bitPerCodon;
        bitPerAA = (static_cast<const SpacedPattern*>(metamerPattern))->bitPerAA;
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
    sp2score.reserve(4096);
    tiedIndices.reserve(4096);
    tempPrioritySpecies.reserve(4096);

    // lowerRankClassification
    cladeCnt.reserve(4096);

    // filterRedundantMatches
    arraySize_filterRedundantMatches = 4096;
    bestMatchForQuotient = new const MatchType*[arraySize_filterRedundantMatches]();
    bestMatchTaxIdForQuotient = new TaxID[arraySize_filterRedundantMatches]();
    minHammingForQuotient = new uint8_t[arraySize_filterRedundantMatches]();
    touchedQuotients.reserve(4096);

    // Output
    taxCounts.reserve(4096);

    if (par.maxEValue > 0) {
        logMaxEValue = std::log(par.maxEValue);
        useEvalueFilter = true;
    } else {
        useEvalueFilter = false;
    }

    if (par.priorityTaxa.empty()) {
        priorityTaxa = {};
    } else {
        std::vector<std::string> taxaStr = Util::split(par.priorityTaxa, ",");
        for (const std::string &taxIdStr : taxaStr) {
            TaxID taxId = taxonomy->getInternalTaxID(stoi(taxIdStr));
            priorityTaxa.push_back(taxId);
        }
    }
}

template <typename MatchType>
Taxonomer<MatchType>::~Taxonomer() {
    delete[] bestMatchForQuotient;
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
        queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
        queryList[currentQuery].hammingDist = speciesScore.hammingDist;
        return;
    }

    // If there are two or more good species level candidates, find the LCA.
    if (speciesScore.LCA) {
        queryList[currentQuery].classification = speciesScore.taxId;
        queryList[currentQuery].idScore = speciesScore.score.idScore;
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


    if (par.printLog) cout << "#" << currentQuery << " " << taxonomy->getOriginalTaxID(speciesScore.taxId) << endl;
    if constexpr (std::is_same<MatchType, MatchWithPos>::value) {
        sp2totalReadLength[speciesScore.taxId] += queryList[currentQuery].queryLength + queryList[currentQuery].queryLength2;
        auto & speciesBins = sp2coverage[speciesScore.taxId];
        if (speciesBins.empty()) { speciesBins.resize(65536, 0); }
        uint8_t * bins = speciesBins.data();
        for (size_t i = speciesScore.matchPathRange.first; i < speciesScore.matchPathRange.second; i ++) {
            const vector<const MatchType*> & currentChain = combinedMatchPaths[i].chain;
            for (size_t j = 0; j < currentChain.size(); j ++) {
                const uint16_t binId = currentChain[j]->posId; 
                if (par.printLog) cout << binId << " ";
                if (bins[binId] < 255) {
                    bins[binId]++;
                }
            }
                if (par.printLog) cout << endl;
        }
    }

    // if constexpr (std::is_same_v<MatchType, MatchWithPos>) {
    //     auto & speciesBins = sp2coverage[speciesScore.taxId];
    //     if (speciesBins.empty()) { speciesBins.resize(65536, 0); }
    //     uint8_t * bins = speciesBins.data();
    //     uint16_t minPos = 65535;
    //     for (size_t i = speciesScore.matchPathRange.first; i < speciesScore.matchPathRange.second; i ++) {
    //         // const vector<const MatchType*> & currentChain = combinedMatchPaths[i].chain;
    //         uint16_t currStartPos = combinedMatchPaths[i].startMatch->posId;
    //         if (currStartPos < minPos) {
    //             minPos = currStartPos;
    //         }
    //         // for (size_t j = 0; j < currentChain.size(); j ++) {
    //         //     const uint16_t binId = currentChain[j]->posId; 
    //         //     if (bins[binId] < 255) {
    //         //         bins[binId]++;
    //         //     }
    //         // }
    //     }
    //     bins[minPos]++;

    // }


    // If score is not enough, classify to the parent of the selected species
    if (speciesScore.score.idScore < par.minSpScore) {
      queryList[currentQuery].classification = taxonomy->taxonNode(
              taxonomy->getTaxIdAtRank(speciesScore.taxId, "species"))->parentTaxId;
      queryList[currentQuery].idScore = speciesScore.score.idScore;
    //   queryList[currentQuery].subScore = speciesScore.score.subScore;
      queryList[currentQuery].eValue = std::exp(speciesScore.score.logE);
      queryList[currentQuery].hammingDist = speciesScore.hammingDist;
      return;
    }
    // Store classification results
    queryList[currentQuery].idScore = speciesScore.score.idScore;
    // queryList[currentQuery].subScore = speciesScore.score.subScore;
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
    touchedQuotients.clear();

    for (size_t i = bestSpeciesRange.first; i < bestSpeciesRange.second; i ++) {
        size_t currQuotient = matchList[i].qKmer.qInfo.pos / dnaShift;
        uint8_t hamming = metamerPattern->hammingDistSum(matchList[i].qKmer.value, matchList[i].tKmer.value);

        if (bestMatchForQuotient[currQuotient] == nullptr) {
            touchedQuotients.push_back(currQuotient);
            bestMatchForQuotient[currQuotient] = matchList + i;
            bestMatchTaxIdForQuotient[currQuotient] = matchList[i].tKmer.tInfo.taxId;
            minHammingForQuotient[currQuotient] = hamming;
        } else {
            if (hamming < minHammingForQuotient[currQuotient]) {
                bestMatchForQuotient[currQuotient] = matchList + i;
                bestMatchTaxIdForQuotient[currQuotient] = matchList[i].tKmer.tInfo.taxId;
                minHammingForQuotient[currQuotient] = hamming;
            } else if (hamming == minHammingForQuotient[currQuotient]) {
                bestMatchTaxIdForQuotient[currQuotient] = taxonomy->LCA(
                    bestMatchTaxIdForQuotient[currQuotient], matchList[i].tKmer.tInfo.taxId);
            }
        }
    }

    for (size_t quotient : touchedQuotients) {
        if (bestMatchTaxIdForQuotient[quotient] != 0) {
            taxCnt[bestMatchTaxIdForQuotient[quotient]]++;
        }
        bestMatchForQuotient[quotient] = nullptr;
        bestMatchTaxIdForQuotient[quotient] = 0;
        minHammingForQuotient[quotient] = std::numeric_limits<uint8_t>::max();
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
    int minSubSpeciesMatch = ((queryLength - 1)/denominator) + (kmerLen > 8) - (par.syncmer == 1);
    if (minSubSpeciesMatch < 0) {
        minSubSpeciesMatch = 0;
    }
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
    sp2score.clear();
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
        size_t speciesStart = i;
        size_t previousPathSize = matchPaths.size();
        // For current species
        while ((i < end + 1) && currentSpecies == matchList[i].tKmer.tInfo.speciesId) {
            uint8_t curFrame = matchList[i].qKmer.qInfo.frame;
            size_t frameStart = i;
            // For current frame
            while ((i < end + 1) && currentSpecies == matchList[i].tKmer.tInfo.speciesId && curFrame == matchList[i].qKmer.qInfo.frame) {
                i ++;
            }
            if (i > frameStart) {
                if (windowSize == kmerLen) {
                    getMatchPaths(matchList + frameStart, i - frameStart, matchPaths, currentSpecies);
                } else {
                    getSpacedMatchPaths(matchList + frameStart, i - frameStart, matchPaths, currentSpecies);
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
            // if (par.printLog==1) {   
            //     cout << "Combined score: " << score.idScore << " " << score.subScore << endl;
            // }
            score.idScore = score.idScore / queryLength;
            score.idScore = min(score.idScore, 1.0f);
            if ((score.idScore < par.minScore) || (useEvalueFilter && score.logE > logMaxEValue)) {
                continue;
            }

            sp2score.emplace_back(currentSpecies, score);
            if (score.idScore > 0.f) {
                meaningfulSpecies++;
            }
            if (score.isLargerThan(bestSpScore, par.scoreMode)) {
                bestSpScore = score;
                bestSpeciesRange = make_pair(speciesStart, i);
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
    tiedIndices.clear();

    const float diff = 0.09f;
    const float myTieRatio = (par.tieRatio - diff) + (bestSpScore.idScore * diff);
    // const float myTieRatio = (windowSize == kmerLen) 
    //                         ? par.tieRatio 
    //                         : (par.tieRatio - diff) + (bestSpScore.idScore * diff);
    
    for (size_t i = 0; i < sp2score.size(); i++) {
        if (sp2score[i].second.isLargerThan(bestSpScore * myTieRatio, par.scoreMode)) {
            maxSpecies.push_back(sp2score[i].first);
            tiedIndices.push_back(i);
            bestScore.score += sp2score[i].second;   
        }
    }

    if (maxSpecies.size() > 1 && !priorityTaxa.empty()) {
        tempPrioritySpecies.clear();
        TaxonScore tempBestScore; 
        
        for (size_t i = 0; i < maxSpecies.size(); ++i) {
            if (taxonomy->isAunderB(maxSpecies[i], priorityTaxa)) {
                tempPrioritySpecies.push_back(maxSpecies[i]);
                tempBestScore.score += sp2score[tiedIndices[i]].second; 
            }
        }
        
        if (!tempPrioritySpecies.empty()) {
            maxSpecies.assign(tempPrioritySpecies.begin(), tempPrioritySpecies.end());
            bestScore.score = tempBestScore.score;
        }
    }

    if (maxSpecies.size() > 1) {
        bestScore.taxId = taxonomy->LCA(maxSpecies)->taxId;
        bestScore.LCA = true;
        bestScore.score.idScore /= maxSpecies.size();
    } else if (maxSpecies.size() == 1) {
        bestScore.taxId = maxSpecies[0];
        bestScore.LCA = false;
        bestScore.matchPathRange = bestMatchPathRange;
    }
    
    bestScore.score.logE = bestSpScore.logE;
    return bestScore;                                 
}

template <typename MatchType>
void Taxonomer<MatchType>::sortMatchPath(std::vector<MatchPath<MatchType>> & matchPaths, size_t i) {
    if (par.scoreMode == 0) {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath<MatchType> &a, const MatchPath<MatchType> &b) {
           if (a.score.idScore != b.score.idScore) {
             return a.score.idScore > b.score.idScore;
           }
           if (a.coveredPosCnt != b.coveredPosCnt) {
             return a.coveredPosCnt < b.coveredPosCnt;
           }
           return a.start > b.start;
         });
    } else if (par.scoreMode == 1) {
        sort(matchPaths.begin() + i, matchPaths.end(),
         [](const MatchPath<MatchType> &a, const MatchPath<MatchType> &b) {
           if (a.score.subScore != b.score.subScore) {
             return a.score.subScore > b.score.subScore;
           }
           if (a.coveredPosCnt != b.coveredPosCnt) {
             return a.coveredPosCnt < b.coveredPosCnt;
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
           if (a.coveredPosCnt != b.coveredPosCnt) {
             return a.coveredPosCnt < b.coveredPosCnt;
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
            score += matchPaths[i].score;
            spanLength += matchPaths[i].end - matchPaths[i].start + 1;
            score.logE = computeLEM_logE(
                queryLength,
                matchPaths[i].end - matchPaths[i].start + 1,
                par.dbTotalLength,
                score.logP);
            combinedMatchPaths.push_back(std::move(matchPaths[i]));
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

                    if (windowSize == kmerLen) {
                        bool trimmed = trimMatchPath(matchPaths[i], combinedMatchPaths[j], overlappedLength);
                        if (!trimmed) {
                            isOverlapped = true;
                            break;
                        }
                        // if (overlappedLength < (windowSize * 3)) { 
                        //     trimMatchPath(matchPaths[i], combinedMatchPaths[j], overlappedLength);
                        // } else {
                        //     isOverlapped = true;
                        //     break;
                        // }
                    } else {
                        bool trimmed = trimSpacedMatchPath(matchPaths[i], combinedMatchPaths[j], overlappedLength);
                        if (!trimmed) {
                            isOverlapped = true;
                            break;
                        }
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

    return score;
}

template <typename MatchType>
bool Taxonomer<MatchType>::isMatchPathOverlapped(
    const MatchPath<MatchType> & matchPath1,
    const MatchPath<MatchType> & matchPath2) {
    return !((matchPath1.end < matchPath2.start) || (matchPath2.end < matchPath1.start));                                       
}


template <typename MatchType>
bool Taxonomer<MatchType>::trimSpacedMatchPath(
    MatchPath<MatchType> & newPath, 
    const MatchPath<MatchType> & existingPath, 
    int overlapLength) 
{
    int overlapAaNum = 0;
    if (overlapLength >= windowSize * 3) {
        overlapAaNum = windowSize;
    } else {
        overlapAaNum = overlapLength / 3;
    }
    int extraPanelty = overlapLength - overlapAaNum * 3;

    uint32_t baseMask = (1U << overlapAaNum) - 1;
    if (newPath.start < existingPath.start) { 
        // Trim right end of newPath
        if (newPath.rightEndTrimmed) return false; 
        newPath.rightEndTrimmed = true;

        newPath.end = existingPath.start - 1;
        uint32_t validPosMask = newPath.endMatch->qKmer.qInfo.frame < 3
            ? (newPath.lastHistoryMask & baseMask)
            : (newPath.lastHistoryMask & (baseMask << (windowSize - overlapAaNum)));
            
        newPath.score -= metamerPattern->calMatchScore(
            newPath.lastAAs,
            newPath.lastCodons_t,
            newPath.lastCodons_q,
            validPosMask);
        
        // newPath.score.idScore -= overlapLength % 3;
        newPath.score.idScore -= extraPanelty;
        
    } else {
        // Trim left end of newPath
        if (newPath.leftEndTrimmed) return false;
        newPath.leftEndTrimmed = true;
        
        newPath.start = existingPath.end + 1;
        uint32_t validPosMask = newPath.startMatch->qKmer.qInfo.frame < 3
            ? (newPath.firstHistoryMask & (baseMask << (windowSize - overlapAaNum)))
            : (newPath.firstHistoryMask & baseMask);

        newPath.score -= metamerPattern->calMatchScore(
            newPath.firstAAs,
            newPath.firstCodons_t,
            newPath.firstCodons_q,
            validPosMask);

        // newPath.score.idScore -= overlapLength % 3;
        newPath.score.idScore -= extraPanelty;
    }

    newPath.score.idScore = max(0.0f, newPath.score.idScore);
    
    return true;
}

template <typename MatchType>
bool Taxonomer<MatchType>::trimMatchPath(
    MatchPath<MatchType> & newPath, 
    const MatchPath<MatchType> & existingPath, 
    int overlapLength) 
{
    int overlapAaNum = 0;
    if (overlapLength >= windowSize * 3) {
        overlapAaNum = windowSize;
    } else {
        overlapAaNum = overlapLength / 3;
    }
    int extraPanelty = overlapLength - overlapAaNum * 3;

    bool trimFromRight = (newPath.endMatch->qKmer.qInfo.frame < 3) ^ (newPath.start >= existingPath.start);

    if (newPath.start < existingPath.start) { 

        if (newPath.rightEndTrimmed) return false;
        newPath.rightEndTrimmed = true;

        newPath.end = existingPath.start - 1;
        newPath.score -= metamerPattern->calMatchScore(
            newPath.endMatch->qKmer.value,
            newPath.endMatch->tKmer.value,
            overlapAaNum,
            trimFromRight);            
        
    } else {
        if (newPath.leftEndTrimmed) return false;
        newPath.leftEndTrimmed = true;

        newPath.start = existingPath.end + 1;
        newPath.score -= metamerPattern->calMatchScore(
            newPath.startMatch->qKmer.value,
            newPath.startMatch->tKmer.value,
            overlapAaNum,
            trimFromRight);
    }
    newPath.score.idScore -= extraPanelty;
    newPath.score.idScore = max(0.0f, newPath.score.idScore);

    return true;
}

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
            true),
        kmerLen,
        metamerPattern->windowSize * 3);
}

template <typename MatchType>
void Taxonomer<MatchType>::addTerminalMatchPath(
    size_t pathIdx,
    const MatchType * matchList,
    vector<MatchPath<MatchType>> & filteredMatchPaths)
{
    MatchPath<MatchType> &path = localMatchPaths[pathIdx];
    if (par.storeKmerPos) {
        path.chain.clear();
        uint64_t currIdx = pathIdx;
        while (currIdx != UINT64_MAX) {
            path.chain.push_back(matchList + currIdx);
            currIdx = localMatchPaths[currIdx].prevMatchIdx;
        }
        std::reverse(path.chain.begin(), path.chain.end());
    }
    filteredMatchPaths.push_back(std::move(path));
}

template <typename MatchType>
void Taxonomer<MatchType>::getMatchPaths_lookbackDP(
    const MatchType * matchList,
    size_t matchNum,
    vector<MatchPath<MatchType>> & filteredMatchPaths,
    TaxID speciesId) 
{
    bool isForward = matchList[0].qKmer.qInfo.frame < 3;
    int MIN_COVERED_POS = taxonomy->IsAncestor(eukaryotaTaxId, speciesId) ? 
                          (int) par.minAaMatchEuk : (int) par.minAaMatch;

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
            localMatchPaths[i].coveredPosCnt = localMatchPaths[bestParentIdx].coveredPosCnt + shift;
            localMatchPaths[i].startMatch = localMatchPaths[bestParentIdx].startMatch;
            localMatchPaths[i].prevMatchIdx = bestParentIdx;   
        }
    }

    // 4. Harvest Terminal Paths
    for (size_t i = 0; i < matchNum; ++i) {
        if (!connectedToNext[i] && localMatchPaths[i].coveredPosCnt >= MIN_COVERED_POS) {
            addTerminalMatchPath(i, matchList, filteredMatchPaths);
        }
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::getMatchPaths(
    const MatchType * matchList,
    size_t matchNum,
    vector<MatchPath<MatchType>> & filteredMatchPaths,
    TaxID speciesId) 
{
    size_t i = 0;
    size_t currPos = matchList[0].qKmer.qInfo.pos;  
    uint64_t frame = matchList[0].qKmer.qInfo.frame;
    int MIN_COVERED_POS = (int) par.minAaMatch;
    if (taxonomy->IsAncestor(eukaryotaTaxId, speciesId)) {
        MIN_COVERED_POS = (int) par.minAaMatchEuk;
    }
        
    connectedToNext.resize(matchNum);
    fill(connectedToNext.begin(), connectedToNext.end(), false);
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    // Track frame direction once to eliminate the massive duplicate if/else branches
    bool isForward = (frame < 3);

    size_t curPosMatchStart = i;        
    while (i < matchNum && matchList[i].qKmer.qInfo.pos == currPos) {
        localMatchPaths[i] = makeMatchPath(matchList + i);
        ++ i;
    }
    size_t curPosMatchEnd = i; // exclusive

    while (i < matchNum) {
        uint32_t nextPos = matchList[i].qKmer.qInfo.pos;
        size_t nextPosMatchStart = i;
        while (i < matchNum  && nextPos == matchList[i].qKmer.qInfo.pos) {
            localMatchPaths[i] = makeMatchPath(matchList + i);
            ++ i;
        }
        size_t nextPosMatchEnd = i; // exclusive

        // Check if current position and next position are consecutive
        int shift = (nextPos - currPos) / 3;
        if (shift > 0 && shift <= maxCodonShift) {
            // Compare curPosMatches and nextPosMatches.
            for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                const MatchPath<MatchType> * bestPath = nullptr;
                uint64_t bestParentIdx = UINT64_MAX;
                MatchScore bestScore;
                for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
                    
                    // Conditionally swap arguments based on frame direction
                    bool isOverlap = isForward 
                        ? metamerPattern->checkOverlap(matchList[curIdx].tKmer.value, matchList[nextIdx].tKmer.value, shift)
                        : metamerPattern->checkOverlap(matchList[nextIdx].tKmer.value, matchList[curIdx].tKmer.value, shift);

                    if (isOverlap) {
                        connectedToNext[curIdx] = true;
                        if (localMatchPaths[curIdx].score.isLargerThan(bestScore, par.scoreMode)) {
                            bestPath = &localMatchPaths[curIdx];
                            bestScore = localMatchPaths[curIdx].score;
                            bestParentIdx = curIdx;
                        }
                    }
                }
                if (bestPath != nullptr) {
                    localMatchPaths[nextIdx].start = bestPath->start;                        
                    localMatchPaths[nextIdx].score = bestPath->score +
                        metamerPattern->calMatchScore(
                            matchList[nextIdx].qKmer.value, 
                            matchList[nextIdx].tKmer.value, 
                            shift,
                            isForward); // Pass the boolean flag dynamically here

                    localMatchPaths[nextIdx].coveredPosCnt = bestPath->coveredPosCnt + shift;
                    localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                    localMatchPaths[nextIdx].prevMatchIdx = bestParentIdx;
                }
            }
        } 
        
        for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
            if (!connectedToNext[curIdx] && localMatchPaths[curIdx].coveredPosCnt >= MIN_COVERED_POS) {
                addTerminalMatchPath(curIdx, matchList, filteredMatchPaths);
            }
        }

        curPosMatchStart = nextPosMatchStart;
        curPosMatchEnd = nextPosMatchEnd;
        currPos = nextPos;    
    }
    
    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
        if (localMatchPaths[curIdx].coveredPosCnt >= MIN_COVERED_POS) {
            addTerminalMatchPath(curIdx, matchList, filteredMatchPaths);
        }
    }
}


template <typename MatchType>
void Taxonomer<MatchType>::makeSpacedMatchPath(
    const MatchType * match,
    size_t index)
{
    localMatchPaths[index] = MatchPath<MatchType>(match, kmerLen, metamerPattern->windowSize * 3);
    localMatchPaths[index].score = metamerPattern->calMatchScore(
                                        match->qKmer.value, 
                                        match->tKmer.value);
 
    localMatchPaths[index].lastHistoryMask  = metamerPattern->spaceMask;
    localMatchPaths[index].lastAAs = disperseBits(
        match->tKmer.value >> metamerPattern->totalDNABits,
        metamerPattern->spaceMask,
        bitPerAA);
    localMatchPaths[index].lastCodons_t = disperseBits(
        match->tKmer.value & metamerPattern->dnaMask,
        metamerPattern->spaceMask,
        bitPerCodon);
    localMatchPaths[index].lastCodons_q = disperseBits(
        match->qKmer.value & metamerPattern->dnaMask,
        metamerPattern->spaceMask,
        bitPerCodon);

                                                
    localMatchPaths[index].firstHistoryMask = localMatchPaths[index].lastHistoryMask;
    localMatchPaths[index].firstAAs         = localMatchPaths[index].lastAAs;
    localMatchPaths[index].firstCodons_t    = localMatchPaths[index].lastCodons_t;
    localMatchPaths[index].firstCodons_q    = localMatchPaths[index].lastCodons_q;
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
    
    int MIN_COVERED_POS = taxonomy->IsAncestor(eukaryotaTaxId, speciesId) ? 
                          (int) par.minAaMatchEuk : (int) par.minAaMatch;

    connectedToNext.assign(matchNum, false);
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    // 1. Initialize all match paths sequentially
    for (size_t i = 0; i < matchNum; ++i) {
        makeSpacedMatchPath(matchList + i, i); 
    }

    // 2. Lookback Dynamic Programming
    for (size_t i = 1; i < matchNum; ++i) {
        uint32_t currentPos = matchList[i].qKmer.qInfo.pos;

        int bestParentIdx = -1;
        MatchScore bestScore;
        uint32_t bestHistoryMask = 0;
        int newCoveredPosCnt = 0;

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
                    ((localMatchPaths[j].lastHistoryMask << shift) & windowMask) : 
                    (localMatchPaths[j].lastHistoryMask >> shift);
                
                uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                
                MatchScore totalScore = localMatchPaths[j].score + 
                    metamerPattern->calMatchScore(
                        matchList[i].qKmer.value, 
                        matchList[i].tKmer.value, 
                        validPosMask);

                if (totalScore.isLargerThan(bestScore, par.scoreMode)) {
                    bestScore = totalScore;
                    bestParentIdx = j;
                    bestHistoryMask = shiftedHistoryMask | metamerPattern->spaceMask;
                    newCoveredPosCnt = __builtin_popcount(static_cast<unsigned int>(validPosMask));
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

            localMatchPaths[i].coveredPosCnt = localMatchPaths[bestParentIdx].coveredPosCnt + newCoveredPosCnt;
            localMatchPaths[i].startMatch = localMatchPaths[bestParentIdx].startMatch;
            localMatchPaths[i].lastHistoryMask = bestHistoryMask;
            localMatchPaths[i].prevMatchIdx = bestParentIdx;   
        }
    }

    // 4. Harvest Terminal Paths (Safely outside the main loop)
    for (size_t i = 0; i < matchNum; ++i) {
        if (!connectedToNext[i] && localMatchPaths[i].coveredPosCnt >= MIN_COVERED_POS) {
            addTerminalMatchPath(i, matchList, filteredMatchPaths);
        }
    }
}

template <typename MatchType>
void Taxonomer<MatchType>::getSpacedMatchPaths(
    const MatchType * matchList,
    size_t matchNum,
    vector<MatchPath<MatchType>> & filteredMatchPaths,
    TaxID speciesId) 
{
    if (matchNum == 0) return;

    size_t i = 0;
    size_t currPos = matchList[0].qKmer.qInfo.pos;  
    uint64_t frame = matchList[0].qKmer.qInfo.frame;
    bool isForward = (frame < 3);

    int MIN_COVERED_POS = taxonomy->IsAncestor(eukaryotaTaxId, speciesId) ? 
                          (int) par.minAaMatchEuk : (int) par.minAaMatch;
        
    connectedToNext.assign(matchNum, false); 
    localMatchPaths.clear();
    localMatchPaths.resize(matchNum);

    size_t curPosMatchStart = i;        
    while (i < matchNum && matchList[i].qKmer.qInfo.pos == currPos) {
        makeSpacedMatchPath(matchList + i, i);
        ++i;
    }
    size_t curPosMatchEnd = i; // exclusive

    while (i < matchNum) {
        uint32_t nextPos = matchList[i].qKmer.qInfo.pos;
        size_t nextPosMatchStart = i;
        while (i < matchNum  && nextPos == matchList[i].qKmer.qInfo.pos) {
            makeSpacedMatchPath(matchList + i, i);
            ++i;
        }
        size_t nextPosMatchEnd = i; // exclusive

        // Check if current position and next position are consecutive
        int shift = (nextPos - currPos) / 3;
        if (shift > 0 && shift <= maxCodonShift) {
            for (size_t nextIdx = nextPosMatchStart; nextIdx < nextPosMatchEnd; nextIdx++) {
                const MatchPath<MatchType> * bestPath = nullptr;
                MatchScore bestScore;
                uint32_t bestHistoryMask = 0;
                uint64_t bestLastCodons = 0;
                uint64_t bestParentIdx = UINT64_MAX;
                int newCoveredPosCnt = 0;
                
                for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {

                    uint32_t shiftedHistoryMask = isForward ? 
                        ((localMatchPaths[curIdx].lastHistoryMask << shift) & windowMask) : 
                        (localMatchPaths[curIdx].lastHistoryMask >> shift);
                    
                    uint32_t checkPos = shiftedHistoryMask & metamerPattern->spaceMask;
                    uint64_t checkPosExtended = stretchBits(checkPos, bitPerCodon);
                    
                    uint64_t shiftedLastCodons = isForward ? 
                        (localMatchPaths[curIdx].lastCodons_t << (shift * bitPerCodon)) : 
                        (localMatchPaths[curIdx].lastCodons_t >> (shift * bitPerCodon));
                    
                    uint64_t shiftedLastCodons_masked = shiftedLastCodons & checkPosExtended;
                    uint64_t currentCodons_masked = localMatchPaths[nextIdx].lastCodons_t & checkPosExtended;
                    
                    if (shiftedLastCodons_masked == currentCodons_masked) {
                        connectedToNext[curIdx] = true;
                        const uint32_t validPosMask = metamerPattern->spaceMask & (~shiftedHistoryMask);
                        const MatchScore totalScore = localMatchPaths[curIdx].score +
                            metamerPattern->calMatchScore(
                                matchList[nextIdx].qKmer.value, 
                                matchList[nextIdx].tKmer.value, 
                                validPosMask);

                        if (totalScore.isLargerThan(bestScore, par.scoreMode)) {
                            bestParentIdx = curIdx;
                            bestPath = &localMatchPaths[curIdx];
                            bestScore = totalScore;
                            bestLastCodons = shiftedLastCodons | localMatchPaths[nextIdx].lastCodons_t;
                            bestHistoryMask = shiftedHistoryMask | metamerPattern->spaceMask;
                            newCoveredPosCnt = __builtin_popcount(static_cast<unsigned int>(validPosMask));
                        }
                    }
                }
                
                if (bestPath != nullptr) {
                    // Store left end of the path
                    int depth = (localMatchPaths[nextIdx].start - bestPath->start) / 3;
                    if (depth < windowSize) {
                        uint32_t shiftedFirstHistroyMask = isForward ? 
                            (localMatchPaths[nextIdx].firstHistoryMask >> depth) :
                            (localMatchPaths[nextIdx].firstHistoryMask << depth);
                        localMatchPaths[nextIdx].firstHistoryMask = shiftedFirstHistroyMask | bestPath->firstHistoryMask;

                        uint64_t shiftedFirstCodons_t = isForward ? 
                            (localMatchPaths[nextIdx].firstCodons_t >> (depth * bitPerCodon)) :
                            (localMatchPaths[nextIdx].firstCodons_t << (depth * bitPerCodon));
                        localMatchPaths[nextIdx].firstCodons_t = shiftedFirstCodons_t | bestPath->firstCodons_t;

                        uint64_t shiftedFirstCodons_q = isForward ? 
                            (localMatchPaths[nextIdx].firstCodons_q >> (depth * bitPerCodon)) :
                            (localMatchPaths[nextIdx].firstCodons_q << (depth * bitPerCodon));
                        localMatchPaths[nextIdx].firstCodons_q = shiftedFirstCodons_q | bestPath->firstCodons_q;

                        uint64_t shiftedFirstAAs = isForward ? 
                            (localMatchPaths[nextIdx].firstAAs >> (depth * bitPerAA)) :
                            (localMatchPaths[nextIdx].firstAAs << (depth * bitPerAA));
                        localMatchPaths[nextIdx].firstAAs = shiftedFirstAAs | bestPath->firstAAs;
                    }
                    localMatchPaths[nextIdx].start = bestPath->start;  
                    localMatchPaths[nextIdx].startMatch = bestPath->startMatch;
                    
                    // Store right end of the path
                    localMatchPaths[nextIdx].lastHistoryMask = bestHistoryMask;
                    localMatchPaths[nextIdx].lastCodons_t = bestLastCodons;

                    uint64_t shiftedLastCodons_q = isForward ? 
                        (bestPath->lastCodons_q << (shift * bitPerCodon)) :
                        (bestPath->lastCodons_q >> (shift * bitPerCodon));
                    localMatchPaths[nextIdx].lastCodons_q = shiftedLastCodons_q | localMatchPaths[nextIdx].lastCodons_q;

                    uint64_t shiftedLastAAs = isForward ? 
                        (bestPath->lastAAs << (shift * bitPerAA)) :
                        (bestPath->lastAAs >> (shift * bitPerAA));
                    localMatchPaths[nextIdx].lastAAs = shiftedLastAAs | localMatchPaths[nextIdx].lastAAs;

                    // Update scores
                    localMatchPaths[nextIdx].score = bestScore;
                    localMatchPaths[nextIdx].coveredPosCnt = bestPath->coveredPosCnt + newCoveredPosCnt; 
                    localMatchPaths[nextIdx].prevMatchIdx = bestParentIdx;                   
                }
            }
        } 
        
        for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
            if (!connectedToNext[curIdx] && localMatchPaths[curIdx].coveredPosCnt >= MIN_COVERED_POS) {
                float gapPenalty = 
                    (localMatchPaths[curIdx].end - localMatchPaths[curIdx].start + 1) // spanned bases
                    - localMatchPaths[curIdx].coveredPosCnt * 3; // AA-matched bases
                localMatchPaths[curIdx].score.idScore = max(0.0f, localMatchPaths[curIdx].score.idScore - gapPenalty);
                
                addTerminalMatchPath(curIdx, matchList, filteredMatchPaths);
            }
        }
        
        curPosMatchStart = nextPosMatchStart;
        curPosMatchEnd = nextPosMatchEnd;
        currPos = nextPos;    
    }

    for (size_t curIdx = curPosMatchStart; curIdx < curPosMatchEnd; ++curIdx) {
        if (localMatchPaths[curIdx].coveredPosCnt >= MIN_COVERED_POS) {
            float gapPenalty = 
                (localMatchPaths[curIdx].end - localMatchPaths[curIdx].start + 1) // spanned bases
                - localMatchPaths[curIdx].coveredPosCnt * 3; // AA-matched bases
            localMatchPaths[curIdx].score.idScore = max(0.0f, localMatchPaths[curIdx].score.idScore - gapPenalty);
            
            addTerminalMatchPath(curIdx, matchList, filteredMatchPaths);
        }
    }    
}

template <typename MatchType>
void Taxonomer<MatchType>::ensureArraySize(size_t newSize) {
    if (newSize > arraySize_filterRedundantMatches) {
        delete[] bestMatchForQuotient;
        delete[] bestMatchTaxIdForQuotient;
        delete[] minHammingForQuotient;
        bestMatchForQuotient = new const MatchType*[newSize]();
        bestMatchTaxIdForQuotient = new TaxID[newSize]();
        minHammingForQuotient = new uint8_t[newSize];
        arraySize_filterRedundantMatches = newSize;
    }
}

template class Taxonomer<Match>;
template class Taxonomer<MatchWithPos>;
