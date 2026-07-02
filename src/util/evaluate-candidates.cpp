#include "CandidateDBReader.h"
#include "Command.h"
#include "LocalParameters.h"
#include "common.h"

#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <unordered_map>
#include <vector>

using namespace std;

struct CandidateCountAtRank {
    int total = 0;
    int FP = 0;
    int TP = 0;
    int FN = 0;
    float precision = 0.0f;
    float sensitivity = 0.0f;
    float f1 = 0.0f;

    void calculate() {
        precision = (TP + FP) == 0 ? 0.0f : static_cast<float>(TP) / static_cast<float>(TP + FP);
        sensitivity = total == 0 ? 0.0f : static_cast<float>(TP) / static_cast<float>(total);
        f1 = (precision + sensitivity) == 0.0f ? 0.0f : 2.0f * precision * sensitivity / (precision + sensitivity);
    }
};

struct CandidateGradeResult {
    unordered_map<string, CandidateCountAtRank> countsAtRanks;
};

struct AlternativeCandidateHit {
    bool found = false;
    size_t candidateIdx = 0;
    TaxID topTaxId = 0;
    TaxID candidateTaxId = 0;
    float topScore = 0.0f;
    float candidateScore = 0.0f;
    float topEValue = 1.0f;
    float candidateEValue = 1.0f;
};

void setGradeDefault(LocalParameters &par);

namespace {

bool loadFileList(const string &listFileName, vector<string> &fileNames) {
    ifstream listFile(listFileName);
    if (!listFile.is_open()) {
        cerr << "Cannot open file list " << listFileName << endl;
        return false;
    }

    string line;
    while (getline(listFile, line)) {
        if (!line.empty()) {
            fileNames.push_back(line);
        }
    }
    return true;
}

bool normalizeReadId(
    string id,
    string &normalizedId,
    const string &testType)
{
    regex gtdbAccessionWithVersion("(GC[AF]_[0-9]+\\.[0-9]+)");
    regex gtdbAccessionNoVersion("(GC[AF]_[0-9]+)");
    smatch assacc;

    if (testType == "gtdb") {
        if (!regex_search(id, assacc, gtdbAccessionWithVersion)) {
            cerr << "Cannot parse GTDB read ID: " << id << endl;
            return false;
        }
        normalizedId = assacc[0];
        size_t pos = normalizedId.find('.');
        if (pos != string::npos) {
            normalizedId = normalizedId.substr(0, pos);
        }
    } else if (testType == "gtdb-amgsim") {
        if (!regex_search(id, assacc, gtdbAccessionNoVersion)) {
            cerr << "Cannot parse GTDB-aMGSIM read ID: " << id << endl;
            return false;
        }
        normalizedId = assacc[0];
    } else if (testType == "hiv" || testType == "hiv-ex") {
        size_t pos = id.find('_');
        normalizedId = id.substr(0, pos);
    } else if (testType == "cami" || testType == "cami-long" || testType == "cami-euk") {
        size_t pos = id.find('/');
        normalizedId = id.substr(0, pos);
    } else if (testType == "over") {
        if (!regex_search(id, assacc, gtdbAccessionWithVersion)) {
            cerr << "Cannot parse overclassification read ID: " << id << endl;
            return false;
        }
        normalizedId = assacc[0];
    } else if (testType == "kapk") {
        const size_t start = id.find("----");
        if (start == string::npos) {
            cerr << "Cannot parse KapK read ID: " << id << endl;
            return false;
        }
        const size_t accessionStart = start + 4;
        const size_t end = id.find("__", accessionStart);
        if (end == string::npos) {
            cerr << "Cannot parse KapK read ID: " << id << endl;
            return false;
        }
        normalizedId = id.substr(accessionStart, end - accessionStart);
    } else {
        normalizedId = id;
    }

    return true;
}

bool loadAnswerSheet(
    const string &mappingFileName,
    const string &testType,
    TaxonomyWrapper &taxonomy,
    unordered_map<string, TaxID> &assacc2taxid)
{
    ifstream map(mappingFileName);
    if (!map.is_open()) {
        cerr << "Cannot open file for answer: " << mappingFileName << endl;
        return false;
    }

    string key;
    string value;
    while (getline(map, key, '\t')) {
        getline(map, value, '\n');
        if (key.empty() || value.empty()) {
            continue;
        }
        if (testType != "kapk") {
            size_t pos = key.find('.');
            if (pos != string::npos) {
                key = key.substr(0, pos);
            }
        }

        TaxID taxId = static_cast<TaxID>(stoi(value));
        if (taxonomy.hasInternalTaxID()) {
            taxId = taxonomy.getInternalTaxID(taxId);
            if (taxId < 0) {
                continue;
            }
        }
        assacc2taxid[key] = taxId;
    }
    return true;
}

char probeTaxonAtRankCami(
    TaxID shot,
    TaxID target,
    const TaxonomyWrapper &taxonomy,
    const string &rank)
{
    if (rank == "subspecies") {
        if (shot == 1 || shot == 0) {
            return 'N';
        }

        const TaxonNode *shotNode = taxonomy.taxonNode(shot);
        if (strcmp(taxonomy.getString(shotNode->rankIdx), "no rank") != 0) {
            return 'N';
        }

        return shot == target ? 'O' : 'X';
    }

    TaxID targetTaxIdAtRank = taxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode *targetNode = taxonomy.taxonNode(targetTaxIdAtRank);
    int rankIdx = taxonomy.findRankIndex2(rank);
    if (taxonomy.findRankIndex2(taxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }

    if (shot == 1 || shot == 0) {
        return 'N';
    }

    TaxID shotTaxIdAtRank = taxonomy.getTaxIdAtRank(shot, rank);
    const TaxonNode *shotNode = taxonomy.taxonNode(shotTaxIdAtRank);
    if (taxonomy.findRankIndex2(taxonomy.getString(shotNode->rankIdx)) > rankIdx) {
        return 'N';
    }

    return shotTaxIdAtRank == targetTaxIdAtRank ? 'O' : 'X';
}

char probeTaxonAtRankCamiEuk(
    TaxID shot,
    TaxID target,
    TaxonomyWrapper &taxonomy,
    const string &rank)
{
    TaxID targetTaxIdAtRank = taxonomy.getTaxIdAtRank(target, rank);
    const TaxonNode *targetNode = taxonomy.taxonNode(targetTaxIdAtRank);
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(taxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }

    if (taxonomy.getTaxIdAtRank(target, "superkingdom") != 2759) {
        return '-';
    }

    if (shot == 1 || shot == 0) {
        return 'N';
    }

    TaxID shotTaxIdAtRank = taxonomy.getTaxIdAtRank(shot, rank);
    const TaxonNode *shotNode = taxonomy.taxonNode(shotTaxIdAtRank);
    if (NcbiTaxonomy::findRankIndex(taxonomy.getString(shotNode->rankIdx)) > rankIdx) {
        return 'N';
    }

    return shotTaxIdAtRank == targetTaxIdAtRank ? 'O' : 'X';
}

char probeTaxonOverclassification(
    TaxID shot,
    TaxID target,
    TaxonomyWrapper &taxonomy,
    const string &rank)
{
    const TaxonNode *targetNode = taxonomy.taxonNode(target);
    int rankIdx = NcbiTaxonomy::findRankIndex(rank);
    if (NcbiTaxonomy::findRankIndex(taxonomy.getString(targetNode->rankIdx)) > rankIdx) {
        return '-';
    }

    if (shot == 1 || shot == 0) {
        return 'N';
    }

    const TaxonNode *shotNode = taxonomy.taxonNode(shot);
    if (NcbiTaxonomy::findRankIndex(taxonomy.getString(shotNode->rankIdx)) > rankIdx) {
        return 'N';
    }

    return shot == target ? 'O' : 'X';
}

char probeTaxonHivExclusion(TaxID shot, TaxID target) {
    if (shot == 1 || shot == 0) {
        return 'N';
    }
    return shot == target ? 'O' : 'X';
}

char probeCandidateAtRank(
    TaxID candidateTaxId,
    TaxID target,
    TaxonomyWrapper &taxonomy,
    const LocalParameters &par,
    const string &rank)
{
    if (par.testType == "over") {
        return probeTaxonOverclassification(candidateTaxId, target, taxonomy, rank);
    } else if (par.testType == "hiv-ex") {
        return probeTaxonHivExclusion(candidateTaxId, 11676);
    } else if (par.testType == "cami-euk") {
        return probeTaxonAtRankCamiEuk(candidateTaxId, target, taxonomy, rank);
    } else {
        return probeTaxonAtRankCami(candidateTaxId, target, taxonomy, rank);
    }
}

char compareCandidateSetAtRank(
    const vector<SpeciesCandidate> &candidates,
    TaxID target,
    TaxonomyWrapper &taxonomy,
    const LocalParameters &par,
    const string &rank)
{
    char targetProbe = 'N';
    if (par.testType == "over") {
        targetProbe = probeTaxonOverclassification(0, target, taxonomy, rank);
    } else if (par.testType == "hiv-ex") {
        targetProbe = probeTaxonHivExclusion(0, 11676);
    } else if (par.testType == "cami-euk") {
        targetProbe = probeTaxonAtRankCamiEuk(0, target, taxonomy, rank);
    } else {
        targetProbe = probeTaxonAtRankCami(0, target, taxonomy, rank);
    }
    if (targetProbe == '-') {
        return '-';
    }

    bool sawFalseCandidate = false;
    for (const SpeciesCandidate &candidate : candidates) {
        char p = probeCandidateAtRank(candidate.speciesId, target, taxonomy, par, rank);
        if (p == 'O') {
            return 'O';
        }
        if (p == 'X') {
            sawFalseCandidate = true;
        }
    }

    return sawFalseCandidate ? 'X' : 'N';
}

AlternativeCandidateHit findAlternativeTruePositiveCandidate(
    const vector<SpeciesCandidate> &candidates,
    TaxID target,
    TaxonomyWrapper &taxonomy,
    const LocalParameters &par,
    const string &rank)
{
    AlternativeCandidateHit hit;
    if (candidates.size() < 2) {
        return hit;
    }

    const SpeciesCandidate &topCandidate = candidates.front();
    hit.topTaxId = topCandidate.speciesId;
    hit.topScore = topCandidate.idScore;
    hit.topEValue = std::exp(topCandidate.logE);

    if (probeCandidateAtRank(topCandidate.speciesId, target, taxonomy, par, rank) == 'O') {
        return hit;
    }

    for (size_t i = 1; i < candidates.size(); ++i) {
        if (probeCandidateAtRank(candidates[i].speciesId, target, taxonomy, par, rank) == 'O') {
            hit.found = true;
            hit.candidateIdx = i;
            hit.candidateTaxId = candidates[i].speciesId;
            hit.candidateScore = candidates[i].idScore;
            hit.candidateEValue = std::exp(candidates[i].logE);
            return hit;
        }
    }

    return hit;
}

TaxID printableTaxId(TaxID taxId, TaxonomyWrapper &taxonomy) {
    return taxonomy.hasInternalTaxID() ? taxonomy.getOriginalTaxID(taxId) : taxId;
}

void updateCount(char p, CandidateCountAtRank &count) {
    if (p == 'O') {
        count.TP++;
        count.total++;
    } else if (p == 'X') {
        count.FP++;
        count.total++;
    } else if (p == 'N') {
        count.FN++;
        count.total++;
    }
}

}

int evaluateCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setGradeDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string candidateDbFileList = par.filenames[0];
    const string mappingFileList = par.filenames[1];
    const string dbDir = par.filenames[2];

    vector<string> ranks;
    if (!par.testRank.empty()) {
        ranks = Util::split(par.testRank, ",");
    } else {
        ranks = {"class", "order", "family", "genus", "species"};
    }

    vector<string> mappingFileNames;
    vector<string> candidateDbFileNames;
    if (!loadFileList(mappingFileList, mappingFileNames)
        || !loadFileList(candidateDbFileList, candidateDbFileNames)) {
        return 1;
    }

    if (mappingFileNames.size() != candidateDbFileNames.size()) {
        cerr << "Candidate DB list and answer sheet list have different sizes." << endl;
        return 1;
    }

    size_t numberOfFiles = candidateDbFileNames.size();
    vector<CandidateGradeResult> results(numberOfFiles);

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

#pragma omp parallel default(none), shared(results, ranks, numberOfFiles, mappingFileNames, candidateDbFileNames, par, cout, cerr, dbDir)
    {
        vector<string> ranksLocal = ranks;
        unique_ptr<TaxonomyWrapper> taxonomy(loadTaxonomy(dbDir, par.taxonomyPath));

#pragma omp for schedule(dynamic)
        for (size_t i = 0; i < numberOfFiles; ++i) {
            unordered_map<string, TaxID> assacc2taxid;
            if (!loadAnswerSheet(mappingFileNames[i], par.testType, *taxonomy, assacc2taxid)) {
                continue;
            }

            CandidateDBReader candidateReader(candidateDbFileNames[i], 1);
            if (!candidateReader.open(DBReader<unsigned int>::SORT_BY_ID)) {
                cerr << "Cannot open candidate DB " << candidateDbFileNames[i] << endl;
                continue;
            }

            size_t totalReads = 0;
            size_t readsWithCandidates = 0;
            CandidateDBEntry entry;
            for (size_t j = 0; j < candidateReader.size(); ++j) {
                if (!candidateReader.getByIndex(j, entry, 0)) {
                    continue;
                }

                string normalizedId;
                if (!normalizeReadId(entry.queryName, normalizedId, par.testType)) {
                    continue;
                }

                auto answerIt = assacc2taxid.find(normalizedId);
                if (answerIt == assacc2taxid.end()) {
                    continue;
                }

                totalReads++;
                if (!entry.candidates.empty()) {
                    readsWithCandidates++;
                }

                for (const string &rank : ranksLocal) {
                    char p = compareCandidateSetAtRank(
                        entry.candidates,
                        answerIt->second,
                        *taxonomy,
                        par,
                        rank);
                    updateCount(p, results[i].countsAtRanks[rank]);
                    if (p == 'O') {
                        AlternativeCandidateHit altHit = findAlternativeTruePositiveCandidate(
                            entry.candidates,
                            answerIt->second,
                            *taxonomy,
                            par,
                            rank);
                        if (altHit.found) {
#pragma omp critical
                            {
                                cout << "ALT_TP\t"
                                     << candidateDbFileNames[i] << "\t"
                                     << entry.queryName << "\t"
                                     << normalizedId << "\t"
                                     << rank << "\t"
                                     << printableTaxId(answerIt->second, *taxonomy) << "\t"
                                     << printableTaxId(altHit.topTaxId, *taxonomy) << "\t"
                                     << altHit.topScore << "\t"
                                     << altHit.topEValue << "\t"
                                     << altHit.candidateIdx + 1 << "\t"
                                     << printableTaxId(altHit.candidateTaxId, *taxonomy) << "\t"
                                     << altHit.candidateScore << "\t"
                                     << altHit.candidateEValue << endl;
                            }
                        }
                    }
                    if (par.verbosity == 3) {
                        cout << entry.queryName << " " << normalizedId << " " << p << endl;
                    }
                }
            }
            candidateReader.close();

            for (const string &rank : ranksLocal) {
                results[i].countsAtRanks[rank].calculate();
            }

#pragma omp critical
            {
                cout << candidateDbFileNames[i] << endl;
                cout << "The number of reads: " << totalReads << endl;
                cout << "The number of reads with candidates: " << readsWithCandidates << endl;
                for (const string &rank : ranksLocal) {
                    CandidateCountAtRank &count = results[i].countsAtRanks[rank];
                    cout << rank << " " << count.total << " "
                         << count.TP + count.FP << " "
                         << count.TP << " " << count.FP << " "
                         << count.precision << " "
                         << count.sensitivity << " " << count.f1 << endl;
                }
                cout << endl;
            }
        }
    }

    cout << "Rank\t";
    for (size_t i = 0; i < results.size(); i++) {
        cout << "Precision\tSensitivity\tF1\t";
    }
    cout << endl;
    for (const string &rank : ranks) {
        cout << rank << "\t";
        for (auto &result : results) {
            CandidateCountAtRank &count = result.countsAtRanks[rank];
            cout << count.precision << "\t" << count.sensitivity << "\t" << count.f1 << "\t";
        }
        cout << endl;
    }

    return 0;
}
