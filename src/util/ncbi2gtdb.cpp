#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include <iostream>
#include "IndexCreator.h"
#include "LocalUtil.h"

using namespace std;

void setNcbi2gtdbDefault(LocalParameters & par){
    par.readIdCol = 1;
    par.taxidCol = 2;
}

int ncbi2gtdb(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    setNcbi2gtdbDefault(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string ncbi2gtdbFile = par.filenames[0];
    string taxonomy = par.filenames[1];

    // Load taxonomy
    NcbiTaxonomy ncbiTaxonomy(taxonomy + "/names.dmp",
                              taxonomy + "/nodes.dmp",
                              taxonomy + "/merged.dmp");

    vector<Query> binnings;
    unordered_map<TaxID, unsigned int> taxonCounts;
    unordered_map<TaxID, TaxID> ncbi2gtdbRes;
    unordered_map<TaxID, vector<TaxID>> ncbi2gtdbEntries;
    // Load old result
    ifstream ncbi2gtdb;
    ncbi2gtdb.open(ncbi2gtdbFile);
    string eachLine;
    vector<string> columns;
    string eachItem;
    int lineCount = 0;
    if (ncbi2gtdb.is_open()) {
        while (getline(ncbi2gtdb, eachLine)) {
            istringstream lineStream(eachLine);
            columns.clear();
            while (getline(lineStream, eachItem, '\t')) {
                columns.push_back(eachItem);
            }
            TaxID ncbiTaxid = stoi(columns[0]);
            TaxID gtdbTaxid = stoi(columns[1]);
            ncbi2gtdbEntries[ncbiTaxid].push_back(gtdbTaxid);
            // if (ncbi2gtdbRes.find(ncbiTaxid) == ncbi2gtdbRes.end()) {
            //     ncbi2gtdbRes[ncbiTaxid] = gtdbTaxid;
            // } else {
            //     TaxID lca = ncbiTaxonomy.LCA(ncbi2gtdbRes[ncbiTaxid], gtdbTaxid);
            //     ncbi2gtdbRes[ncbiTaxid] = lca;
            // }
            lineCount++;
        }
    } else {
        cerr << "Cannot open the file of binning result" << endl;
    }
    ncbi2gtdb.close();

    // Iterate ncbi2gtdbEntries
    int mapped2subspecies = 0;
    int mapped2species = 0;
    int mapped2genus = 0;
    int mapped2family = 0;
    int mapped2order = 0;
    int mapped2class = 0;
    int mapped2phylum = 0;
    int mapped2higher = 0;

    for (auto &item : ncbi2gtdbEntries) {
        TaxID ncbiTaxid = item.first;
        vector<TaxID> gtdbTaxids = item.second;
        const TaxonNode * lca = ncbiTaxonomy.LCA(gtdbTaxids);
        string gtdbRank = ncbiTaxonomy.getString(lca->rankIdx);
        if (gtdbRank == "subspecies") {
            mapped2subspecies += gtdbTaxids.size();
        } else if (gtdbRank == "species") {
            mapped2species += gtdbTaxids.size();
        } else if (gtdbRank == "genus") {
            mapped2genus += gtdbTaxids.size();
        } else if (gtdbRank == "family") {
            mapped2family += gtdbTaxids.size();
        } else if (gtdbRank == "order") {
            mapped2order += gtdbTaxids.size();
        } else if (gtdbRank == "class") {
            mapped2class += gtdbTaxids.size();
        } else if (gtdbRank == "phylum") {
            mapped2phylum += gtdbTaxids.size();
        } else {
            mapped2higher += gtdbTaxids.size();
        }
        ncbi2gtdbRes[ncbiTaxid] = lca->taxId;
    }
    
    // Iterate ncbi2gtdbRes
    int subspeciesCnt = 0;
    int speciesCnt = 0;
    int genusCnt = 0;
    int familyCnt = 0;
    int orderCnt = 0;
    int classCnt = 0;
    int phylumCnt = 0;
    int higherCnt = 0;

    for (auto &item : ncbi2gtdbRes) {
        TaxID gtdbTaxid = item.second;
        const TaxonNode *gtdbNode = ncbiTaxonomy.taxonNode(gtdbTaxid);
        string gtdbRank = ncbiTaxonomy.getString(gtdbNode->rankIdx);
        if (gtdbRank == "subspecies") {
            subspeciesCnt++;
        } else if (gtdbRank == "species") {
            speciesCnt++;
        } else if (gtdbRank == "genus") {
            genusCnt++;
        } else if (gtdbRank == "family") {
            familyCnt++;
        } else if (gtdbRank == "order") {
            orderCnt++;
        } else if (gtdbRank == "class") {
            classCnt++;
        } else if (gtdbRank == "phylum") {
            phylumCnt++;
        } else {
            higherCnt++;
        }
    }

    cout << "Processed " << lineCount << " lines" << endl;
    cout << "Subspecies: " << mapped2subspecies << "\t" << subspeciesCnt << endl;
    cout << "Species: " << mapped2species << "\t" << speciesCnt << endl;
    cout << "Genus: " << mapped2genus << "\t" << genusCnt << endl;
    cout << "Family: " << mapped2family << "\t" << familyCnt << endl;
    cout << "Order: " << mapped2order << "\t" << orderCnt << endl;
    cout << "Class: " << mapped2class << "\t" << classCnt << endl;
    cout << "Phylum: " << mapped2phylum << "\t" << phylumCnt << endl;
    cout << "Higher: " << mapped2higher << "\t" << higherCnt << endl;
    
  


    LocalUtil::writeMappingFile_text(ncbi2gtdbRes, ncbi2gtdbFile + ".processed");
    return 0;
}

