#include "IndexCreator.h"
#include <iostream>
#include <istream>
#include <string>
#include <vector>
#include "report.h"
#include "FileUtil.h"

using namespace std;

void setDefaults_filterByGenus(LocalParameters & par) {
    par.taxidCol = 3;
}

int filterByGenus(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_filterByGenus(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    string results = par.filenames[0];
    string genus_list = par.filenames[1];
    string taxonomy = par.filenames[2];

    // Load taxonomy
    NcbiTaxonomy ncbiTaxonomy(taxonomy + "/names.dmp",
                              taxonomy + "/nodes.dmp",
                              taxonomy + "/merged.dmp");

    // Load genus list
    vector<int> genusList;
    ifstream genusListFile(genus_list);
    if (!genusListFile.is_open()) {
        cerr << "Error: could not open genus list file " << genus_list << endl;
        return 1;
    }
    string line;
    while (getline(genusListFile, line)) {
        genusList.push_back(stoi(line));
    }

    vector<TaxID> observedSp;
    for(size_t i = 0; i < genusList.size(); i++){
        TaxID sp = ncbiTaxonomy.getTaxIdAtRank(genusList[i], "species");
        if(find(observedSp.begin(), observedSp.end(), sp) == observedSp.end()){
            observedSp.push_back(sp);
        }
    }
    cout << "Observed species: " << observedSp.size() << endl;

    // Read each result line
    ifstream resultsFile(results);
    if (!resultsFile.is_open()) {
        cerr << "Error: could not open results file " << results << endl;
        return 1;
    }

    vector<string> tokens;
    while (getline(resultsFile, line)) {
        // Parse line
        tokens.clear();
        tokens = Util::split(line, "\t");
        int taxid = stoi(tokens[par.taxidCol - 1]);

        if (taxid == 0) {
            continue;
        }

        int genus = ncbiTaxonomy.getTaxIdAtRank(taxid, "genus");
        if (genus == 0 || genus == 1) {
            cerr << "Error: could not find genus for taxid " << taxid << endl;
        }
        // Check if genus is in list
        bool found = false;
        for (int g: genusList) {
            if (g == genus) {
                found = true;
                break;
            }
        }
        if (!found) {
            continue;
        }
        // Print line
        cout << line << endl;
    }
}