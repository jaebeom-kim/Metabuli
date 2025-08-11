#include "NcbiTaxonomy.h"
#include "LocalParameters.h"
#include <Command.h>
#include <string>
#include "KSeqWrapper.h"
#include <iostream>
#include "IndexCreator.h"
#include "FileUtil.h"
#include "common.h"
#include <regex>

using namespace std;

void setDefaults_addToLibrary(LocalParameters & par){
    par.taxonomyPath = "DBDIR/taxonomy/" ;
    par.libraryPath = "DBDIR/library/";
}

// Group sequences by species
int addToLibrary(int argc, const char **argv, const Command &command){
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_addToLibrary(par);
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);

    const string fileList = par.filenames[0];
    const string mappingFileName = par.filenames[1];
    const string dbDir = par.filenames[2];
    if (par.taxonomyPath == "DBDIR/taxonomy/") par.taxonomyPath = dbDir + "/taxonomy/";
    if (par.libraryPath == "DBDIR/library/") par.libraryPath = dbDir + "/library/";

    if (!FileUtil::directoryExists(par.libraryPath.c_str())) {
        FileUtil::makeDir(par.libraryPath.c_str());
    }
    
    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);

    unordered_map<string, TaxID> accession2taxid;
    vector<string> fileNames;
    getObservedAccessionList(fileList, fileNames, accession2taxid);
    fillAcc2TaxIdMap(accession2taxid, mappingFileName);

    if(!par.assembly) {
        vector<string> unmapped;
        for (size_t i = 0; i < fileNames.size(); ++i) {
            KSeqWrapper* kseq = KSeqFactory(fileNames[i].c_str());
            while (kseq->ReadEntry()) {
                const KSeqWrapper::KSeqEntry & e = kseq->entry;
                char* pos = strchr(e.name.s, '.'); 
                if (pos != nullptr) {
                    *pos = '\0';
                }

                TaxID taxId = accession2taxid[e.name.s];
                if (taxId == 0) {
                    cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
                         " is not found in the mapping file. It is skipped.\n";
                    unmapped.push_back(e.name.s);
                    continue;
                }

                int speciesTaxID = taxonomy->getTaxIdAtRank(taxId, "species");
                if (speciesTaxID == 0) {
                    cout << "During processing " << fileNames[i] << ", accession " << e.name.s <<
                         " is not matched to any species. It is skipped.\n";
                    unmapped.push_back(e.name.s);
                    continue;
                }

                // Write each sequence to file with species taxID as file name
                FILE *file = fopen((dbDir + "/library/" + to_string(speciesTaxID) + ".fna").c_str(), "a");
                fprintf(file, ">%s %s\n", e.name.s, e.comment.s);
                fprintf(file, "%s\n", e.sequence.s);
                fclose(file);
            }
            delete kseq;
        }

        // Write unmapped accession to file
        FILE *file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
        for (const auto & i : unmapped) {
            fprintf(file, "%s\n", i.c_str());
        }
        fclose(file);
    }
    else // ASSEMBLY
    {
        // Load mapping file
        cout << "Load mapping from assembly accession to taxonomy ID" << endl;
        unordered_map<string, int> assembly2taxid;
        unordered_map<string, int > acc2taxid;
        if (FILE *mappingFile = fopen(mappingFileName.c_str(), "r")) {
            char buffer[512];
            int taxID;
            while (fscanf(mappingFile, "%s\t%d", buffer, &taxID) == 2) {
                // Remove the version number
                string accession = string(buffer);
                size_t pos = accession.find('.');
                if (pos != string::npos) { accession = accession.substr(0, pos); }
                assembly2taxid[accession] = taxID;
            }
        } else {
            cerr << "Cannot open the mapping from assembly accession to tax ID" << endl;
        }

        // Process each file
        vector<string> unmapped;
        regex regex1("(GC[AF]_[0-9]+\\.[0-9]+)");

        for (size_t i = 0; i < fileNames.size(); ++i) {
            // Get assembly accession from file name using regex and remove the version number
            smatch match;
            regex_search(fileNames[i], match, regex1);
            string assemblyID = match[0];
            size_t pos = assemblyID.find('.');
            if (pos != string::npos) { assemblyID = assemblyID.substr(0, pos); }

            // Skip if current assembly accession is not in the mapping file
            if (assembly2taxid.find(assemblyID) == assembly2taxid.end()) {
                cout << "During processing " << fileNames[i] << ", accession " << assemblyID <<
                     " is not found in the mapping file. It is skipped.\n";
                unmapped.push_back(assemblyID);
                continue;
            }

            // Get species taxID
            int speciesTaxID = taxonomy->getTaxIdAtRank(assembly2taxid[assemblyID], "species");
            if (speciesTaxID == 0) {
                cout << "During processing " << fileNames[i] << ", accession " << assemblyID <<
                     " is not matched to any species. It is skipped.\n";
                continue;
            }

            KSeqWrapper* kseq = KSeqFactory(fileNames[i].c_str());
            while (kseq->ReadEntry()){
                const KSeqWrapper::KSeqEntry & e = kseq->entry;
                // Extract accession
                string accession = string(e.name.s);
                acc2taxid[accession] = assembly2taxid[assemblyID];

                // Write to file
                FILE *file = fopen((dbDir + "/library/" + to_string(speciesTaxID) + ".fna").c_str(), "a");
                fprintf(file, ">%s %s\n", e.name.s, e.comment.s);
                fprintf(file, "%s\n", e.sequence.s);
                fclose(file);
            }
            delete kseq;
        }

        // Write unmapped accession to file
        FILE *file = fopen((dbDir + "/unmapped.txt").c_str(), "w");
        for (const auto & i : unmapped) {
            fprintf(file, "%s\n", i.c_str());
        }
        fclose(file);

        // Write mapping file
        cout << "Write mapping from accession to taxonomy ID" << endl;
        file = fopen((dbDir + "/my.accession2taxid").c_str(), "w");
        fprintf(file, "accession\taccession.version\ttaxid\tgi");
        for (auto & it : acc2taxid) {
            // Get accession without a version number
            string accession = it.first;
            size_t pos = accession.find('.');
            if (pos != string::npos) { accession = accession.substr(0, pos);}
            fprintf(file, "\n%s\t%s\t%d\t0", accession.c_str(), it.first.c_str(), it.second);
        }
        fclose(file);
    }
    delete taxonomy;
    return EXIT_SUCCESS;
}