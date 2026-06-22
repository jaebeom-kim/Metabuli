#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include "Reporter.h"

void setExtractDefaults(LocalParameters & par){
    par.taxonomyPath = "" ;
    par.outputDir = "";
    par.seqMode = 2;
    par.targetTaxId = 0;
    par.excludeTaxid = "";
    par.extractMode = 0;
}

void extractBaseNameAndExtension(const string &queryFileName, string & dirPath, string &baseName, string &extension) {
    size_t lastDotPos = queryFileName.find_last_of('.');
    if (lastDotPos == string::npos) {
        // No extension
        baseName = queryFileName;
        extension = "";
    } else if (queryFileName.substr(lastDotPos) == ".gz") {
        // Compressed file
        size_t secondLastDotPos = queryFileName.substr(0, lastDotPos).find_last_of('.');
        if (secondLastDotPos == string::npos) {
            baseName = queryFileName.substr(0, lastDotPos);
            extension = "";
        } else {
            baseName = queryFileName.substr(0, secondLastDotPos);
            extension = queryFileName.substr(secondLastDotPos, lastDotPos - secondLastDotPos);
        }
    } else {
        // Regular file with an extension
        baseName = queryFileName.substr(0, lastDotPos);
        extension = queryFileName.substr(lastDotPos);
    }
    size_t lastSlashPos = baseName.find_last_of('/');
    if (lastSlashPos != string::npos) {
        dirPath = baseName.substr(0, lastSlashPos);
        baseName = baseName.substr(lastSlashPos + 1);
    } else {
        dirPath = "";
    }

}

int extract(int argc, const char **argv, const Command& command) {
    LocalParameters & par = LocalParameters::getLocalInstance();
    setExtractDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.seqMode == 2) {
        // Check if the second argument is a directory
        if (FileUtil::directoryExists(par.filenames[2].c_str())) {
            cerr << "For '--seq-mode 2', please provide two query files." << endl;
            exit(1);
        }
    } else {
        // Check if the second argument is file
        if (!FileUtil::directoryExists(par.filenames[2].c_str())) {
            cerr << "For '--seq-mode 1' and '--seq-mode 3', please provide one query file." << endl;
            exit(1);
        }
    }

    const bool excludeMode = !par.excludeTaxid.empty();
    if (par.targetTaxId == 0 && !excludeMode) {
        cerr << "Please provide a target taxon ID with --tax-id or --exclude-taxid parameter." << endl;
        exit(1);
    }
    if (par.targetTaxId != 0 && excludeMode) {
        cerr << "Please provide only one of --tax-id or --exclude-taxid." << endl;
        exit(1);
    }

    if (!par.outputDir.empty()) {
        if (!FileUtil::directoryExists(par.outputDir.c_str())) {
            FileUtil::makeDir(par.outputDir.c_str());
        }
    }

    

    string classificationFileName = par.filenames[1 + (par.seqMode == 2)];
    string dbDir = par.filenames[2 + (par.seqMode == 2)];
    TaxID externalTaxID = par.targetTaxId;
    if (excludeMode) {
        if (par.excludeTaxid.find(',') != string::npos) {
            cerr << "Please provide a single target taxon ID with --exclude-taxid for extract." << endl;
            exit(1);
        }
        try {
            externalTaxID = stoi(par.excludeTaxid);
        } catch (...) {
            cerr << "Invalid taxon ID provided with --exclude-taxid: " << par.excludeTaxid << endl;
            exit(1);
        }
        if (externalTaxID == 0) {
            cerr << "Please provide a non-zero target taxon ID with --exclude-taxid." << endl;
            exit(1);
        }
        par.targetTaxId = externalTaxID;
    }
    
    TaxonomyWrapper *taxonomy = loadTaxonomy(dbDir, par.taxonomyPath);
    Reporter reporter(par, taxonomy);

    TaxID targetTaxID;
    if (externalTaxID == -1) {
        targetTaxID = -1;
    } else {
        targetTaxID = taxonomy->getInternalTaxID(externalTaxID);
        if (targetTaxID == -1) {
            cerr << "Error: Taxon ID " << externalTaxID << " not found in the taxonomy database." << endl;
            exit(1);
        }
    }

    vector<size_t> readIdxs;
    
    if (excludeMode) {
        cout << "Extracting reads not classified under taxon " << externalTaxID << " ... " << flush;
        reporter.getReadsNotClassifiedToClade(targetTaxID, classificationFileName, readIdxs);
    } else {
        cout << "Extracting reads classified under taxon " << externalTaxID << " ... " << flush;
        reporter.getReadsClassifiedToClade(targetTaxID, classificationFileName, readIdxs);
    }
    cout << "done." << endl;

    for(size_t i = 0; i < 10; i++) {
        if (i < readIdxs.size()) {
            cout << "  Read index " << readIdxs[i] << endl;
        }
    }

    string queryFileName = par.filenames[0];
    string outdirPath, baseName, extension;

    extractBaseNameAndExtension(queryFileName, outdirPath, baseName, extension);
    cout << "Base name       : " << baseName << endl;
    if (!par.outputDir.empty()) {
        outdirPath = par.outputDir + "/";
    } else {
        outdirPath = outdirPath + "/";
    }
    string outputTaxIdSuffix = excludeMode ? "not_" + to_string(externalTaxID) : to_string(externalTaxID);
    string outFileName = outdirPath + baseName + "_" + outputTaxIdSuffix;
    reporter.printSpecifiedReads(readIdxs, queryFileName, outFileName);
    cout << "Extracted file  : " << outFileName << endl;
    
    if (par.seqMode == 2) {
        cout << "Processing the second file ... " << endl;
        queryFileName = par.filenames[1];
        extractBaseNameAndExtension(queryFileName, outdirPath, baseName, extension);
        if (!par.outputDir.empty()) {
            outdirPath = par.outputDir + "/";
        } else {
            outdirPath = outdirPath + "/";
        }
        outFileName = outdirPath + baseName + "_" + outputTaxIdSuffix;
        reporter.printSpecifiedReads(readIdxs, queryFileName, outFileName);
        cout << "Extracted file 2: " << outFileName << endl;
    }

    delete taxonomy;
    return 0;
}
