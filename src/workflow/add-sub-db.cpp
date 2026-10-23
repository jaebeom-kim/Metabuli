#include "IndexCreator.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "accession2taxid.h"
#include "editNames.h"
#include "fasta_validate.h"

#include <cstdio>
#include <fstream>

// add-sub-db builds a k-mer list from new sequences and stores it as a
// subdirectory of an existing database instead of merging it into the
// database's main list (what updateDB does). The sub-DB is built against the
// database's taxonomyDB, so every list under the database directory shares one
// internal tax-ID numbering and classify/create-candidates can search any
// subset of them (--sub-dbs).

void setDefaults_addSubDb(LocalParameters & par){
    par.dbTotalLength = 0;
    par.noMaskTaxa = "";
    par.makeLibrary = 0;
    par.gtdb = 0;
    par.validateInput = 0;
    par.validateDb = 0;
    par.kmerFormat = 1;
    par.skipRedundancy = 0;
    par.ramUsage = LocalParameters::defaultRamUsage();
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    par.maskProb = 0.9;
    par.maskMode = 1;
    par.packInfo = 1;
    // Get current date
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);

    // Get random alphanumeric string fore dbName from current time
    srand(time(NULL));
    string randStr = to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int addSubDb(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_addSubDb(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string dbDir = par.filenames[0];
    const string subDbName = par.filenames[3];

    // The name becomes a --sub-dbs selector and a directory name.
    if (subDbName.empty() || subDbName == "." || subDbName == ".."
        || subDbName.find('/') != string::npos
        || subDbName.find(',') != string::npos) {
        cout << "Error: invalid sub-DB name \"" << subDbName << "\"." << endl;
        cout << "       It must not contain '/' or ','." << endl;
        return 1;
    }
    // A pre-kmer-lists database keeps its top-level list addressable as
    // "base", so that selector is taken.
    if (subDbName == SUBDB_BASE_NAME && fileExist(dbDir + "/diffIdx")) {
        cout << "Error: \"" << SUBDB_BASE_NAME << "\" is the selector of this database's top-level k-mer list." << endl;
        cout << "       Choose another name." << endl;
        return 1;
    }
    const string subDir = dbDir + "/" + SUBDB_DIR_NAME + "/" + subDbName;
    if (fileExist(subDir + "/diffIdx")) {
        cout << "Error: sub-DB \"" << subDbName << "\" already exists in " << dbDir << "." << endl;
        return 1;
    }
    if (!fileExist(dbDir + "/taxonomyDB")) {
        cout << "Error: " << dbDir << " has no taxonomyDB file." << endl;
        cout << "       Sub-DBs share the internal tax-ID numbering of a serialized taxonomyDB." << endl;
        return 1;
    }

    // One or more accession2taxid files, comma-separated; plain or gzipped.
    {
        bool missing = false;
        for (const string & f : Util::split(par.filenames[2], ",")) {
            if (!FileUtil::fileExists(f.c_str())) {
                cout << "Error: accession2taxid file " << f << " does not exist." << endl;
                missing = true;
            }
        }
        if (missing) {
            return 1;
        }
    }

    if (par.validateInput == 1) {
        const string & fnaListFileName = par.filenames[1];
        // Read the file line by line
        ifstream file(fnaListFileName);
        if (!file.is_open()) {
            cout << "Error: Unable to open file " << fnaListFileName << endl;
            return 1;
        }
        string singleFastaName;
        while (getline(file, singleFastaName)) {
            if (!LocalUtil::isFasta(singleFastaName)) {
                cout << "Error: " << singleFastaName << " has a unsupported extension." << endl;
                cout << "       Allowed extensions are .fna, .fasta, .fa, and their gzip versions (e.g., .fna.gz)" << endl;
                file.close();
                return 1;
            }
            if (!FileUtil::fileExists(singleFastaName.c_str())) {
                cout << "Error: " << singleFastaName << " does not exist." << endl;
                file.close();
                return 1;
            }
            int validateRes = validate_fasta_file(singleFastaName.c_str(), 1);
            if (validateRes != 0) {
                cout << "Error: " << singleFastaName << " is not a valid FASTA file." << endl;
                cout << "       Please check the entries in the file" << endl;
                file.close();
                return 1;
            }
        }
    }

    if (par.gtdb == 1 && par.filenames[2].find(',') != string::npos) {
        cout << "Error: --gtdb takes a single taxid map file, not a comma-separated list." << endl;
        return 1;
    }
    if (par.gtdb == 1) {
        accession2taxid(par.filenames[1], par.filenames[2]);
        par.filenames[2] = par.filenames[2].substr(0, par.filenames[2].find_last_of('.')) + ".accession2taxid";
    }

    // Load the database's taxonomy
    Debug(Debug::INFO) << "Loading taxonomy DB from " << dbDir << " ... ";
    TaxonomyWrapper * taxonomy = loadTaxonomy(dbDir);
    Debug(Debug::INFO) << "done.\n";

    // Make a new taxonomy DB if new taxa are added. Taxa are only appended;
    // existing internal IDs never change, so other sub-DBs stay valid.
    if (!par.newTaxa.empty()) {
        Debug(Debug::INFO) << "Adding new taxa to the taxonomy DB.\n";
        taxonomy->checkNewTaxa(par.newTaxa);
        std::vector<NewTaxon> newTaxaList;
        TaxonomyWrapper::getListOfTaxa(par.newTaxa, newTaxaList);
        TaxonomyWrapper * newTaxonomy = taxonomy->addNewTaxa(newTaxaList);
        delete taxonomy;
        taxonomy = newTaxonomy;
        Debug(Debug::INFO) << "New taxonomy generated.\n";
    }

    loadDbParameters(par, dbDir);
    // loadDbParameters imported the database's Total_seq_length, which
    // IndexCreator would use as the starting total. The sub-DB must record
    // only its own sequences; classify sums the lengths of the selected
    // sub-DBs itself.
    par.dbTotalLength = 0;
    // Databases from Metabuli versions with accession-level support carry one
    // synthetic taxon per accession; their numbering conventions predate
    // sub-DBs, so refuse to extend them.
    if (readDbParameter(dbDir, "Accession_level") == "1") {
        cout << "Error: " << dbDir << " is an accession-level database (a removed feature)." << endl;
        cout << "       add-sub-db does not support accession-level databases." << endl;
        return 1;
    }

    // Lineage ID: identifies the internal numbering all sub-DBs must share.
    // Minted once per database; verified against each sub-DB at load time.
    string lineage = readDbParameter(dbDir, SUBDB_KEY_LINEAGE);
    if (lineage.empty()) {
        lineage = randomHexString(16);
        appendDbParameter(dbDir, SUBDB_KEY_LINEAGE, lineage);
    }

    const string listsDir = dbDir + "/" + SUBDB_DIR_NAME;
    if (!FileUtil::directoryExists(listsDir.c_str())) {
        FileUtil::makeDir(listsDir.c_str());
    }
    if (!FileUtil::directoryExists(subDir.c_str())) {
        FileUtil::makeDir(subDir.c_str());
    }

    // IndexCreator derives every output path from filenames[0].
    par.filenames[0] = subDir;
    IndexCreator idxCre(par, taxonomy, par.kmerFormat);
    idxCre.setIsUpdating(true);
    // Unlike updateDB, no pre-existing k-mer lists are merged in, so the info
    // index only ever holds this run's writable tax IDs; let it pack to their
    // maximum instead of the taxonomy-wide maxTaxID.
    idxCre.setCapInfoIdsToObserved(true);
    idxCre.createIndex();

    // Finalize the sub-DB's k-mer list from its flushed chunks. Unlike
    // updateDB, the database's existing lists are NOT merged in.
    cout << "Finalize sub-DB files" << endl;
    idxCre.updateTaxId2SpeciesTaxId(subDir + "/taxID_list");
    idxCre.printFilesToMerge();
    idxCre.setMergedFileNames(subDir + "/diffIdx", subDir + "/info", subDir + "/split");
    idxCre.mergeTargetFiles<FilterMode::DB_CREATION>();

    // Persist the (possibly extended) taxonomy. Write-to-temp + rename so a
    // crash never leaves the database without a valid taxonomyDB, and the
    // mmapped old file is never truncated in place.
    if (!par.newTaxa.empty() && !taxonomy->IsExternalData()) {
        const string tmpTaxonomy = dbDir + "/taxonomyDB.tmp";
        taxonomy->writeTaxonomyDB(tmpTaxonomy);
        if (rename(tmpTaxonomy.c_str(), (dbDir + "/taxonomyDB").c_str()) != 0) {
            cout << "Error: could not replace " << dbDir << "/taxonomyDB" << endl;
            return 1;
        }
    }

    // Move the new accession mappings (external tax IDs) up to the
    // database-level acc2taxid.map; it is shared by all sub-DBs.
    {
        ifstream newMappings(subDir + "/acc2taxid.map");
        if (newMappings.is_open()) {
            ofstream dbMappings(dbDir + "/acc2taxid.map", std::ios::app);
            dbMappings << newMappings.rdbuf();
            dbMappings.close();
            newMappings.close();
            remove((subDir + "/acc2taxid.map").c_str());
        }
    }

    // Guards checked by classify/create-candidates before searching this
    // sub-DB together with others.
    appendDbParameter(subDir, SUBDB_KEY_LINEAGE, lineage);
    appendDbParameter(subDir, SUBDB_KEY_MAX_TAXID, to_string(taxonomy->getMaxTaxID()));

    delete taxonomy;
    cout << "Sub-DB \"" << subDbName << "\" created in " << dbDir << "." << endl;
    cout << "Search it with --sub-dbs " << subDbName << " (or together with others; default searches all)." << endl;

    if (par.validateDb) {
        cout << "Note: --validate-db is not supported for sub-DBs yet; skipping." << endl;
    }

    return 0;
}
