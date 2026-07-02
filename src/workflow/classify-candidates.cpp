#include "Classifier.h"
#include "Command.h"
#include "FileUtil.h"
#include "LocalParameters.h"
#include "validateDatabase.h"

void setClassifyDefaults(LocalParameters &par);

int classifyCandidates(int argc, const char **argv, const Command& command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setClassifyDefaults(par);
    par.seqMode = 1;
    par.storeKmerPos = 0;
    par.topSpecies = 0;
    par.mappingOutput = "";
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string &candidateDb = par.filenames[0];
    const std::string &dbDir = par.filenames[1];
    const std::string &outDir = par.filenames[2];

    if (!FileUtil::fileExists(candidateDb.c_str())) {
        std::cout << "Error: " << candidateDb << " is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::fileExists((candidateDb + ".index").c_str())) {
        std::cout << "Error: " << candidateDb + ".index" << " is not found." << std::endl;
        return 1;
    }
    if (!FileUtil::directoryExists(outDir.c_str())) {
        FileUtil::makeDir(outDir.c_str());
    }

    if (par.validateDb) {
        if (validateDatabase(dbDir) != 0) {
            std::cout << "Error: Database validation failed." << std::endl;
            return 1;
        }
    } else {
        if (!FileUtil::directoryExists(dbDir.c_str())) {
            std::cout << "Error: " << dbDir << " is not found." << std::endl;
            return 1;
        }

        const std::string taxonomyDb = dbDir + "/taxonomyDB";
        const std::string taxonomyDir = dbDir + "/taxonomy";
        if (!FileUtil::fileExists(taxonomyDb.c_str())
            && par.taxonomyPath.empty()
            && !FileUtil::fileExists((taxonomyDir + "/merged.dmp").c_str())) {
            std::cout << "Error: taxonomy files are not found." << std::endl;
            std::cout << "       One of the followings should be provided:" << std::endl;
            std::cout << "       1. File: " << taxonomyDb << std::endl;
            std::cout << "       2. Dir : " << taxonomyDir << std::endl;
            std::cout << "       3. Specify --taxonomy-path" << std::endl;
            return 1;
        }
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Classifier classifier(par);
    return classifier.classifyCandidates(candidateDb) ? 0 : 1;
}
