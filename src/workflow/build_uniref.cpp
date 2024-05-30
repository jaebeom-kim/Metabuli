#include "ProteinDbIndexer.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"

void setDefaults_build_uniref(LocalParameters & par){
    par.reducedAA = 0;
    par.spaceMask = "11111111";
    par.taxonomyPath = "" ;
    par.splitNum = 4096;
    par.maskProb = 0.9;
    par.maskMode = 1;
    par.bufferSize = 1'000'000'000;
    par.accessionLevel = 0;
    // Get current date
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = to_string(1900 + ltm->tm_year) + "-" + to_string(1 + ltm->tm_mon) + "-" + to_string(ltm->tm_mday);
    
    // Get random alphanumeric string fore dbName from current time
    srand(time(NULL));
    string randStr = to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int build_uniref(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build_uniref(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
  
    // If dbDirectory does not exist, create it
    if (!FileUtil::directoryExists(par.filenames[0].c_str())) {
        FileUtil::makeDir(par.filenames[0].c_str());
    }

    // Create index
    par.proteinDB = par.filenames[1];
    ProteinDbIndexer idxCre(par);
    idxCre.index();

    return 0;
}
