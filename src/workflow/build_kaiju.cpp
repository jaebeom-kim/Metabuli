#include "KaijuWrapper.h"
#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "LocalUtil.h"

void setDefaults_build_kaiju(LocalParameters & par){
    par.reducedAA = 0;
    par.spaceMask = "11111111";
    par.splitNum = 4096;
    par.maskProb = 0.9;
    par.maskMode = 1;
    par.bufferSize = 1'000'000'000;
    par.accessionLevel = 0;
    // Get current date
    time_t now = time(0);
    tm *ltm = localtime(&now);
    par.dbDate = std::to_string(1900 + ltm->tm_year) + "-" + std::to_string(1 + ltm->tm_mon) + "-" + std::to_string(ltm->tm_mday);
    
    // Get random alphanumeric string fore dbName from current time
    srand(time(NULL));
    std::string randStr = std::to_string(rand());
    par.dbName = randStr.substr(0, 32);
}

int build_kaiju(int argc, const char **argv, const Command &command){
    // Load parameters
    LocalParameters &par = LocalParameters::getLocalInstance();
    setDefaults_build_kaiju(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);
  
    // If dbDirectory does not exist, create it
    if (!FileUtil::directoryExists(par.filenames[0].c_str())) {
        FileUtil::makeDir(par.filenames[0].c_str());
    }

    // Create index
    std::string outdir = par.filenames[0];
    std::string proteinFastaName = par.filenames[1];

    std::vector<std::string> unirefIds;
    std::vector<TaxID> unirefTaxIds;
    KaijuWrapper kaijuWrapper(par, unirefIds, unirefTaxIds);
    kaijuWrapper.buildKaijuIndex(proteinFastaName, outdir);

    // kaijuWrapper.makeBWT(proteinFastaName, outdir + "/kaiju_mtbl");

    // // print unirefIds
    // for (auto &id : unirefIds) {
    //     std::cout << id << std::endl;
    // }

    // LocalUtil::exportVector(unirefIds, outdir + "/kaiju_mtbl.UniRefIds");



    return 0;
}
