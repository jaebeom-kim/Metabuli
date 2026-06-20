#include "ClassificationResultRefiner.h"
#include "Debug.h"
#include "FileUtil.h"
#include "LocalParameters.h"

void classifiedRefinerDefault(LocalParameters &par) {
    par.removeUnclassified = false;
    par.excludeTaxid = "";
    par.selectTaxid = "";
    par.selectColumns = "";
    par.report = false;
    par.rank = "";
    par.higherRankFile = 0;
    par.minScore = 0.0;
    par.maxEValue = 0;
}

int classifiedRefiner(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    classifiedRefinerDefault(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    const std::string &classifiedFile = par.filenames[0];
    const std::string &taxonomyDir = par.filenames[1];

    if (!FileUtil::fileExists(classifiedFile.c_str())) {
        Debug(Debug::INFO) << "Classified file" << classifiedFile << " is NOT exists.\n";
        return 0;
    }

    if (!FileUtil::directoryExists(taxonomyDir.c_str())) {
        Debug(Debug::INFO) << "Taxonomy dump" << taxonomyDir << " is NOT exists.\n";
        return 0;
    }

    ClassificationResultRefiner refiner(par);
    return refiner.runRefineResult(classifiedFile);
}
