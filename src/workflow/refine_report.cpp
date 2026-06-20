#include "ClassificationResultRefiner.h"
#include "LocalParameters.h"

void setRefineReportDefault(LocalParameters &par) {
    par.minAvgScore = 0;
    par.minCladeCount = 0;
    par.minCladeProportion = 0;
    par.outFilteredResults = "";
    par.maxEValue = 0;
    par.excludeTaxid = "";
    par.minScore = 0;
}

int refine_report(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    ClassificationResultRefiner refiner(par);
    return refiner.runRefineReport(par.filenames[0], par.filenames[1], par.filenames[2]);
}
