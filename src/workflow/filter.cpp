#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "QueryFilter.h"

void setFilterDefaults(LocalParameters & par){
    par.seqMode = 2;
    par.reducedAA = 0;
    par.minScore = 0.7;
    par.minCoverage = 0;
    par.minSpScore = 0;
    par.spaceMask = "11111111";
    par.hammingMargin = 0;
    par.verbosity = 3;
    par.ramUsage = 128;
    par.minCoveredPos = 4;
    par.printLog = 0;
    par.maxGap = 0;
    par.taxonomyPath = "DBDIR/taxonomy/" ;
    par.minConsCnt = 4;
    par.minConsCntEuk = 9;
    par.eukaryotaTaxId = 2759;
    par.maskMode = 0;
    par.maskProb = 0.9;
    par.matchPerKmer = 4;
}

int filter(int argc, const char **argv, const Command& command)
{
    LocalParameters & par = LocalParameters::getLocalInstance();
    setFilterDefaults(par);
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    if (par.seqMode == 2) {
        if (!FileUtil::directoryExists(par.filenames[3].c_str())) {
            FileUtil::makeDir(par.filenames[3].c_str());
        }
    } else {
        if (!FileUtil::directoryExists(par.filenames[2].c_str())) {
            FileUtil::makeDir(par.filenames[2].c_str());
        }
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    cout << "Number of threads: " << par.threads << endl;

    QueryFilter * queryFilter = new QueryFilter(par);

    queryFilter->startClassify(par);

    delete classifier;
    return 0;
}