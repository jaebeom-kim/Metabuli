#include "Classifier.h"
#include "Parameters.h"
#include "LocalParameters.h"
#include "FileUtil.h"
#include "common.h"
#include "LocalUtil.h"
#include "fastq.h"
#include "validateDatabase.h"

// Defined in fastq_info.cpp (compiled as part of classify.cpp). Declared here
// to avoid including the .cpp twice, which would produce duplicate symbols.
FASTQ_FILE* validate_single_fastq_file(const char *f);

void setClassifyDefaults(LocalParameters &par);

int createCandidates(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    setClassifyDefaults(par);
    par.candidateOnly = true;
    par.storeKmerPos = 0;
    par.em = false;
    par.topSpecies = 5; // default number of candidate species per read; overridable
    par.parseParameters(argc, argv, command, true, Parameters::PARSE_ALLOW_EMPTY, 0);

    // Interleaved paired-end reads: single file, paired processing.
    if (par.interleaved) {
        par.seqMode = 2;
    }

    if (par.topSpecies <= 0) {
        std::cout << "Warning: --top-species must be >= 1 for create-candidates. Using 5." << std::endl;
        par.topSpecies = 5;
    }

    // Positional layout:
    //   --seq-mode 2 (paired) : <read_1> <read_2> <DB dir> <candidate DB>
    //   --seq-mode 1/3        : <read>            <DB dir> <candidate DB>
    const int dbIdx = 1 + par.pairedFileInput();
    const int outIdx = 2 + par.pairedFileInput();
    const std::string dbDir = par.filenames[dbIdx];
    const std::string candidateDbOut = par.filenames[outIdx];
    par.mappingOutput = candidateDbOut;

    // Validate query file(s)
    if (!LocalUtil::isValidQueryFile(par.filenames[0])) {
        std::cout << "Error: " << par.filenames[0] << " is not a valid query file." << std::endl;
        std::cout << "       Allowed extensions are .fna, .fasta, .fa, .fq, .fastq, and their gzip versions (e.g., .fna.gz)" << std::endl;
        return 1;
    }
    if (par.pairedFileInput()) {
        if (FileUtil::directoryExists(par.filenames[1].c_str())) {
            std::cout << "Error: " << par.filenames[1] << " is a directory. Please specify a query file name." << std::endl;
            std::cout << "       For '--seq-mode 2', please provide two query files." << std::endl;
            return 1;
        }
        if (!LocalUtil::isValidQueryFile(par.filenames[1])) {
            std::cout << "Error: " << par.filenames[1] << " is not a valid query file." << std::endl;
            std::cout << "       Allowed extensions are .fna, .fasta, .fa, .fq, .fastq, and their gzip versions (e.g., .fna.gz)" << std::endl;
            return 1;
        }
    }

    if (par.validateInput) {
        const std::vector<std::string> queryFiles = par.pairedFileInput()
            ? std::vector<std::string>{par.filenames[0], par.filenames[1]}
            : std::vector<std::string>{par.filenames[0]};
        for (const std::string &queryFile : queryFiles) {
            if (LocalUtil::isFasta(queryFile)) {
                std::cout << "Validating FASTA file: " << queryFile << std::endl;
                if (validate_fasta_file(queryFile.c_str(), 1) != 0) {
                    std::cout << "Error: " << queryFile << " is not a valid FASTA file." << std::endl;
                    return 1;
                }
            } else if (LocalUtil::isFastq(queryFile)) {
                std::cout << "Validating FASTQ file: " << queryFile << std::endl;
                FASTQ_FILE *fd = validate_single_fastq_file(queryFile.c_str());
                fastq_destroy(fd);
            }
        }
    }

    // Validate database directory
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
        if (!FileUtil::fileExists(taxonomyDb.c_str())
            && par.taxonomyPath.empty()
            && !FileUtil::fileExists((dbDir + "/taxonomy/merged.dmp").c_str())) {
            std::cout << "Error: taxonomy files are not found." << std::endl;
            std::cout << "       One of the followings should be provided:" << std::endl;
            std::cout << "       1. File: " << taxonomyDb << std::endl;
            std::cout << "       2. Dir : " << dbDir + "/taxonomy" << std::endl;
            std::cout << "       3. Specify --taxonomy-path" << std::endl;
            return 1;
        }
    }

    // Make sure the output directory exists
    const std::string outParent = FileUtil::dirName(candidateDbOut);
    if (!outParent.empty() && !FileUtil::directoryExists(outParent.c_str())) {
        FileUtil::makeDir(outParent.c_str());
    }

#ifdef OPENMP
    omp_set_num_threads(par.threads);
#endif

    Classifier *classifier = new Classifier(par);
    std::cout << "Generating species-candidate DB ..." << std::endl;
    classifier->generateCandidates();
    delete classifier;
    return 0;
}
