#ifndef ADCLASSIFIER2_LOCALPARAMETERS_H
#define ADCLASSIFIER2_LOCALPARAMETERS_H

#include "Parameters.h"

const int CITATION_SPACEPHARER = CITATION_END;

class LocalParameters : public Parameters {
public:
    static const int DBTYPE_METABULI = 100;

    // static void initInstance() {
    //     new LocalParameters;
    // }

    LocalParameters();
    static int defaultRamUsage();

    // True when the query is provided as two separate paired-end files.
    // Interleaved input is still paired processing but comes from a single
    // file, so it uses the single-file positional layout.
    bool pairedFileInput() const { return seqMode == 2 && interleaved == 0; }
    static LocalParameters& getLocalInstance() {
        if (instance == NULL) {
            initParameterSingleton();
        }
        return static_cast<LocalParameters&>(LocalParameters::getInstance());
    }

    std::vector<MMseqsParameter*> classify;
    std::vector<MMseqsParameter*> classifyCandidates;
    std::vector<MMseqsParameter*> createCandidates;
    std::vector<MMseqsParameter*> filterCandidates;
    std::vector<MMseqsParameter*> viewCandidates;
    std::vector<MMseqsParameter*> groupGeneration;
    std::vector<MMseqsParameter*> extract;
    std::vector<MMseqsParameter*> filter;
    std::vector<MMseqsParameter*> exclusiontest_hiv;
    std::vector<MMseqsParameter*> seqHeader2TaxId;
    std::vector<MMseqsParameter*> grade;
    std::vector<MMseqsParameter*> addToLibrary;
    std::vector<MMseqsParameter*> build;
    std::vector<MMseqsParameter*> updateDB;
    std::vector<MMseqsParameter*> applyThreshold;
    std::vector<MMseqsParameter*> binning2report;
    std::vector<MMseqsParameter*> filterByGenus;
    std::vector<MMseqsParameter*> databaseReport;
    std::vector<MMseqsParameter*> mapping2taxon;
    std::vector<MMseqsParameter*> printInfo;
    std::vector<MMseqsParameter*> query2reference;
    std::vector<MMseqsParameter*> expand_diffidx;
    std::vector<MMseqsParameter*> taxdump;
    std::vector<MMseqsParameter*> accession2taxid;
    std::vector<MMseqsParameter*> editNames;
    std::vector<MMseqsParameter*> createnewtaxalist;
    std::vector<MMseqsParameter*> classifiedRefiner;
    std::vector<MMseqsParameter*> validateDatabase;
    std::vector<MMseqsParameter*> makeBenchmarkSet;
    std::vector<MMseqsParameter*> buildUnirefDb;
    std::vector<MMseqsParameter*> buildUnirefTree;
    std::vector<MMseqsParameter*> assignUniref;
    std::vector<MMseqsParameter*> createCommonKmerList;
    std::vector<MMseqsParameter*> mergeAssemblyFiles;
    std::vector<MMseqsParameter*> createTaxDb;
    std::vector<MMseqsParameter*> refineReport;

    // UniRef
    PARAMETER(UNIREF_NUMBERS)
    std::string unirefNumbers;

    // Superkingdom taxonomy id
    PARAMETER(VIRUS_TAX_ID)
    PARAMETER(BACTERIA_TAX_ID)
    PARAMETER(ARCHAEA_TAX_ID)

    // DB and classify
    PARAMETER(SKIP_REDUNDANCY)
    PARAMETER(VALIDATE_DB)
    PARAMETER(SYNCMER)
    PARAMETER(SMER_LEN)
    PARAMETER(KMER_FORMAT)
    PARAMETER(UNIREF_XML)
    PARAMETER(PARAM_CUSTOM_METAMER)
    PARAMETER(SPACE_MASK)

    // Classify
    PARAMETER(SEQ_MODE)
    PARAMETER(INTERLEAVED)
    PARAMETER(PRECISION_MODE)
    PARAMETER(MIN_SCORE)
    PARAMETER(HAMMING_MARGIN)
    PARAMETER(MIN_SP_SCORE)
    PARAMETER(TINFO_PATH)
    PARAMETER(RAM_USAGE)
    PARAMETER(PRINT_LOG)
    PARAMETER(MIN_AA_MATCH)
    PARAMETER(MIN_AA_MATCH_EUK)
    PARAMETER(MATCH_PER_KMER)
    PARAMETER(MIN_SS_MATCH)
    PARAMETER(TIE_RATIO)
    PARAMETER(PRINT_LINEAGE)
    PARAMETER(MAX_SHIFT)
    PARAMETER(GAP_PENALTY)
    PARAMETER(EM)
    PARAMETER(NEIGHBOR_KMERS)
    PARAMETER(PMD_KMER)
    PARAMETER(DISABLE_TRIMMING)
    PARAMETER(SCORE_MODE)
    PARAMETER(DB_TOTAL_LENGTH)
    PARAMETER(MAX_E_VALUE)
    PARAMETER(USE_ALL_MATCHES)
    PARAMETER(MAX_HDIST)
    PARAMETER(TIE_BRAKER)
    PARAMETER(TOP_SPECIES)
    PARAMETER(MAPPING_OUTPUT)
    PARAMETER(UNCLASSIFIED)
    PARAMETER(QUERY_FILE)

    // classify || refine-report
    PARAMETER(MIN_AVG_SCORE)
    // filter-candidates: coverage-based species filtering
    PARAMETER(MIN_ADJ_EVENNESS)
    PARAMETER(COV_USE_ALL_HITS)
    PARAMETER(MIN_COUNT)
    // filter-candidates: filter method selection + best-evidence/uniqueness filter
    PARAMETER(FILTER_METHOD)
    PARAMETER(MIN_STRONG_SCORE)
    PARAMETER(MIN_STRONG_READS)
    PARAMETER(MIN_UNIQUE_READS)
    PARAMETER(MIN_CLADE_COUNT)
    PARAMETER(MIN_CLADE_PROPORTION)
    PARAMETER(PRINT_FILTERED_RESULTS)

    // extract
    PARAMETER(TARGET_TAX_ID)
    PARAMETER(EXTRACT_MODE)
    PARAMETER(PARAM_OUTDIR)

    // Group generation
    PARAMETER(MIN_EDGE_WEIGHT)
    PARAMETER(MIN_VOTE_SCORE)
    PARAMETER(SCORE_COL)
    PARAMETER(WEIGHT_MODE)
    int weightMode;


    // DB build parameters
    PARAMETER(LIBRARY_PATH)
    PARAMETER(TAXONOMY_PATH)
    PARAMETER(IS_ASSEMBLY)
    PARAMETER(SPLIT_NUM)
    PARAMETER(BUFFER_SIZE)
    PARAMETER(ACCESSION_LEVEL)
    PARAMETER(DB_NAME)
    PARAMETER(DB_DATE)
    PARAMETER(CDS_INFO)
    PARAMETER(MAKE_LIBRARY)
    PARAMETER(GTDB)
    PARAMETER(VALIDATE_INPUT)
    PARAMETER(READING_FRAME)
    PARAMETER(STORE_KMER_POS)
    PARAMETER(REP_GENOME_LIST)
    PARAMETER(NO_MASK_TAXA)

    // DB updated parameters
    PARAMETER(NEW_TAXA)

    //  parameters
    PARAMETER(TEST_RANK)
    PARAMETER(TEST_TYPE)
    PARAMETER(READID_COL)
    PARAMETER(TAXID_COL)

    PARAMETER(PRINT_COLUMNS)
    PARAMETER(CLADE_RANK)
    PARAMETER(SKIP_SECONDARY)
    PARAMETER(TOP_HIT_ONLY)

    // Filter
    PARAMETER(PRINT_MODE)
    PARAMETER(CONTAM_LIST)

    // printInfo
    PARAMETER(INFO_BEGIN)
    PARAMETER(INFO_END)

    // expand_diffidx
    PARAMETER(KMER_BEGIN)
    PARAMETER(KMER_END)

    // classifiedRefiner
    PARAMETER(REMOVE_UNCLASSIFIED)
    PARAMETER(EXCLUDE_TAXID)
    PARAMETER(SELECT_TAXID)
    PARAMETER(SELECT_COLUMNS)
    PARAMETER(REPORT)
    PARAMETER(RANK)
    PARAMETER(HIGHER_RANK_FILE)

    // benchmark set
    PARAMETER(RANDOM_SEED)
    PARAMETER(ASSACC2TAXID)

    // Superkingdom taxonomy id
    int virusTaxId;
    int bacteriaTaxId;
    int archaeaTaxId;

    // DB and classify
    int skipRedundancy;
    int validateDb;
    int syncmer;
    int smerLen;
    int kmerFormat;
    std::string unirefXml;
    std::string customMetamer;

    // Classify
    int seqMode;
    int interleaved = 0;   // 1 = query file holds interleaved paired-end reads (implies paired processing)
    int precisionMode;
    float minScore;
    std::string spaceMask;
    uint8_t hammingMargin;
    float minSpScore;
    int ramUsage;
    int printLog;
    int matchPerKmer;
    int minSSMatch;
    float tieRatio;
    float thresholdK;
    float minVoteScr;
    int minEdgeWeight;
    int neighborKmers;
    int printLineage;
    int maxShift;
    int gapPenalty;
    bool em;
    int pmdKmer;
    int disableTrimming;
    int scoreMode;
    size_t dbTotalLength;
    double maxEValue;
    int useAllMatches;
    int maxHdist;
    int minAaMatch;
    int minAaMatchEuk;
    std::string priorityTaxa;
    int topSpecies;
    std::string mappingOutput;
    int unclassified = 0;      // --unclassified: also write out the reads left unclassified
    std::string queryFile;     // --query-file: original FASTA/Q read file(s) for --unclassified with classify-candidates
    bool candidateOnly; // create-candidates: write species-candidate DB and skip classification
    float minAvgScore;
    float minAdjEvenness;  // filter-candidates: remove species with adjustedEvenness below this
    int covUseAllHits;     // filter-candidates: 1 = aggregate coverage over all candidate hits, 0 = top hit per read only
    int minCount;          // filter-candidates (method 0): min top-hit reads to keep a species (0 = disabled)
    int filterMethod;      // filter-candidates: 0 = score+coverage, 1 = best-evidence+uniqueness
    float minStrongScore;  // filter-candidates (method 1): per-read idScore counted as strong evidence
    int minStrongReads;    // filter-candidates (method 1): min strong reads to keep a species
    int minUniqueReads;    // filter-candidates (method 1): min unique-top reads to keep a species
    int minCladeCount;
    float minCladeProportion;
    std::string outFilteredResults;

    // Extract
    int targetTaxId;
    int extractMode;
    std::string outputDir;

    // Database creation
    std::string tinfoPath;
    std::string libraryPath;
    std::string taxonomyPath;
    std::string dbName;
    std::string dbDate;
    int splitNum;
    size_t bufferSize;
    int accessionLevel;
    std::string cdsInfo;
    int makeLibrary;
    std::string assAcc2taxid;
    int gtdb;
    int validateInput;
    int readingFrame;
    int storeKmerPos;
    std::string repGenomeList;
    std::string noMaskTaxa;

    // DB updated parameters
    std::string newTaxa;

    // Test parameters
    std::string testRank;
    std::string testType;
    std::string printColumns;
    int readIdCol;
    int taxidCol;
    int scoreCol;
    int topHitOnly;  // evaluate-candidates: evaluate only the top-scoring candidate per read
    std::string cladeRank;
    int skipSecondary;

    // Add to library
    bool assembly;

    // Filter
    int printMode;
    std::string contamList;

    // printInfo
    size_t infoBegin;
    size_t infoEnd;
    size_t kmerBegin;
    size_t kmerEnd;

    bool removeUnclassified;
    std::string excludeTaxid;
    std::string selectTaxid;
    std::string selectColumns;
    bool report;
    std::string rank;
    int higherRankFile;
   
    // benchmark set
    int randomSeed;
    std::string assacc2taxid;

    void printParameters(const std::string &module, int argc,
                         const char* pargv[],
                         const std::vector<MMseqsParameter*> &par);
    
    void parseParameters(int argc, const char *pargv[], const Command &command, bool printPar, int parseFlags,
                        int outputFlags);

private:

    LocalParameters(LocalParameters const&);
    ~LocalParameters() {};
    void operator=(LocalParameters const&);
};

#endif //ADCLASSIFIER2_LOCALPARAMETERS_H
