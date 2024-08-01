#ifndef ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#define ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
#include "Command.h"

extern int build(int argc, const char **argv, const Command& command);
extern int updataDB(int argc, const char **argv, const Command& command);
extern int classify(int argc, const char **argv, const Command& command);
extern int filter(int argc, const char **argv, const Command& command);
extern int grade(int argc, const char **argv, const Command& command);
extern int gradeByCladeSize(int argc, const char **argv, const Command& command);
extern int seqHeader2TaxId(int argc, const char **argv, const Command& command);
extern int addToLibrary(int argc, const char **argv, const Command& command);
extern int applyThreshold(int argc, const char **argv, const Command& command);
extern int binning2report(int argc, const char **argv, const Command& command);
extern int filterByGenus(int argc, const char **argv, const Command& command);
extern int databaseReport(int argc, const char **argv, const Command& command);
extern int mapping2taxon(int argc, const char **argv, const Command& command);
extern int build_kaiju(int argc, const char **argv, const Command& command);
extern int expand_diffidx(int argc, const char **argv, const Command& command);
extern int expand_diffidx_func(int argc, const char **argv, const Command& command);
extern int ncbi2gtdb(int argc, const char **argv, const Command& command);
extern int build_func(int argc, const char **argv, const Command& command);

#endif //ADCLASSIFIER2_LOCALCOMMANDDECLARATIONS_H
