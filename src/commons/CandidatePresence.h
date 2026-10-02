#ifndef METABULI_CANDIDATE_PRESENCE_H
#define METABULI_CANDIDATE_PRESENCE_H

#include "common.h"
#include "CandidateDBReader.h"
#include "LocalParameters.h"

#include <cstdint>
#include <string>
#include <unordered_map>
#include <unordered_set>

// Load speciesId -> genome size from <dbDir>/species2genomeSize.tsv
// (col0 = speciesId, col2 = genome size), matching Classifier::parseSp2GenomeSize.
// Returns an empty map (and warns) if the file is absent, which disables the
// coverage-based gate.
std::unordered_map<TaxID, uint64_t> loadGenomeSizes(const std::string &dbDir);

// Method-0 presence detection: iterate top-hit average score + read count +
// unique-top ratio (and, when genome sizes are available, genome-coverage
// evenness) to a fixpoint over a species-candidate DB, returning the set of
// species judged truly present (the survivors). Shared by filter-candidates and
// classify-candidates (--community-refine) so the "is this species really here?"
// decision is computed in exactly one place.
std::unordered_set<TaxID> filterByScoreAndCoverage(
    CandidateDBReader &reader,
    size_t entryCount,
    const std::unordered_map<TaxID, uint64_t> &sp2genomeSize,
    const LocalParameters &par);

#endif // METABULI_CANDIDATE_PRESENCE_H
