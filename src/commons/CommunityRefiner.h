#ifndef METABULI_COMMUNITY_REFINER_H
#define METABULI_COMMUNITY_REFINER_H

#include "common.h"
#include "CandidateDBReader.h"
#include "TaxonomyWrapper.h"
#include "LocalParameters.h"

#include <cstdint>
#include <memory>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Weighted, undirected graph over DB species. An edge (i,j) with weight w means w
// reads tied species i and j together (within the shared tie margin). The graph is
// built only from reads whose true species is likely NOT in the DB (top score
// below a ceiling), so its dense clusters are homology footprints of unknown
// organisms rather than present DB siblings.
struct SpeciesGraph {
    std::unordered_set<TaxID> nodes;
    std::unordered_map<TaxID, std::unordered_map<TaxID, double>> adj; // symmetric
    void addEdge(TaxID a, TaxID b, double w);
};

// Build the species tie graph from a species-candidate DB.
//   scoreCeiling  : skip a read whose top idScore >= this (its true species is
//                   likely in the DB, so it should not shape unknown footprints)
//   tieRatio/minScore : define each read's tie set via the shared tieMarginFactor
//   maxTieBreadth : skip reads tying more than this many species (conserved-region
//                   noise); 0 = no cap
SpeciesGraph buildSpeciesTieGraph(
    CandidateDBReader &reader,
    size_t entryCount,
    float scoreCeiling,
    float tieRatio,
    float minScore,
    size_t maxTieBreadth);

// ---- Community detection (pluggable via --community-method) ----
class CommunityDetector {
public:
    virtual ~CommunityDetector() = default;
    // Partition graph nodes into communities (each a list of >= 2 species).
    // Species in no surviving edge are not returned.
    virtual std::vector<std::vector<TaxID>> detect(const SpeciesGraph &graph) const = 0;
};

// method 0 (v1): keep edges with weight >= minEdgeWeight, return the connected
// components of size >= 2. Future methods (1, 2, ...) can add weighted
// community detection if the simple grouping over-merges through hub species.
std::unique_ptr<CommunityDetector> makeCommunityDetector(int method, double minEdgeWeight);

// ---- Community labeling (pluggable via --community-label) ----
struct CommunityMember {
    TaxID speciesId;
    double weight; // read support for this member (used by weighted labelers)
};

class CommunityLabeler {
public:
    virtual ~CommunityLabeler() = default;
    // Representative (anchor) taxID for a community.
    virtual TaxID label(const std::vector<CommunityMember> &members,
                        TaxonomyWrapper &taxonomy) const = 0;
};

// label 0 (v1): plain LCA of all members (weights ignored).
// Future label 1: weighted-majority-covering LCA (NcbiTaxonomy::weightedMajorityLCA),
// which uses the per-member weights to drop low-support outliers for a lower label.
std::unique_ptr<CommunityLabeler> makeCommunityLabeler(int labelMethod);

// ---- Orchestration (classify-candidates --community-refine) ----
// Detect present species (shared method-0 presence gates), build the species tie
// graph, detect communities, and rewrite each read's classification in place:
//   - tie set has >=1 present species  -> classify among the present ones (LCA if tied)
//   - tie set all removed, top in a community -> community label (recovers reads that
//        would otherwise be a false species call / discarded)
//   - tie set all removed, >=2 species, no community -> LCA of the tie set
//   - uniquely mapped to one removed species, no community -> unclassified
// queryList must carry speciesCandidates (as chooseBestTaxonFromCandidates leaves them).
// Reads a candidate DB (the original create-candidates DB) for presence + the graph.
void refineClassificationsWithCommunities(
    std::vector<Query> &queryList,
    CandidateDBReader &reader,
    size_t entryCount,
    const std::string &dbDir,
    TaxonomyWrapper &taxonomy,
    const LocalParameters &par);

#endif // METABULI_COMMUNITY_REFINER_H
