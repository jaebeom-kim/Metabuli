#include "CommunityRefiner.h"
#include "CandidatePresence.h"

#include <algorithm>
#include <iostream>
#include <queue>

// ---------------------------------------------------------------------------
// SpeciesGraph
// ---------------------------------------------------------------------------
void SpeciesGraph::addEdge(TaxID a, TaxID b, double w) {
    if (a == b) {
        return;
    }
    nodes.insert(a);
    nodes.insert(b);
    adj[a][b] += w;
    adj[b][a] += w;
}

SpeciesGraph buildSpeciesTieGraph(
    CandidateDBReader &reader,
    size_t entryCount,
    float scoreCeiling,
    float tieRatio,
    float minScore,
    size_t maxTieBreadth)
{
    SpeciesGraph graph;
    std::vector<TaxID> tieSet;
    CandidateDBEntry entry;
    for (size_t i = 0; i < entryCount; ++i) {
        if (!reader.getByIndex(i, entry, 0)) {
            continue;
        }

        // Top hit = first candidate at/above minScore (candidates are best-first).
        float bestScore = -1.0f;
        for (const SpeciesCandidate &c : entry.candidates) {
            if (c.idScore >= minScore) { bestScore = c.idScore; break; }
        }
        if (bestScore < 0.0f) {
            continue; // no candidate clears minScore
        }
        // Skip reads whose true species is likely in the DB: only moderate-score
        // reads shape unknown-organism footprints.
        if (bestScore >= scoreCeiling) {
            continue;
        }

        // Collect the tie set with the same scaled margin classify uses.
        const float threshold = bestScore * tieMarginFactor(bestScore, tieRatio);
        tieSet.clear();
        for (const SpeciesCandidate &c : entry.candidates) {
            if (c.idScore < minScore) {
                continue;
            }
            if (c.idScore >= threshold) {
                tieSet.push_back(c.speciesId);
            } else {
                break; // best-first: the rest are below the margin
            }
        }

        if (tieSet.size() < 2) {
            continue; // unique mapping -> no tie edge
        }
        if (maxTieBreadth > 0 && tieSet.size() > maxTieBreadth) {
            continue; // conserved-region read tying many species -> noise
        }

        // One read co-ties every pair in its tie set (+1 per pair).
        for (size_t a = 0; a < tieSet.size(); ++a) {
            for (size_t b = a + 1; b < tieSet.size(); ++b) {
                graph.addEdge(tieSet[a], tieSet[b], 1.0);
            }
        }
    }
    return graph;
}

// ---------------------------------------------------------------------------
// Community detection
// ---------------------------------------------------------------------------
namespace {

// method 0: edges with weight >= minEdgeWeight define connected components;
// components of size >= 2 are returned as communities.
class ThresholdConnectedDetector : public CommunityDetector {
public:
    explicit ThresholdConnectedDetector(double minEdgeWeight)
        : minEdgeWeight(minEdgeWeight) {}

    std::vector<std::vector<TaxID>> detect(const SpeciesGraph &graph) const override {
        std::vector<std::vector<TaxID>> communities;
        std::unordered_set<TaxID> visited;
        for (const TaxID start : graph.nodes) {
            if (visited.count(start) != 0) {
                continue;
            }
            // BFS over edges that pass the weight threshold.
            std::vector<TaxID> component;
            std::queue<TaxID> frontier;
            frontier.push(start);
            visited.insert(start);
            while (!frontier.empty()) {
                const TaxID cur = frontier.front();
                frontier.pop();
                component.push_back(cur);
                const auto adjIt = graph.adj.find(cur);
                if (adjIt == graph.adj.end()) {
                    continue;
                }
                for (const auto &nb : adjIt->second) {
                    if (nb.second < minEdgeWeight) {
                        continue; // edge too weak
                    }
                    if (visited.insert(nb.first).second) {
                        frontier.push(nb.first);
                    }
                }
            }
            if (component.size() >= 2) {
                communities.push_back(std::move(component));
            }
        }
        return communities;
    }

private:
    double minEdgeWeight;
};

// label 0: plain LCA of all members.
class LcaLabeler : public CommunityLabeler {
public:
    TaxID label(const std::vector<CommunityMember> &members,
                TaxonomyWrapper &taxonomy) const override {
        if (members.empty()) {
            return 0;
        }
        std::vector<TaxID> ids;
        ids.reserve(members.size());
        for (const CommunityMember &m : members) {
            ids.push_back(m.speciesId);
        }
        const TaxonNode *node = taxonomy.LCA(ids);
        return node == nullptr ? 0 : node->taxId;
    }
};

} // namespace

std::unique_ptr<CommunityDetector> makeCommunityDetector(int method, double minEdgeWeight) {
    switch (method) {
        case 0:
        default:
            return std::unique_ptr<CommunityDetector>(new ThresholdConnectedDetector(minEdgeWeight));
    }
}

std::unique_ptr<CommunityLabeler> makeCommunityLabeler(int labelMethod) {
    switch (labelMethod) {
        case 0:
        default:
            return std::unique_ptr<CommunityLabeler>(new LcaLabeler());
    }
}

// ---------------------------------------------------------------------------
// Orchestration
// ---------------------------------------------------------------------------
void refineClassificationsWithCommunities(
    std::vector<Query> &queryList,
    CandidateDBReader &reader,
    size_t entryCount,
    const std::string &dbDir,
    TaxonomyWrapper &taxonomy,
    const LocalParameters &par)
{
    std::cout << "Community refinement (--community-refine)" << std::endl;

    // 1) Present species (shared method-0 gates).
    const std::unordered_map<TaxID, uint64_t> sp2genomeSize = loadGenomeSizes(dbDir);
    const std::unordered_set<TaxID> survivors =
        filterByScoreAndCoverage(reader, entryCount, sp2genomeSize, par);

    // 2) Species tie graph -> communities.
    const size_t maxBreadth = par.communityMaxBreadth < 0 ? 0 : static_cast<size_t>(par.communityMaxBreadth);
    SpeciesGraph graph = buildSpeciesTieGraph(
        reader, entryCount, par.novelScoreCeiling, par.tieRatio, par.minScore, maxBreadth);
    std::unique_ptr<CommunityDetector> detector =
        makeCommunityDetector(par.communityMethod, static_cast<double>(par.communityMinEdge));
    const std::vector<std::vector<TaxID>> communities = detector->detect(graph);
    std::unique_ptr<CommunityLabeler> labeler = makeCommunityLabeler(par.communityLabel);

    // 3) speciesId -> community label (taxID). Weights are read support per member.
    std::unordered_map<TaxID, double> memberWeight; // # reads where a species is top/tied in the graph
    std::unordered_map<TaxID, TaxID> sp2communityLabel;
    size_t labelledCommunities = 0;
    for (const std::vector<TaxID> &members : communities) {
        std::vector<CommunityMember> cm;
        cm.reserve(members.size());
        for (const TaxID m : members) {
            cm.push_back({m, memberWeight.count(m) ? memberWeight[m] : 1.0});
        }
        const TaxID lbl = labeler->label(cm, taxonomy);
        if (lbl == 0) {
            continue;
        }
        for (const TaxID m : members) {
            sp2communityLabel[m] = lbl;
        }
        ++labelledCommunities;
    }

    // 4) Per-read refined assignment (queryList carries speciesCandidates).
    size_t recovered = 0, rescuedLca = 0, droppedFalse = 0;
    for (Query &q : queryList) {
        const std::vector<SpeciesCandidate> &cands = q.speciesCandidates;

        float best = -1.0f;
        for (const SpeciesCandidate &c : cands) {
            if (c.idScore >= par.minScore) { best = c.idScore; break; }
        }
        if (best < 0.0f) {
            continue; // nothing passes minScore; leave the baseline call
        }
        const float threshold = best * tieMarginFactor(best, par.tieRatio);

        TaxID topSp = 0;
        std::vector<TaxID> tieAll, tieSurvivors;
        for (const SpeciesCandidate &c : cands) {
            if (c.idScore < par.minScore) {
                continue;
            }
            if (c.idScore >= threshold) {
                if (topSp == 0) { topSp = c.speciesId; }
                tieAll.push_back(c.speciesId);
                if (survivors.count(c.speciesId) != 0) {
                    tieSurvivors.push_back(c.speciesId);
                }
            } else {
                break; // best-first: the rest are below the margin
            }
        }
        if (tieAll.empty()) {
            continue;
        }

        // Only intervene when the read has NO present species in its tie set.
        // Reads with present support keep their baseline call, which preserves
        // any finer (e.g. strain-level) resolution the baseline achieved.
        if (!tieSurvivors.empty()) {
            continue;
        }

        // Tie set is entirely removed (artifact) species -> recover or drop.
        const auto it = sp2communityLabel.find(topSp);
        if (it != sp2communityLabel.end()) {
            q.classification = it->second;              // unknown-organism community label
            ++recovered;
        } else if (tieAll.size() >= 2) {
            q.classification = taxonomy.LCA(tieAll)->taxId; // tied removed species -> their LCA
            ++rescuedLca;
        } else {
            q.classification = 0;                       // uniquely mapped to one removed species, no community
            ++droppedFalse;
        }
    }

    std::cout << "  Present species                                       : " << survivors.size() << std::endl;
    std::cout << "  Communities                                           : " << labelledCommunities << std::endl;
    std::cout << "  Reads recovered to a community                        : " << recovered << std::endl;
    std::cout << "  Reads rescued to a removed-tie LCA                    : " << rescuedLca << std::endl;
    std::cout << "  Reads dropped (uniquely mapped to a removed species)  : " << droppedFalse << std::endl;
}
