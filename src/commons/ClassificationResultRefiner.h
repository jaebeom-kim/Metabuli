#ifndef METABULI_CLASSIFICATION_RESULT_REFINER_H
#define METABULI_CLASSIFICATION_RESULT_REFINER_H

#include "LocalParameters.h"
#include "TaxonomyWrapper.h"

#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

class ClassificationResultRefiner {
public:
    explicit ClassificationResultRefiner(const LocalParameters &par);

    int runRefineReport(const std::string &resultFileName,
                        const std::string &dbDir,
                        const std::string &outFileName) const;

    int runRefineResult(const std::string &classifiedFile) const;

private:
    struct ParsedClassification {
        bool isClassified = false;
        std::string readId;
        TaxID taxonomyId = 0;
        int effectiveReadLength = 0;
        float dnaIdentityScore = 0.0f;
        double evalue = -1.0;
        std::string classificationRank;
        std::string taxIdKmerCounts;
        std::string fullLineage;
    };

    const LocalParameters &par;

    static ParsedClassification parseFields(const std::vector<std::string> &fields,
                                            bool hasEvalue,
                                            bool hasLineage);

    static bool isInAnyClade(TaxonomyWrapper *taxonomy,
                             const std::vector<TaxID> &cladeTaxIds,
                             TaxID taxId);

    static std::unordered_set<TaxID> getSelfAndChildren(
        const std::vector<TaxID> &taxIds,
        const std::unordered_map<TaxID, std::vector<TaxID>> &parentToChildren);
};

#endif
