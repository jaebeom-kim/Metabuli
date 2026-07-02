#ifndef METABULI_CANDIDATE_DB_READER_H
#define METABULI_CANDIDATE_DB_READER_H

#include "DBReader.h"
#include "common.h"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

struct CandidateDBEntry {
    uint32_t queryId = 0;
    uint32_t queryLength = 0;
    std::string queryName;
    std::vector<SpeciesCandidate> candidates;
};

class CandidateDBReader {
public:
    CandidateDBReader(
        const std::string &dataFileName,
        int threads);

    CandidateDBReader(
        const std::string &dataFileName,
        const std::string &indexFileName,
        int threads);

    ~CandidateDBReader();

    bool open(int sortMode = DBReader<unsigned int>::SORT_BY_ID);
    void close();

    size_t size() const;
    bool getByIndex(size_t index, CandidateDBEntry &entry, int threadIdx = 0);
    bool getByQueryId(uint32_t queryId, CandidateDBEntry &entry, int threadIdx = 0);

    const std::string &getDataFileName() const {
        return dataFileName;
    }

    const std::string &getIndexFileName() const {
        return indexFileName;
    }

private:
    std::string dataFileName;
    std::string indexFileName;
    int threads;
    std::unique_ptr<DBReader<unsigned int>> reader;

    static std::string defaultIndexFileName(const std::string &dataFileName);
    static bool deserializeEntry(
        uint32_t queryId,
        const char *data,
        size_t dataSize,
        CandidateDBEntry &entry);
};

#endif
