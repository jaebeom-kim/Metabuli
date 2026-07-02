#ifndef METABULI_CANDIDATE_DB_WRITER_H
#define METABULI_CANDIDATE_DB_WRITER_H

#include "DBWriter.h"
#include "common.h"

#include <cstdint>
#include <memory>
#include <string>

class CandidateDBWriter {
public:
    CandidateDBWriter(
        const std::string &dataFileName,
        unsigned int threads,
        size_t mode = 0);

    CandidateDBWriter(
        const std::string &dataFileName,
        const std::string &indexFileName,
        unsigned int threads,
        size_t mode = 0);

    ~CandidateDBWriter();

    void open(size_t bufferSize = SIZE_MAX);
    void close(bool merge = true, bool needsSort = true);

    void writeQuery(
        uint32_t queryId,
        const Query &query,
        unsigned int threadIdx = 0);

    bool isClosed() const;

    const std::string &getDataFileName() const {
        return dataFileName;
    }

    const std::string &getIndexFileName() const {
        return indexFileName;
    }

private:
    std::string dataFileName;
    std::string indexFileName;
    unsigned int threads;
    size_t mode;
    std::unique_ptr<DBWriter> writer;

    static std::string defaultIndexFileName(const std::string &dataFileName);
    static void serializeQuery(const Query &query, std::string &buffer);
};

#endif
