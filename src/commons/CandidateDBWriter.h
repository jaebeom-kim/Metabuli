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

    // seqMode is the --seq-mode the candidates were generated with (1/2/3);
    // pass 0 when unknown (e.g. rewriting records from a pre-v3 DB).
    void writeQuery(
        uint32_t queryId,
        const Query &query,
        unsigned int threadIdx,
        uint32_t seqMode);

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
    static void serializeQuery(const Query &query, uint32_t seqMode, std::string &buffer);
};

#endif
