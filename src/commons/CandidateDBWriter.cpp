#include "CandidateDBWriter.h"

#include "Parameters.h"

namespace {
constexpr uint32_t CANDIDATE_RECORD_MAGIC = 0x444e4143; // CAND
constexpr uint32_t CANDIDATE_RECORD_VERSION = 1;

template <typename T>
void appendPod(std::string &buffer, const T &value) {
    const char *data = reinterpret_cast<const char *>(&value);
    buffer.append(data, sizeof(T));
}
}

CandidateDBWriter::CandidateDBWriter(
    const std::string &dataFileName,
    unsigned int threads,
    size_t mode)
    : CandidateDBWriter(dataFileName, defaultIndexFileName(dataFileName), threads, mode) {}

CandidateDBWriter::CandidateDBWriter(
    const std::string &dataFileName,
    const std::string &indexFileName,
    unsigned int threads,
    size_t mode)
    : dataFileName(dataFileName),
      indexFileName(indexFileName),
      threads(threads == 0 ? 1 : threads),
      mode(mode) {}

CandidateDBWriter::~CandidateDBWriter() {
    close();
}

void CandidateDBWriter::open(size_t bufferSize) {
    writer.reset(new DBWriter(
        dataFileName.c_str(),
        indexFileName.c_str(),
        threads,
        mode,
        Parameters::DBTYPE_GENERIC_DB));
    writer->open(bufferSize);
}

void CandidateDBWriter::close(bool merge, bool needsSort) {
    if (writer && !writer->isClosed()) {
        writer->close(merge, needsSort);
    }
}

bool CandidateDBWriter::isClosed() const {
    return !writer || writer->isClosed();
}

void CandidateDBWriter::writeQuery(
    uint32_t queryId,
    const Query &query,
    unsigned int threadIdx)
{
    if (!writer || writer->isClosed()) {
        return;
    }

    std::string payload;
    serializeQuery(query, payload);
    writer->writeData(
        payload.data(),
        payload.size(),
        queryId,
        threadIdx,
        false,
        true);
}

std::string CandidateDBWriter::defaultIndexFileName(const std::string &dataFileName) {
    return dataFileName + ".index";
}

void CandidateDBWriter::serializeQuery(const Query &query, std::string &buffer) {
    buffer.clear();

    const uint32_t queryLength = static_cast<uint32_t>(query.queryLength + query.queryLength2);
    const uint32_t nameLength = static_cast<uint32_t>(query.name.size());
    const uint32_t candidateCount = static_cast<uint32_t>(query.speciesCandidates.size());

    appendPod(buffer, CANDIDATE_RECORD_MAGIC);
    appendPod(buffer, CANDIDATE_RECORD_VERSION);
    appendPod(buffer, queryLength);
    appendPod(buffer, nameLength);
    appendPod(buffer, candidateCount);
    buffer.append(query.name.data(), query.name.size());

    for (const SpeciesCandidate &candidate : query.speciesCandidates) {
        appendPod(buffer, candidate.speciesId);
        appendPod(buffer, candidate.idScore);
        appendPod(buffer, candidate.subScore);
        appendPod(buffer, candidate.logE);

        const uint32_t taxCountSize = static_cast<uint32_t>(candidate.taxCnt.size());
        appendPod(buffer, taxCountSize);
        for (const auto &taxCount : candidate.taxCnt) {
            appendPod(buffer, taxCount.first);
            appendPod(buffer, taxCount.second);
        }
    }
}
