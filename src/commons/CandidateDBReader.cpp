#include "CandidateDBReader.h"

#include <climits>
#include <cstring>
#include <utility>

namespace {
constexpr uint32_t CANDIDATE_RECORD_MAGIC = 0x444e4143; // CAND
constexpr uint32_t CANDIDATE_RECORD_VERSION = 1;

template <typename T>
bool readPod(const char *data, size_t dataSize, size_t &offset, T &value) {
    if (offset + sizeof(T) > dataSize) {
        return false;
    }
    std::memcpy(&value, data + offset, sizeof(T));
    offset += sizeof(T);
    return true;
}
}

CandidateDBReader::CandidateDBReader(
    const std::string &dataFileName,
    int threads)
    : CandidateDBReader(dataFileName, defaultIndexFileName(dataFileName), threads) {}

CandidateDBReader::CandidateDBReader(
    const std::string &dataFileName,
    const std::string &indexFileName,
    int threads)
    : dataFileName(dataFileName),
      indexFileName(indexFileName),
      threads(threads <= 0 ? 1 : threads) {}

CandidateDBReader::~CandidateDBReader() {
    close();
}

bool CandidateDBReader::open(int sortMode) {
    reader.reset(new DBReader<unsigned int>(
        dataFileName.c_str(),
        indexFileName.c_str(),
        threads,
        DBReader<unsigned int>::USE_INDEX | DBReader<unsigned int>::USE_DATA));
    reader->open(sortMode);
    return true;
}

void CandidateDBReader::close() {
    if (reader) {
        reader->close();
        reader.reset();
    }
}

size_t CandidateDBReader::size() const {
    return reader ? reader->getSize() : 0;
}

bool CandidateDBReader::getByIndex(size_t index, CandidateDBEntry &entry, int threadIdx) {
    if (!reader || index >= reader->getSize()) {
        return false;
    }

    char *data = reader->getData(index, threadIdx);
    const size_t dataSize = reader->getEntryLen(index);
    return deserializeEntry(reader->getDbKey(index), data, dataSize, entry);
}

bool CandidateDBReader::getByQueryId(uint32_t queryId, CandidateDBEntry &entry, int threadIdx) {
    if (!reader) {
        return false;
    }

    const size_t index = reader->getId(queryId);
    if (index == UINT_MAX) {
        return false;
    }

    return getByIndex(index, entry, threadIdx);
}

std::string CandidateDBReader::defaultIndexFileName(const std::string &dataFileName) {
    return dataFileName + ".index";
}

bool CandidateDBReader::deserializeEntry(
    uint32_t queryId,
    const char *data,
    size_t dataSize,
    CandidateDBEntry &entry)
{
    size_t offset = 0;
    uint32_t magic = 0;
    uint32_t version = 0;
    uint32_t queryLength = 0;
    uint32_t nameLength = 0;
    uint32_t candidateCount = 0;

    if (!readPod(data, dataSize, offset, magic)
        || !readPod(data, dataSize, offset, version)
        || !readPod(data, dataSize, offset, queryLength)
        || !readPod(data, dataSize, offset, nameLength)
        || !readPod(data, dataSize, offset, candidateCount)) {
        return false;
    }

    if (magic != CANDIDATE_RECORD_MAGIC || version != CANDIDATE_RECORD_VERSION) {
        return false;
    }

    if (offset + nameLength > dataSize) {
        return false;
    }

    CandidateDBEntry parsed;
    parsed.queryId = queryId;
    parsed.queryLength = queryLength;
    parsed.queryName.assign(data + offset, nameLength);
    offset += nameLength;
    parsed.candidates.reserve(candidateCount);

    for (uint32_t candidateIdx = 0; candidateIdx < candidateCount; ++candidateIdx) {
        SpeciesCandidate candidate;
        uint32_t taxCountSize = 0;
        if (!readPod(data, dataSize, offset, candidate.speciesId)
            || !readPod(data, dataSize, offset, candidate.idScore)
            || !readPod(data, dataSize, offset, candidate.subScore)
            || !readPod(data, dataSize, offset, candidate.logE)
            || !readPod(data, dataSize, offset, taxCountSize)) {
            return false;
        }

        candidate.taxCnt.reserve(taxCountSize);
        for (uint32_t taxCountIdx = 0; taxCountIdx < taxCountSize; ++taxCountIdx) {
            TaxID taxId = 0;
            uint32_t count = 0;
            if (!readPod(data, dataSize, offset, taxId)
                || !readPod(data, dataSize, offset, count)) {
                return false;
            }
            candidate.taxCnt.emplace_back(taxId, count);
        }
        parsed.candidates.push_back(std::move(candidate));
    }

    if (offset != dataSize) {
        return false;
    }

    entry = std::move(parsed);
    return true;
}
