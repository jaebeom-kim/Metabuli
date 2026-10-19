#ifndef METABULI_DELTAIDXREADER_H
#define METABULI_DELTAIDXREADER_H

#include <cstdint>
#include <cstddef>
#include <vector>
#include <string>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <memory>
#include <algorithm>
#include <unistd.h>
#include <fcntl.h>

#include "Kmer.h"
#include "common.h"
#include "InfoIndex.h"

#define MEM_SIZE_16MB ((size_t) (16 * 1024 * 1024))
#define MEM_SIZE_32MB ((size_t) (32 * 1024 * 1024))

class KmerDbReader {
private:
    std::string fileName;

    size_t valueBufferSize;
    uint64_t * valueBuffer;
    size_t valueBufferIdx = 0;
    size_t valueCnt;
    uint64_t lastValue;

    ReadBuffer<uint16_t> deltaIdxBuffer;
    bool fileCompleted = false;
    bool valueBufferCompleted = false;

    void fillValueBuffer() {
        for (; valueCnt < valueBufferSize; ++valueCnt) {
            valueBuffer[valueCnt] = getNextMetamer();
            if (unlikely(valueBuffer[valueCnt] == UINT64_MAX)) {
                fileCompleted = true;
                break;
            }
        }
    }

    uint64_t getNextMetamer() {
        if (deltaIdxBuffer.end < deltaIdxBuffer.p + 7) {
            size_t readCnt = deltaIdxBuffer.loadBuffer(deltaIdxBuffer.end - deltaIdxBuffer.p);
            if (readCnt == 0) {
                return UINT64_MAX; // No more values
            }
        }
        uint64_t diffIn64bit = 0;
        while ((*deltaIdxBuffer.p & 0x8000) == 0) {
            diffIn64bit = (diffIn64bit << 15) | *deltaIdxBuffer.p;
            ++deltaIdxBuffer.p;
        }
        diffIn64bit = (diffIn64bit << 15) | (*deltaIdxBuffer.p & 0x7FFF);
        ++deltaIdxBuffer.p;
        this->lastValue = diffIn64bit + this->lastValue;
        return this->lastValue;
    }

public:
    KmerDbReader(
        std::string fileName,
        size_t valueBufferSize = 32768, 
        size_t readBufferSize = 8192)
        : fileName(fileName),
        valueBufferSize(valueBufferSize), 
        deltaIdxBuffer(fileName, readBufferSize)
    {
        lastValue = 0;
        valueCnt = 0;
        valueBuffer = new uint64_t[valueBufferSize];
        fillValueBuffer();
    }

    ~KmerDbReader() {
        delete[] valueBuffer;
    }

    uint64_t getLastValue() const {
        return lastValue;
    }

    bool isCompleted() const {
        return fileCompleted && valueBufferCompleted;
    }

    uint64_t next() {
        if (unlikely(valueBufferIdx >= valueCnt)) {
            valueCnt = 0;
            valueBufferIdx = 0;
            fillValueBuffer();
        }
        if (unlikely(valueCnt == 0)) {
            valueBufferCompleted = true;
            return uint64_t(); // Return an empty k-mer
        }
        return valueBuffer[valueBufferIdx++];
    }

    void setReadPosition(DiffIdxSplit offset) {
        deltaIdxBuffer.loadBufferAt(offset.diffIdxOffset);
        if (offset.ADkmer == 0 && offset.diffIdxOffset == 0 && offset.infoIdxOffset == 0) {
            valueCnt = 0;
            lastValue = 0;
        } else {
            lastValue = offset.ADkmer;
            valueBuffer[0] = lastValue;
            valueCnt = 1;
        }
        valueBufferIdx = 0;
        fillValueBuffer();
    }

    uint64_t * getValueBuffer() {
        return valueBuffer;
    }

    size_t getValueCnt() const {
        return valueCnt;
    }

};
class DeltaIdxReader {
private:
    std::string deltaIdxFileName;
    std::string infoFileName;
    std::string posFileName;
    size_t totalValueNum;
    bool hasPos = false;

    // To manage values
    size_t valueBufferSize;
    Kmer * valueBuffer;
    size_t valueBufferIdx = 0;
    size_t valueCnt;
    uint64_t lastValue;
    
    // To manage delta indices    
    ReadBuffer<uint16_t> deltaIdxBuffer;
    InfoIndexMetadata infoMetadata;
    bool infoIsPacked = false;
    // Only one info buffer is active. Plain DBs keep the old direct TaxID
    // pointer path; packed DBs decode from uint64 words.
    std::unique_ptr<ReadBuffer<TaxID>> plainInfoBuffer;
    std::unique_ptr<ReadBuffer<uint64_t>> packedInfoBuffer;
    uint8_t packedIdBits = 32;
    uint8_t packedValuesPerWord = 1;
    uint64_t packedIdMask = UINT32_MAX;
    uint64_t packedCurrentWord = 0;
    uint8_t packedLane = 0;
    size_t packedReadCount = 0;
    // Optional per-k-mer positions, read in lockstep with the info index
    // (whether that index is plain or packed). Populated only via the
    // posFileName constructor.
    ReadBuffer<uint16_t> posBuffer;
    bool fileCompleted = false;
    bool valueBufferCompleted = false;

    // --- Lazy info mode (classify-time matchKmers) --------------------------
    // The metamer stream must be decoded fully (deltas are sequential), but
    // tax IDs are only needed for AA-matched candidate groups, which are
    // contiguous ordinal ranges of the sorted k-mer list. In lazy mode next()
    // skips the info stream entirely; fetchInfos() reads a requested range
    // through a small seekable window, so info bytes of never-matched regions
    // are never read from disk.
    bool lazyInfo = false;
    size_t decodeOrdinal = 0;      // ordinal of the next metamer to decode
    size_t bufferBaseOrdinal = 0;  // ordinal of valueBuffer[0]
    std::unique_ptr<ReadBuffer<TaxID>> fetchPlainBuffer;
    std::unique_ptr<ReadBuffer<uint64_t>> fetchWordBuffer;
    size_t fetchCursorOrdinal = SIZE_MAX; // ordinal after the last fetched ID
    uint64_t fetchWord = 0;
    uint8_t fetchLane = 0;
    size_t fetchBufStartElem = SIZE_MAX;  // file element index of the window start
    size_t fetchedIdCnt = 0;
    size_t fetchSeekCnt = 0;
    size_t fetchCallCnt = 0;
public:
    // log2-bucket histogram of candidate group sizes:
    // [1, 2-3, 4-7, 8-15, 16-31, 32-63, 64-127, 128-255, 256-511, 512+]
    static constexpr int GROUP_HIST_BUCKETS = 10;
private:
    size_t fetchGroupCnt[GROUP_HIST_BUCKETS] = {0};  // groups per bucket
    size_t fetchGroupIds[GROUP_HIST_BUCKETS] = {0};  // candidates per bucket
    static constexpr size_t FETCH_SEEK_READ = 512; // elements per cold seek (4 KB)

    void fillValueBufferLazy() {
        if (valueCnt == 0) {
            bufferBaseOrdinal = decodeOrdinal;
        }
        for (; valueCnt < valueBufferSize && decodeOrdinal < totalValueNum; ++valueCnt) {
            const uint64_t v = getNextMetamer();
            if (unlikely(v == UINT64_MAX)) {
                fileCompleted = true;
                break;
            }
            valueBuffer[valueCnt].value = v;
            valueBuffer[valueCnt].id = 0; // filled on demand via fetchInfos()
            ++decodeOrdinal;
        }
        if (decodeOrdinal >= totalValueNum) {
            fileCompleted = true;
        }
    }

    // Reads one position in lockstep with the info stream, refilling as needed.
    // Positions are 1:1 with info IDs regardless of whether the info is packed.
    inline uint16_t nextPos() {
        if (unlikely(posBuffer.p == posBuffer.end)) {
            posBuffer.loadBuffer();
        }
        return *posBuffer.p++;
    }

    void fillValueBuffer() {
        if (lazyInfo) {
            fillValueBufferLazy();
            return;
        }
        // Format is chosen once per refill, not once per ID. This keeps legacy
        // uint32 databases on the same hot path they used before bitpacking.
        if (unlikely(infoIsPacked)) {
            switch (packedIdBits) {
                case 8:  fillValueBufferPacked<8>(); break;
                case 10: fillValueBufferPacked<10>(); break;
                case 12: fillValueBufferPacked<12>(); break;
                case 16: fillValueBufferPacked<16>(); break;
                case 21: fillValueBufferPacked<21>(); break;
                default: fillValueBufferPackedGeneric(); break;
            }
            return;
        }
        fillValueBufferPlain();
    }

    void fillValueBufferPlain() {
        for (; valueCnt < valueBufferSize; ++valueCnt) {
            if (unlikely(plainInfoBuffer->p == plainInfoBuffer->end)) {
                if (hasPos && posBuffer.p != posBuffer.end) {
                    std::cerr << "Error: Info buffer is exhausted but position buffer is not. This should not happen." << std::endl;
                    exit(EXIT_FAILURE);
                }
                size_t readCnt = plainInfoBuffer->loadBuffer();
                if (readCnt == 0) {
                    fileCompleted = true;
                    break;
                }
            }

            valueBuffer[valueCnt].tInfo.taxId = *plainInfoBuffer->p++;
            valueBuffer[valueCnt].value = getNextMetamer();

            if (hasPos) {
                valueBuffer[valueCnt].tInfo.pos = nextPos();
            }
        }
    }

    template <uint8_t ID_BITS>
    void fillValueBufferPacked() {
        static constexpr uint8_t VALUES_PER_WORD = InfoIndex::WORD_BITS / ID_BITS;
        static constexpr uint64_t ID_MASK = (uint64_t{1} << ID_BITS) - 1;

        while (valueCnt < valueBufferSize && packedReadCount < totalValueNum) {
            if (unlikely(packedLane >= VALUES_PER_WORD)) {
                if (unlikely(!loadPackedWord())) {
                    break;
                }
            }

            // Decode all remaining lanes from the current word before touching
            // the packed input buffer again. ID_BITS is a compile-time constant
            // in the common formats, so shifts and masks optimize well.
            while (valueCnt < valueBufferSize &&
                   packedLane < VALUES_PER_WORD &&
                   packedReadCount < totalValueNum) {
                valueBuffer[valueCnt].tInfo.taxId =
                    static_cast<TaxID>((packedCurrentWord >> (packedLane * ID_BITS)) & ID_MASK);
                valueBuffer[valueCnt].value = getNextMetamer();
                if (hasPos) {
                    valueBuffer[valueCnt].tInfo.pos = nextPos();
                }
                ++valueCnt;
                ++packedLane;
                ++packedReadCount;
            }
        }
        if (unlikely(packedReadCount >= totalValueNum)) {
            fileCompleted = true;
        }
    }

    void fillValueBufferPackedGeneric() {
        while (valueCnt < valueBufferSize && packedReadCount < totalValueNum) {
            if (unlikely(packedLane >= packedValuesPerWord)) {
                if (unlikely(!loadPackedWord())) {
                    break;
                }
            }
            while (valueCnt < valueBufferSize &&
                   packedLane < packedValuesPerWord &&
                   packedReadCount < totalValueNum) {
                valueBuffer[valueCnt].tInfo.taxId =
                    static_cast<TaxID>((packedCurrentWord >> (packedLane * packedIdBits)) & packedIdMask);
                valueBuffer[valueCnt].value = getNextMetamer();
                if (hasPos) {
                    valueBuffer[valueCnt].tInfo.pos = nextPos();
                }
                ++valueCnt;
                ++packedLane;
                ++packedReadCount;
            }
        }
        if (unlikely(packedReadCount >= totalValueNum)) {
            fileCompleted = true;
        }
    }

    inline bool loadPackedWord() {
        if (unlikely(packedInfoBuffer->p >= packedInfoBuffer->end)) {
            if (packedInfoBuffer->loadBuffer() == 0) {
                fileCompleted = true;
                return false;
            }
        }
        packedCurrentWord = *packedInfoBuffer->p++;
        packedLane = 0;
        return true;
    }

    void loadPackedInfoAt(size_t logicalOffset) {
        packedReadCount = logicalOffset;
        const size_t wordOffset = logicalOffset / packedValuesPerWord;
        const uint8_t targetLane = logicalOffset % packedValuesPerWord;
        packedInfoBuffer->loadBufferAt(wordOffset);

        // Starting exactly on a word boundary lets the normal refill path load
        // the word. Starting inside a word needs that word cached immediately.
        packedCurrentWord = 0;
        packedLane = packedValuesPerWord;
        if (targetLane != 0) {
            loadPackedWord();
            packedLane = targetLane;
        }
    }

    inline bool readNextPackedIdGeneric(uint32_t &taxId) {
        if (unlikely(packedReadCount >= totalValueNum)) {
            return false;
        }
        if (unlikely(packedLane >= packedValuesPerWord && !loadPackedWord())) {
            return false;
        }
        taxId = static_cast<uint32_t>((packedCurrentWord >> (packedLane * packedIdBits)) & packedIdMask);
        ++packedLane;
        ++packedReadCount;
        return true;
    }

    template <uint8_t ID_BITS>
    inline bool readNextPackedId(uint32_t &taxId) {
        static constexpr uint8_t VALUES_PER_WORD = InfoIndex::WORD_BITS / ID_BITS;
        static constexpr uint64_t ID_MASK = (uint64_t{1} << ID_BITS) - 1;
        if (unlikely(packedReadCount >= totalValueNum)) {
            return false;
        }
        if (unlikely(packedLane >= VALUES_PER_WORD && !loadPackedWord())) {
            return false;
        }
        taxId = static_cast<uint32_t>((packedCurrentWord >> (packedLane * ID_BITS)) & ID_MASK);
        ++packedLane;
        ++packedReadCount;
        return true;
    }

    bool readNextPackedId(uint32_t &taxId) {
        // This is used only after a split seek to seed valueBuffer[0]. Bulk
        // streaming uses fillValueBufferPacked() above.
        switch (packedIdBits) {
            case 8:  return readNextPackedId<8>(taxId);
            case 10: return readNextPackedId<10>(taxId);
            case 12: return readNextPackedId<12>(taxId);
            case 16: return readNextPackedId<16>(taxId);
            case 21: return readNextPackedId<21>(taxId);
            default: return readNextPackedIdGeneric(taxId);
        }
    }

    uint64_t getNextMetamer() {
        if (unlikely(deltaIdxBuffer.end - deltaIdxBuffer.p < 5)) {
            size_t readCnt = deltaIdxBuffer.loadBuffer(deltaIdxBuffer.end - deltaIdxBuffer.p);
            if (readCnt == 0) {
                return UINT64_MAX; // No more values
            }
        }

        uint16_t* p = deltaIdxBuffer.p;
        uint64_t diffIn64bit = 0;
        uint16_t word;

        do {
            word = *p++;
            diffIn64bit = (diffIn64bit << 15) | (word & 0x7FFF);
        } while ((word & 0x8000) == 0);

        deltaIdxBuffer.p = p;
        lastValue += diffIn64bit;
        return lastValue;
    }

public:
    DeltaIdxReader(
        std::string deltaIdxFileName,
        std::string infoFileName,
        size_t valueBufferSize = 32768,
        size_t readBufferSize = 8192,
        bool lazyInfo = false)
        : deltaIdxFileName(deltaIdxFileName),
        infoFileName(infoFileName),
        valueBufferSize(valueBufferSize),
        deltaIdxBuffer(deltaIdxFileName, readBufferSize),
        infoMetadata(InfoIndex::loadMetadata(infoFileName)),
        lazyInfo(lazyInfo)
    {
        lastValue = 0;
        valueCnt = 0;
        valueBuffer = new Kmer[valueBufferSize];
        infoIsPacked = infoMetadata.isPacked();
        // The on-demand fetch window is deliberately small: a seek costs one
        // window read, so dense regions stream as before while sparse regions
        // skip whole stretches of the info file.
        constexpr size_t FETCH_WINDOW_WORDS = 32768;  // 256 KB
        constexpr size_t FETCH_WINDOW_IDS = 65536;    // 256 KB
        if (infoIsPacked) {
            packedIdBits = infoMetadata.idBits;
            packedValuesPerWord = InfoIndex::idsPerWord(packedIdBits);
            packedIdMask = InfoIndex::maskForBits(packedIdBits);
            totalValueNum = infoMetadata.idCount;
            if (lazyInfo) {
                fetchWordBuffer = std::make_unique<ReadBuffer<uint64_t>>(infoFileName,
                                                                         FETCH_WINDOW_WORDS);
            } else {
                // readBufferSize is historically a count of 32-bit IDs. A half-size
                // uint64 buffer keeps roughly the same byte budget for packed data.
                const size_t wordBufferSize = std::max<size_t>(1, readBufferSize / 2);
                packedInfoBuffer = std::make_unique<ReadBuffer<uint64_t>>(infoFileName,
                                                                          wordBufferSize);
                packedLane = packedValuesPerWord;
            }
        } else {
            totalValueNum = FileUtil::getFileSize(infoFileName) / sizeof(TaxID);
            if (lazyInfo) {
                fetchPlainBuffer = std::make_unique<ReadBuffer<TaxID>>(infoFileName,
                                                                      FETCH_WINDOW_IDS);
            } else {
                plainInfoBuffer = std::make_unique<ReadBuffer<TaxID>>(infoFileName,
                                                                      readBufferSize);
            }
        }
        fillValueBuffer();
    }

    // Position-index variant: a plain (unpacked) info index read in lockstep
    // with a separate per-k-mer position file.
    DeltaIdxReader(
        std::string deltaIdxFileName,
        std::string infoFileName,
        std::string posFileName,
        size_t valueBufferSize = 32768,
        size_t readBufferSize = 8192)
        : deltaIdxFileName(deltaIdxFileName),
        infoFileName(infoFileName),
        posFileName(posFileName),
        valueBufferSize(valueBufferSize),
        deltaIdxBuffer(deltaIdxFileName, readBufferSize),
        infoMetadata(InfoIndex::loadMetadata(infoFileName)),
        posBuffer(posFileName, readBufferSize)
    {
        lastValue = 0;
        valueCnt = 0;
        valueBuffer = new Kmer[valueBufferSize];
        hasPos = true;
        // Position DBs pack their info index like any other DB; positions live
        // in a separate file read in lockstep with the (plain or packed) IDs.
        infoIsPacked = infoMetadata.isPacked();
        if (infoIsPacked) {
            packedIdBits = infoMetadata.idBits;
            packedValuesPerWord = InfoIndex::idsPerWord(packedIdBits);
            packedIdMask = InfoIndex::maskForBits(packedIdBits);
            const size_t wordBufferSize = std::max<size_t>(1, readBufferSize / 2);
            packedInfoBuffer = std::make_unique<ReadBuffer<uint64_t>>(infoFileName,
                                                                      wordBufferSize);
            packedLane = packedValuesPerWord;
            totalValueNum = infoMetadata.idCount;
        } else {
            plainInfoBuffer = std::make_unique<ReadBuffer<TaxID>>(infoFileName,
                                                                  readBufferSize);
            totalValueNum = FileUtil::getFileSize(infoFileName) / sizeof(TaxID);
        }
        fillValueBuffer();
    }

    ~DeltaIdxReader() {
        delete[] valueBuffer;
    }

    // void setPosBuffer(const std::string & posFileName) {
    //     this->posFileName = posFileName;
    //     posBuffer = ReadBuffer<uint16_t>(posFileName, readBufferSize);
    //     hasPos = true;
    // }


    uint64_t getLastValue() const {
        return lastValue;
    }

    // Copy values <= maxValue to the provided buffer
    size_t getValues(Kmer * largeBuffer, uint64_t maxValue) {
        size_t n = 0;
        while (n < valueCnt && valueBuffer[n].value <= maxValue) {
            ++n;
        }
        if (n > 0) {
            std::memcpy(largeBuffer, valueBuffer, n * sizeof(Kmer));
            std::memmove(valueBuffer, valueBuffer + n, (valueCnt - n) * sizeof(Kmer));
            valueCnt -= n;
            if (!fileCompleted) fillValueBuffer();
            if (valueCnt == 0) {
                valueBufferCompleted = true;
            }
        }
        return n;
    }

    bool isCompleted() const {
        return fileCompleted && valueBufferCompleted;
    }

    size_t getTotalValueNum() const {
        return totalValueNum;
    }

    Kmer next() {
        if (unlikely(valueBufferIdx >= valueCnt)) {
            valueCnt = 0;
            valueBufferIdx = 0;
            fillValueBuffer();
        }
        if (unlikely(valueCnt == 0)) {
            valueBufferCompleted = true;
            return Kmer(); // Return an empty k-mer
        }
        return valueBuffer[valueBufferIdx++];
    }

    Kmer current() {
        if (unlikely(valueBufferIdx >= valueCnt)) {
            valueCnt = 0;
            valueBufferIdx = 0;
            fillValueBuffer();
        }
        if (unlikely(valueCnt == 0)) {
            valueBufferCompleted = true;
            return {UINT64_MAX, UINT32_MAX}; // Return a dummy k-mer
        }
        return valueBuffer[valueBufferIdx];
    }

    void setReadPosition(DiffIdxSplit offset) {
        // The same reader instance can be reused for multiple query splits.
        // Seeking must clear EOF state from any previous scan.
        fileCompleted = false;
        valueBufferCompleted = false;
        deltaIdxBuffer.loadBufferAt(offset.diffIdxOffset);
        const size_t infoOffset = offset.infoIdxOffset - (offset.ADkmer != 0);
        if (lazyInfo) {
            // Only ordinal bookkeeping; the info file is touched by
            // fetchInfos() alone.
            if (offset.ADkmer == 0 && offset.diffIdxOffset == 0 && offset.infoIdxOffset == 0) {
                valueCnt = 0;
                lastValue = 0;
                decodeOrdinal = 0;
                bufferBaseOrdinal = 0;
            } else {
                lastValue = offset.ADkmer;
                valueBuffer[0].value = lastValue;
                valueBuffer[0].id = 0;
                valueCnt = 1;
                bufferBaseOrdinal = infoOffset;  // the anchor k-mer's ordinal
                decodeOrdinal = infoOffset + 1;
            }
            valueBufferIdx = 0;
            fillValueBuffer();
            return;
        }
        if (infoIsPacked) {
            // Split offsets are logical ID counts. Packed files translate that
            // once at seek time, then stream IDs from the chosen word/lane.
            loadPackedInfoAt(infoOffset);
        } else {
            plainInfoBuffer->loadBufferAt(infoOffset);
        }
        if (hasPos) {
            // Positions are 1:1 with info entries, so the same logical offset applies.
            posBuffer.loadBufferAt(infoOffset);
        }
        if (offset.ADkmer == 0 && offset.diffIdxOffset == 0 && offset.infoIdxOffset == 0) {
            valueCnt = 0;
            lastValue = 0;
        } else {
            lastValue = offset.ADkmer;
            valueBuffer[0].value = lastValue;
            if (infoIsPacked) {
                uint32_t taxId = 0;
                readNextPackedId(taxId);
                valueBuffer[0].tInfo.taxId = static_cast<TaxID>(taxId);
            } else {
                valueBuffer[0].tInfo.taxId = *plainInfoBuffer->p++;
            }
            if (hasPos) {
                valueBuffer[0].tInfo.pos = *posBuffer.p++;
            }
            valueCnt = 1;
        }
        valueBufferIdx = 0;
        fillValueBuffer();
    }

    // Ordinal of the k-mer returned by the most recent next() (lazy mode).
    size_t lastOrdinal() const {
        return bufferBaseOrdinal + valueBufferIdx - 1;
    }

    // Lazy mode: read n tax IDs starting at startOrdinal. Candidate groups are
    // contiguous and visited in ascending order, so consecutive calls either
    // continue the window (dense regions, streamed as before) or seek forward
    // (sparse regions, skipping unread info bytes).
    void fetchInfos(size_t startOrdinal, size_t n, uint32_t *out) {
        fetchedIdCnt += n;
        ++fetchCallCnt;
        const int bucket = std::min(GROUP_HIST_BUCKETS - 1,
                                    63 - __builtin_clzll(static_cast<unsigned long long>(n)));
        ++fetchGroupCnt[bucket];
        fetchGroupIds[bucket] += n;
        if (infoIsPacked) {
            if (startOrdinal != fetchCursorOrdinal) {
                const size_t wordOffset = startOrdinal / packedValuesPerWord;
                const uint8_t targetLane = startOrdinal % packedValuesPerWord;
                const size_t bufWords = fetchWordBuffer->end - fetchWordBuffer->start;
                if (fetchBufStartElem != SIZE_MAX
                    && wordOffset >= fetchBufStartElem
                    && wordOffset < fetchBufStartElem + bufWords) {
                    // Target is already buffered: reposition without a syscall.
                    fetchWordBuffer->p = fetchWordBuffer->start + (wordOffset - fetchBufStartElem);
                } else {
                    // Cold seek: read only a small chunk; a dense region keeps
                    // extending via sequential refills below.
                    ++fetchSeekCnt;
                    const size_t lastWord = (startOrdinal + n - 1) / packedValuesPerWord;
                    const size_t wordsNeeded = lastWord - wordOffset + 1;
                    fetchWordBuffer->loadBufferAt(wordOffset,
                                                  std::max(wordsNeeded, FETCH_SEEK_READ));
                    fetchBufStartElem = wordOffset;
                }
                fetchLane = packedValuesPerWord; // force a word load
                if (targetLane != 0) {
                    fetchWord = *fetchWordBuffer->p++;
                    fetchLane = targetLane;
                }
            }
            for (size_t k = 0; k < n; ++k) {
                if (fetchLane >= packedValuesPerWord) {
                    if (unlikely(fetchWordBuffer->p >= fetchWordBuffer->end)) {
                        fetchBufStartElem += fetchWordBuffer->end - fetchWordBuffer->start;
                        fetchWordBuffer->loadBuffer();
                    }
                    fetchWord = *fetchWordBuffer->p++;
                    fetchLane = 0;
                }
                out[k] = static_cast<uint32_t>(
                    (fetchWord >> (fetchLane * packedIdBits)) & packedIdMask);
                ++fetchLane;
            }
        } else {
            if (startOrdinal != fetchCursorOrdinal) {
                const size_t bufIds = fetchPlainBuffer->end - fetchPlainBuffer->start;
                if (fetchBufStartElem != SIZE_MAX
                    && startOrdinal >= fetchBufStartElem
                    && startOrdinal < fetchBufStartElem + bufIds) {
                    fetchPlainBuffer->p = fetchPlainBuffer->start + (startOrdinal - fetchBufStartElem);
                } else {
                    ++fetchSeekCnt;
                    fetchPlainBuffer->loadBufferAt(startOrdinal, std::max(n, FETCH_SEEK_READ));
                    fetchBufStartElem = startOrdinal;
                }
            }
            for (size_t k = 0; k < n; ++k) {
                if (unlikely(fetchPlainBuffer->p >= fetchPlainBuffer->end)) {
                    fetchBufStartElem += fetchPlainBuffer->end - fetchPlainBuffer->start;
                    fetchPlainBuffer->loadBuffer();
                }
                out[k] = static_cast<uint32_t>(*fetchPlainBuffer->p++);
            }
        }
        fetchCursorOrdinal = startOrdinal + n;
    }

    size_t getFetchedIdCnt() const { return fetchedIdCnt; }
    size_t getFetchSeekCnt() const { return fetchSeekCnt; }
    size_t getFetchCallCnt() const { return fetchCallCnt; }
    const size_t * getFetchGroupCnt() const { return fetchGroupCnt; }
    const size_t * getFetchGroupIds() const { return fetchGroupIds; }

    Kmer  * getValueBuffer() {
        return valueBuffer;
    }

    size_t getValueCnt() const {
        return valueCnt;
    }
};
#endif // METABULI_DELTAIDXREADER_H
