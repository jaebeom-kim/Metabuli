#ifndef ADCLASSIFIER2_COMMON_H
#define ADCLASSIFIER2_COMMON_H
#include <cstddef>
#include <utility>

#include <iostream>
#include <unordered_set>
#include "FileUtil.h"
// #include "Match.h"
#include <cstdint>
#include <bitset>

#include <queue>
#include <mutex>
#include <condition_variable>
#include <cmath>
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <cstring>
#include <thread>

#include "LocalParameters.h"
#include "TaxonomyWrapper.h"
#include "printBinary.h"


#define likely(x) __builtin_expect((x),1)
#define unlikely(x) __builtin_expect((x),0)
#define AA(kmer) ((kmer) & ~16777215)

extern const std::string atcg;
extern const std::string iRCT;

struct CovMetric {
    double evenness;
    double coverage;
    double adjustedEvenness;
    double unifiedScore;
    double macroCoverage;
};

// Compute genome-coverage metrics for one species from its per-bin hit counts.
// `bins` is a genome-position histogram (posId -> count), `readCnt` the number of
// reads contributing to it, `totalReadLength` their summed length, `genomeSize`
// the species genome size. Shared by the classify and filter-candidates paths.
CovMetric computeCoverageMetric(
    const std::vector<uint8_t> &bins,
    uint64_t readCnt,
    uint64_t totalReadLength,
    uint64_t genomeSize);


struct MappingRes {
    MappingRes(uint32_t queryId, TaxID speicesId, float score) 
        : queryId(queryId), speciesId(speicesId), score(score) {}
    MappingRes() : queryId(0), speciesId(0), score(0.0f) {}
    uint32_t queryId;
    TaxID speciesId;
    float score; 
};

struct SpeciesCandidate {
    TaxID speciesId = 0;
    float idScore = 0.0f;
    float subScore = 0.0f;
    float logE = 0.0f;
    std::vector<std::pair<TaxID, uint32_t>> taxCnt;
    // Unique genome bins (k-mer positions, posId) this read touched for this
    // species. Populated only when the DB stores k-mer positions. Used later to
    // estimate per-species genome coverage from a candidate DB.
    std::vector<uint16_t> kmerPositions;
};

// Tie-margin factor shared by classify-candidates (chooseBestTaxonFromCandidates)
// and filter-candidates. A candidate ties with the top hit when its score is
// strictly greater than bestScore * this factor. The factor scales with the best
// idScore so higher-confidence reads use a tighter margin. Every tie decision in
// both commands must go through this so "survives filtering" and "resolves to one
// species at classification" use the exact same rule.
inline float tieMarginFactor(float bestIdScore, float tieRatio) {
    const float diff = 0.09f;
    return (tieRatio - diff) + (bestIdScore * diff);
}

struct Assembly {
    std::string name;
    TaxID taxid;
    TaxID speciesId;
    TaxID genusId;
    TaxID familyId;
    TaxID orderId;
    Assembly(std::string name) : name(name) {}
    Assembly() : name(""), taxid(0), speciesId(0), genusId(0), familyId(0), orderId(0) {}
    

    void print () const {
        std::cout << "Assembly: " << name << ", TaxID: " << taxid
                  << ", SpeciesID: " << speciesId
                  << ", GenusID: " << genusId
                  << ", FamilyID: " << familyId
                  << ", OrderID: " << orderId << std::endl;
    }
};

struct KmerCnt {
    KmerCnt(size_t length, size_t kmerCnt, size_t totalCnt) : length(length), kmerCnt(kmerCnt), totalCnt(totalCnt) {}
    KmerCnt() : length(0), kmerCnt(0), totalCnt(0){}
    size_t length;
    size_t kmerCnt;
    size_t totalCnt;
};

struct CDSinfo{
    uint32_t protId; //4,294,967,295 counted from 0
    int frame;
    bool isComplement;
    std::vector<std::pair<size_t, size_t>> loc;
    CDSinfo() = default;
    CDSinfo(uint32_t protId, int frame) : protId(protId), frame(frame) {}
};


struct SequenceBlock {
    SequenceBlock(int start, int end, int strand) : start(start), end(end), strand(strand) {}

    void printSequenceBlock() {
        std::cout << strand << " " << start << " " << end << std::endl;
    }

    int start;  // 0-based inclusive
    int end;    // 0-based inclusive
    int strand; // 1: forward, -1: reverse
};

struct Classification {
    Classification() : taxId(0), length(0), score(0.0) {}
    Classification(const std::string & name, TaxID taxId, int length, double score)
        : name(std::move(name)), taxId(taxId), length(length), score(score) {}
    Classification(const std::string & name, int length) 
        : name(std::move(name)), taxId(0), length(length), score(0.0) {}
    std::string name;
    TaxID taxId;
    int length;
    double score;
};




struct Query {
    int classification;
    float idScore;
    float subScore;
    float eValue;
    int hammingDist;
    int queryLength;
    int queryLength2;
    int kmerCnt;
    int kmerCnt2;
    TaxID topSpeciesId; 

    std::string name;
    std::map<TaxID,int> taxCnt; // 8 byte per element
    std::vector<std::pair<TaxID, float>> species2Score;
    std::vector<SpeciesCandidate> speciesCandidates;

    Query()
        : classification(0), idScore(0), subScore(0), eValue(-1), hammingDist(0),
          queryLength(0), queryLength2(0), kmerCnt(0), kmerCnt2(0), topSpeciesId(0) {}
    Query(int queryLength, int kmerCnt, const std::string & name)
        : classification(0), idScore(0), subScore(0), eValue(-1), hammingDist(0),
          queryLength(queryLength), queryLength2(0), kmerCnt(kmerCnt), kmerCnt2(0),
          topSpeciesId(0), name(name) {}
};

struct ProteinQuery {
    uint32_t queryId;
    uint32_t length;
    uint32_t kmerCnt;
    uint32_t result;
    uint32_t kmerMatchCnt;
    float score;

    std::string name;

    ProteinQuery(uint32_t queryId, uint32_t length, uint32_t kmerCnt, uint32_t result, float score, std::string name)
        : queryId(queryId), length(length), kmerCnt(kmerCnt), result(result), score(score), name(std::move(name)) {}
    
    ProteinQuery() : queryId(0), length(0), kmerCnt(0), result(0), kmerMatchCnt(0), score(0.0) {}

};

template<typename T>
struct Buffer {
    T *buffer;
    size_t startIndexOfReserve;
    size_t bufferSize;

    explicit Buffer(size_t sizeOfBuffer=100) {
        buffer = (T *) calloc(sizeOfBuffer, sizeof(T));
        bufferSize = sizeOfBuffer;
        startIndexOfReserve = 0;
    };

    ~Buffer() {
        if (buffer) {
            free(buffer);
        }
    };

    size_t reserveMemory(size_t numOfKmer) {
        return __sync_fetch_and_add(&startIndexOfReserve, numOfKmer);
    };

    void reallocateMemory(size_t sizeOfBuffer) {
        buffer = (T *) realloc(buffer, sizeof(T) * sizeOfBuffer);
        bufferSize = sizeOfBuffer;
    };

    void init() {
        startIndexOfReserve = 0;
        if (buffer) {
            memset(buffer, 0, sizeof(T) * bufferSize);
        }
    }

    static size_t calculateBufferSize(size_t maxRam, size_t threads, size_t bytePerKmer) {
        constexpr double GIGABYTE = 1024.0 * 1024.0 * 1024.0;
        constexpr double MEGABYTE = 1024.0 * 1024.0;
        constexpr double MEMORY_PER_THREAD_MB = 50.0;
        
        float c = 0.7f;
        if (maxRam < 16) {
            c = 0.5f;
        } else if (maxRam <= 32) {
            c = 0.6f;
        }
        double totalRamBytes = maxRam * GIGABYTE * c;
        double memoryForThreads = threads * MEMORY_PER_THREAD_MB * MEGABYTE;
        double availableMemory = totalRamBytes - memoryForThreads;

        if (availableMemory <= 0.0) {
            std::cerr << "Not enough memory to create index" << std::endl;
            std::cerr << "Please increase the RAM usage or decrease the number of threads" << std::endl;
            exit(EXIT_FAILURE);
        }

        return static_cast<size_t>(availableMemory / bytePerKmer);
    }

    static bool moveSmallToLarge(Buffer<T> * small, Buffer<T> * large) {
        size_t posToWrite = large->reserveMemory(small->startIndexOfReserve);
        if (unlikely(posToWrite + small->startIndexOfReserve >= large->bufferSize)) {
            __sync_fetch_and_sub(&large->startIndexOfReserve, small->startIndexOfReserve);
            return false;
        }
        memcpy(large->buffer + posToWrite, small->buffer, sizeof(T) * small->startIndexOfReserve);
        small->startIndexOfReserve = 0;
        return true;
    } 

    bool afford(size_t num) {
        return (startIndexOfReserve + num < bufferSize);
    }
};

// --- Optional I/O profiling (METABULI_IO_PROF=1) -----------------------------
// Accumulates wall time spent blocked in ReadBuffer freads and bytes read,
// summed across threads. KmerMatcher prints per-search deltas when the env
// var is set; otherwise the only overhead is a cached-bool check per refill.
inline bool ioProfEnabled() {
    static const bool enabled = getenv("METABULI_IO_PROF") != nullptr;
    return enabled;
}
inline std::atomic<uint64_t> & ioProfNanos() {
    static std::atomic<uint64_t> v{0};
    return v;
}
inline std::atomic<uint64_t> & ioProfBytes() {
    static std::atomic<uint64_t> v{0};
    return v;
}
// Chunks fetched ahead by prefetch workers; shows whether adaptive read-ahead
// decided to overlap (printed per search with METABULI_IO_PROF=1).
inline std::atomic<uint64_t> & ioPrefetchChunks() {
    static std::atomic<uint64_t> v{0};
    return v;
}
// --- Adaptive read-ahead (METABULI_PREFETCH) ---------------------------------
// Eligible ReadBuffers can switch to a double-buffered mode: while the
// consumer drains one half, a background worker freads the next chunk into
// the other half, overlapping disk reads with decoding/compute. Overlap only
// pays when reads actually block on the device; on page-cache hits it costs
// handshakes, cross-core cache transfers, and extra memcpy CPU. By default
// each stream therefore decides for itself from its measured fread speed
// (auto). METABULI_PREFETCH=1 forces overlap on, =0 forces it off.
inline int ioPrefetchMode() { // 0 = forced off, 1 = forced on, 2 = auto
    static const int mode = []() {
        const char * v = getenv("METABULI_PREFETCH");
        if (v == nullptr) {
            return 2;
        }
        return (strcmp(v, "0") == 0) ? 0 : 1;
    }();
    return mode;
}

inline size_t profFread(void * buf, size_t size, size_t count, FILE * fp) {
    if (likely(!ioProfEnabled())) {
        return fread(buf, size, count, fp);
    }
    const auto t0 = std::chrono::steady_clock::now();
    const size_t readCount = fread(buf, size, count, fp);
    ioProfNanos().fetch_add(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t0).count(),
        std::memory_order_relaxed);
    ioProfBytes().fetch_add(readCount * size, std::memory_order_relaxed);
    return readCount;
}

template<typename T>
struct ReadBuffer {
    FILE * fp;
    T * p;
    size_t size;
    size_t capacity;
    T * start;
    T * end;
    // Double-buffered read-ahead (constructor-eligible + METABULI_PREFETCH=1):
    // two halves alternate between "being drained" and "being filled by a
    // background fread". Each half keeps a small headroom in front so
    // loadBuffer(unused) can place the carried-over tail right before the
    // prefetched data without moving it. A PERSISTENT worker thread serves the
    // reads - spawning a task per chunk costs millions of thread creations per
    // batch and dominates once reads hit the page cache.
    static constexpr size_t PREFETCH_HEADROOM = 16;
    T * allocBase = nullptr;
    T * half[2] = {nullptr, nullptr};
    int active = 0;
    bool prefetching = false;
    enum class PfState : uint8_t { IDLE, REQUESTED, READY, QUIT };
    std::thread pfWorker;
    std::mutex pfMutex;
    std::condition_variable pfCv;
    PfState pfState = PfState::IDLE;
    T * pfDst = nullptr;
    size_t pfReadCount = 0;
    double pfReadNs = 0.0;
    // Fetch ahead only after two consecutive full sequential reads, so
    // seek-heavy access does not keep reading chunks it then discards.
    uint8_t pfSeqRuns = 0;
    // Auto mode: overlap only while reads are device-bound. Every fread (sync
    // or worker) is timed; the EMA is updated by the consumer thread alone.
    // Thresholds are PER-THREAD speeds under full concurrency, where all
    // search threads share the device/DRAM bandwidth: measured ~34 MB/s per
    // thread on cold NVMe streams vs ~3.4 GB/s per thread on page-cache hits
    // (128 threads). Hysteresis sits in that two-orders-of-magnitude gap:
    // slower than 250 MB/s turns overlap on, faster than 1 GB/s turns it off.
    static constexpr double PF_ON_NS_PER_BYTE = 4.0;  // < 250 MB/s
    static constexpr double PF_OFF_NS_PER_BYTE = 1.0; // > 1 GB/s
    uint8_t pfMode = 0;          // resolved METABULI_PREFETCH (0/1/2=auto)
    bool pipelining = false;     // auto-mode decision; starts synchronous
    double emaNsPerByte = 0.0;   // 0 = no sample yet

    void updateEma(double ns, size_t bytes) {
        if (bytes == 0) {
            return;
        }
        const double nsPerByte = ns / double(bytes);
        emaNsPerByte = (emaNsPerByte == 0.0)
                           ? nsPerByte
                           : 0.875 * emaNsPerByte + 0.125 * nsPerByte;
    }

    // Timed synchronous read on the consumer thread (feeds the EMA).
    size_t timedFread(T * dst, size_t count) {
        const auto t0 = std::chrono::steady_clock::now();
        const size_t n = profFread(dst, sizeof(T), count, fp);
        const double ns = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - t0).count();
        updateEma(ns, n * sizeof(T));
        return n;
    }

    bool wantPipeline() {
        if (pfMode != 2) {
            return pfMode == 1;
        }
        if (emaNsPerByte > PF_ON_NS_PER_BYTE) {
            pipelining = true;
        } else if (emaNsPerByte != 0.0 && emaNsPerByte < PF_OFF_NS_PER_BYTE) {
            pipelining = false;
        }
        return pipelining;
    }

    ReadBuffer() : fp(nullptr), p(nullptr), size(0), capacity(0), start(nullptr), end(nullptr) {}

    explicit ReadBuffer(std::string file, size_t sizeOfBuffer=100, bool prefetchEligible=false) {
        fp = fopen(file.c_str(), "rb");
        if (!fp) {
            std::cerr << "Error opening file: " << file << std::endl;
            exit(EXIT_FAILURE);
        }
        capacity = sizeOfBuffer;
        size = 0;
        pfMode = (uint8_t) ioPrefetchMode();
        prefetching = prefetchEligible && pfMode != 0;
        if (prefetching) {
            const size_t halfStride = PREFETCH_HEADROOM + capacity;
            // Small tail pad keeps wide (8-byte) decoder loads near EOF inside
            // the allocation when the last half ends flush with it. The worker
            // thread is spawned lazily on the first kick, so a stream that
            // stays synchronous costs no extra thread.
            allocBase = (T *) calloc(2 * halfStride + PREFETCH_HEADROOM, sizeof(T));
            half[0] = allocBase;
            half[1] = allocBase + halfStride;
            start = half[0] + PREFETCH_HEADROOM;
        } else {
            allocBase = (T *) calloc(sizeOfBuffer, sizeof(T));
            start = allocBase;
        }
        p = start;
        end = start;
    };

    ~ReadBuffer() {
        if (pfWorker.joinable()) {
            {
                std::unique_lock<std::mutex> lock(pfMutex);
                // Let an in-flight read finish; never free/close under it.
                pfCv.wait(lock, [this] { return pfState != PfState::REQUESTED; });
                pfState = PfState::QUIT;
            }
            pfCv.notify_all();
            pfWorker.join();
        }
        if (fp) {
            fclose(fp);
        }
        if (allocBase) {
            free(allocBase);
        }
    };

    void pfLoop() {
        std::unique_lock<std::mutex> lock(pfMutex);
        for (;;) {
            pfCv.wait(lock, [this] {
                return pfState == PfState::REQUESTED || pfState == PfState::QUIT;
            });
            if (pfState == PfState::QUIT) {
                return;
            }
            T * dst = pfDst;
            lock.unlock();
            const auto t0 = std::chrono::steady_clock::now();
            const size_t n = profFread(dst, sizeof(T), capacity, fp);
            const double ns = (double) std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - t0).count();
            lock.lock();
            pfReadCount = n;
            pfReadNs = ns;
            pfState = PfState::READY;
            pfCv.notify_all();
        }
    }

    void kickPrefetch() {
        if (!pfWorker.joinable()) {
            pfWorker = std::thread(&ReadBuffer::pfLoop, this);
        }
        {
            std::lock_guard<std::mutex> lock(pfMutex);
            pfDst = half[1 - active] + PREFETCH_HEADROOM;
            pfState = PfState::REQUESTED;
        }
        pfCv.notify_all();
        ioPrefetchChunks().fetch_add(1, std::memory_order_relaxed);
    }

    // If a prefetch is in flight or done, wait for it and claim its result.
    bool reapPrefetch(size_t & readCount) {
        if (!pfWorker.joinable()) { // worker never started: nothing in flight
            return false;
        }
        std::unique_lock<std::mutex> lock(pfMutex);
        if (pfState == PfState::IDLE) {
            return false;
        }
        pfCv.wait(lock, [this] { return pfState == PfState::READY; });
        readCount = pfReadCount;
        pfState = PfState::IDLE;
        lock.unlock();
        updateEma(pfReadNs, readCount * sizeof(T)); // consumer-thread-only EMA
        return true;
    }

    size_t loadBuffer(size_t unused = 0) {
        if (prefetching) {
            size_t readCount;
            T * dst = half[1 - active] + PREFETCH_HEADROOM;
            if (reapPrefetch(readCount)) {
                // Fresh data sits in the other half; move the carried tail
                // over and swap. (The carry fits in the headroom.)
                if (unused > 0) {
                    std::memcpy(dst - unused, end - unused, unused * sizeof(T));
                }
                active = 1 - active;
            } else {
                // Synchronous read (cold start, sync mode, after a seek, or
                // past EOF): reuse the CURRENT half so a stream that stays
                // synchronous touches the same memory as non-prefetch mode.
                // Save the carry into the headroom before overwriting it
                // (end == dst after an EOF tail: already in place, skip).
                dst = half[active] + PREFETCH_HEADROOM;
                if (unused > 0 && end != dst) {
                    std::memcpy(dst - unused, end - unused, unused * sizeof(T));
                }
                readCount = timedFread(dst, capacity);
            }
            start = dst - unused;
            p = start;
            end = dst + readCount;
            size = readCount + unused;
            if (readCount == capacity) { // short read = EOF: nothing left to fetch
                if (pfSeqRuns < 2) {
                    ++pfSeqRuns;
                }
                if (pfSeqRuns >= 2 && wantPipeline()) {
                    kickPrefetch();
                }
            } else {
                pfSeqRuns = 0;
            }
            return readCount + unused;
        }
        memmove(start,                   // dest
                start + (size - unused), // src
                unused * sizeof(T));     // bytes
        size_t readCount = profFread(start + unused, sizeof(T), capacity - unused, fp);
        size = readCount + unused;
        end = start + readCount + unused;
        p = start;
        return readCount + unused;
    }

    inline T getNext() {
        if (p >= end) {
            if(loadBuffer() == 0) {
               return T(); // Return default value if no more data
            }
        }
        return *p++;
    }

    size_t loadBufferAt(size_t offset) {
        return loadBufferAt(offset, capacity);
    }

    // Seek and read at most maxCount elements (capped by capacity). Lets a
    // caller doing sparse random reads avoid paying a full-capacity read per
    // seek.
    size_t loadBufferAt(size_t offset, size_t maxCount) {
        if (prefetching) {
            size_t discarded;
            reapPrefetch(discarded); // repositioning: discard any fetched chunk
            pfSeqRuns = 0;
            fseek(fp, offset * sizeof(T), SEEK_SET);
            T * dst = half[active] + PREFETCH_HEADROOM;
            size_t readCount = timedFread(dst, std::min(maxCount, capacity));
            size = readCount;
            start = dst;
            end = dst + readCount;
            p = start;
            return readCount;
        }
        fseek(fp, offset * sizeof(T), SEEK_SET);
        size_t readCount = profFread(start, sizeof(T), std::min(maxCount, capacity), fp);
        size = readCount;
        end = start + readCount;
        p = start;
        return readCount;
    }
};

template<typename T>
struct WriteBuffer {
    FILE * fp;
    T * buffer;
    size_t capacity;
    size_t size;
    size_t writeCnt;

    explicit WriteBuffer(std::string file, size_t sizeOfBuffer=100) {
        if (file.empty()) {
            fp = nullptr;
            buffer = nullptr;
            capacity = 0;
            size = 0;
            writeCnt = 0;
            return;
        }
        fp = fopen(file.c_str(), "wb");
        if (!fp) {
            std::cerr << "Error opening file: " << file << std::endl;
            exit(EXIT_FAILURE);
        }
        buffer = (T *) calloc(sizeOfBuffer, sizeof(T));
        capacity = sizeOfBuffer;
        size = 0;
        writeCnt = 0;
    };

    ~WriteBuffer() {
        close();
        if (buffer) {
            free(buffer);
        }
    };

    void flush() {
        fwrite(buffer, sizeof(T), size, fp);
        size = 0;
    }

    void close() {
        if (fp) {
            if (size > 0) {
                flush();
            }
            fclose(fp);
            fp = nullptr;
        }
    }

    void write(T *data, size_t dataNum = 1) {
        if (size + dataNum > capacity) {
            flush();
        }
        memcpy(buffer + size, data, dataNum * sizeof(T));
        size += dataNum;
        writeCnt += dataNum;
    }
};

inline bool fileExist(const std::string& name) {
    if (FILE *file = fopen(name.c_str(), "r")) {
        fclose(file);
        return true;
    } else {
        return false;
    }
}

void process_mem_usage(double& vm_usage, double& resident_set);

TaxonomyWrapper * loadTaxonomy(const std::string & dbDir, const std::string & taxonomyDir = "");

TaxonomyWrapper * loadTaxonomyDB(const std::string & taxDbFile);

int loadDbParameters(LocalParameters & par, const std::string & dbDir);

// --- Sub-DB layout -----------------------------------------------------------
// A database directory holds its k-mer lists (diffIdx/info/split/taxID_list)
// in <DB dir>/kmer-lists/<name>/ — build creates kmer-lists/base, add-sub-db
// adds siblings. Databases from before this layout keep their single list at
// the top level; it resolves as sub-DB "base". Either way, all lists share the
// database-level taxonomyDB and therefore one internal tax-ID numbering.

// db.parameters keys guarding the shared internal tax-ID numbering.
#define SUBDB_KEY_LINEAGE "Taxonomy_lineage"
#define SUBDB_KEY_MAX_TAXID "Taxonomy_max_id"
// Directory that holds the k-mer lists.
#define SUBDB_DIR_NAME "kmer-lists"
// Default sub-DB name; also the selector of a legacy top-level k-mer list.
#define SUBDB_BASE_NAME "base"

struct SubDbEntry {
    std::string name; // user-facing selector (SUBDB_BASE_NAME for the top-level list)
    std::string dir;  // directory containing diffIdx/info/split/taxID_list
};

// Enumerate the k-mer lists inside dbDir and apply the user's selection.
// subDbList is a comma-separated list of sub-DB names ("" selects all).
// Exits with an error message if no list is found or a name is unknown.
std::vector<SubDbEntry> resolveSubDbs(const std::string & dbDir, const std::string & subDbList);

// Read one key from <dir>/db.parameters; returns "" if the file or key is absent.
std::string readDbParameter(const std::string & dir, const std::string & key);

// Append "key\tvalue" to <dir>/db.parameters (creates the file if needed).
void appendDbParameter(const std::string & dir, const std::string & key, const std::string & value);

// Verify that every selected sub-DB was built against the currently loaded
// taxonomy: same lineage ID and a build-time max tax ID not exceeding the
// loaded taxonomy's. Also checks that k-mer-space parameters agree with the
// umbrella. Exits on violation. No-op for a classic single-list database.
void validateSubDbs(const std::string & dbDir, const std::vector<SubDbEntry> & subDbs, int maxTaxID);

// Sum of Total_seq_length over the selected sub-DBs; falls back to
// readDbSize() for lists without the key. Species shared between sub-DBs are
// counted once per sub-DB (slight overestimation).
size_t totalDbSize(const std::vector<SubDbEntry> & subDbs);

// Random lowercase-hex string (used to mint taxonomy lineage IDs).
std::string randomHexString(size_t length);

// Number of logical CPUs (SMT threads included) the process may run on.
// Honors the cgroup/cpuset affinity (SLURM `-c N`, taskset, containers).
int availableCpus();

// Seed par.threads with availableCpus(). Call BEFORE parseParameters so the
// value appears in the parameter dump and an explicit --threads overrides it.
void setCpuThreadDefault(LocalParameters & par);

int searchAccession2TaxID(const std::string & name, const std::unordered_map<std::string, int> & acc2taxid);

template <typename T>
size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t number, int cnt) {
    bufferIdx = 0;
    fseek(fp, cnt * sizeof(T), SEEK_CUR);
    return fread(buffer, sizeof(T), number, fp);
}

template <typename T>
size_t loadBuffer(FILE *fp, T *buffer, size_t &bufferIdx, size_t number) {
    bufferIdx = 0;
    return fread(buffer, sizeof(T), number, fp);
}


// template <typename T>
// size_t loadBuffer2(FILE *fp, T *buffer, size_t number, int cnt) {
//     fseek(fp, cnt * sizeof(T), SEEK_CUR);
//     return fread(buffer, sizeof(T), number, fp);
// }

template <typename T>
size_t loadBuffer2(FILE *fp, T *buffer, size_t number, size_t unused) {
    memmove(buffer,
            buffer + (number - unused),
            unused * sizeof(T));
    size_t readCount = fread(buffer + unused, sizeof(T), number - unused, fp);
    return readCount + unused;
}

template <typename T>
size_t loadBuffer2(FILE *fp, T *buffer, size_t number) {
    return fread(buffer, sizeof(T), number, fp);
}


template <typename T>
inline T getElement(
    size_t bufferSize,
    FILE *kmerInfoFp,
    T *infoBuffer,
    size_t &infoBufferIdx) 
{
    if (unlikely(infoBufferIdx >= bufferSize)) {
        loadBuffer(kmerInfoFp, infoBuffer, infoBufferIdx, bufferSize,
                static_cast<int>(infoBufferIdx - bufferSize));
    }
    return infoBuffer[infoBufferIdx];
}

void getObservedAccessionList(const std::string & fnaListFileName,
                              std::vector<std::string> & fastaList,
                              std::unordered_map<std::string, TaxID> & acc2taxid);

void fillAcc2TaxIdMap(std::unordered_map<std::string, TaxID> & acc2taxid,
                      const std::string & acc2taxidFileName);
                                


struct SeqEntry {
    std::string name;
    std::string s;
    SeqEntry () = default;
    SeqEntry (const std::string & name, const std::string & sequence) 
        : name(name), s(sequence) {}
};
struct WorkItem {
    int bufferIndex;
    size_t sequenceCount;
    size_t writePos;
    uint32_t sequenceIdOffset;
};

// A simple thread-safe queue for passing chunks of reads
template<typename T>
class BlockingQueue {
private:
    std::mutex mtx;
    std::condition_variable cv;
    std::queue<T> q;

public:
    void push(T value) {
        {
            std::lock_guard<std::mutex> lock(mtx);
            q.push(std::move(value));
        }
        cv.notify_one();
    }

    // Tries to pop a value. Returns true on success, false if the queue is empty.
    bool try_pop(T& value) {
        std::lock_guard<std::mutex> lock(mtx);
        if (q.empty()) {
            return false;
        }
        value = std::move(q.front());
        q.pop();
        return true;
    }

    // Waits until a value is available and then pops it.
    void wait_and_pop(T& value) {
        std::unique_lock<std::mutex> lock(mtx);
        cv.wait(lock, [this]{ return !q.empty(); });
        value = std::move(q.front());
        q.pop();
    }
};

inline uint64_t extract_bits(uint64_t value, unsigned int start_bit, unsigned int block_size) {
    // Edge case handling: If block_size is 64, the mask is ~0ULL (all bits set).
    // (1ULL << 64) is Undefined Behavior, so we handle 64 separately.
    // if (block_size >= 64) {
    //     return value;
    // }

    // // Generalized Mask Creation: (2^L - 1) << start_bit
    // // 1. (1ULL << block_size) - 1  -> Creates the L-bit mask (0b11...1)
    // // 2. (...) << start_bit      -> Shifts the mask to the correct position
    // uint64_t full_mask = ((1ULL << block_size) - 1) << start_bit;

    // // Apply the mask and shift down (normalization)
    // return (value & full_mask) >> start_bit;
    return (value & (((1ULL << block_size) - 1) << start_bit)) >> start_bit;
}

inline unsigned bitsNeeded(unsigned N) {
    return (N <= 1) ? 0 : static_cast<unsigned>(std::ceil(std::log2(N)));
}

int hammingDist(const std::string & codon1, const std::string & codon2);


size_t readDbSize(const std::string& dbDir);

inline double computeLEM_logE(
    int queryLengthNt,
    int spanLengthNt,
    size_t dbSizeNt,
    double logP)
{
    const double aaQueryLen = queryLengthNt / 3.0;    
    const double aaMatchLen = spanLengthNt / 3.0;
    double lenDiff = std::max(0.0, aaQueryLen - aaMatchLen);
    return std::log((lenDiff + 1) * 6.0) + 
        std::log(dbSizeNt / 3.0) +
        logP;
}

uint32_t parseMask(const char* s);

uint32_t safe_right_shift_32(uint32_t value, unsigned int shift);
uint32_t safe_left_shift_32(uint32_t value, unsigned int shift);

std::string reverseComplement(std::string &read);

char *reverseComplement(char *read, size_t length);
int getFirstOneAfterFirstZero(uint32_t mask);

uint64_t disperseBits(uint64_t source, uint64_t pattern, int chunk_size);
uint64_t stretchBits(uint64_t pattern, int repeat_count);

#endif //ADCLASSIFIER2_COMMON_H
