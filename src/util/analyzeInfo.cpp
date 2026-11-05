#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "Kmer.h"
#include "Mmap.h"
#include "DeltaIdxReader.h"
#include "InfoIndex.h"
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <algorithm>
#include <vector>
#ifdef OPENMP
#include <omp.h>
#endif

using namespace std;

// Hidden developer tool: walks diffIdx and info in lockstep and sizes the
// info-index reduction levers discussed for the dedup format:
//   1. ID frequency distribution -> entropy bound and a frequency-ranked
//      class-coding projection (requires renumbering by k-mer count).
//   2. Within-group ID deltas (IDs sharing one k-mer value, sorted) -> a
//      group-delta varint projection under the CURRENT numbering (upper
//      bound; a taxonomy-aware renumbering would shrink it further).
//   3. HyperLogLog estimate of DISTINCT ID sets -> whether a set dictionary
//      is worth pursuing.
// Parallelized over split anchors like analyze-diffidx; identical-value runs
// crossing a range start are stitched at merge time.
//
// Usage: metabuli analyze-info <DB dir with diffIdx/info/split>

namespace {

uint64_t mix64(uint64_t x) {
    x ^= x >> 30; x *= 0xBF58476D1CE4E5B9ULL;
    x ^= x >> 27; x *= 0x94D049BB133111EBULL;
    x ^= x >> 31;
    return x;
}

// HyperLogLog, p = 14 (16K registers, ~0.8% standard error).
struct Hll {
    static constexpr int P = 14;
    vector<uint8_t> reg;
    Hll() : reg(1u << P, 0) {}
    void add(uint64_t h) {
        const uint32_t idx = (uint32_t)(h >> (64 - P));
        const uint64_t rest = h << P;
        const uint8_t rank = rest == 0 ? (uint8_t)(64 - P + 1)
                                       : (uint8_t)(__builtin_clzll(rest) + 1);
        if (rank > reg[idx]) {
            reg[idx] = rank;
        }
    }
    void merge(const Hll &o) {
        for (size_t i = 0; i < reg.size(); ++i) {
            if (o.reg[i] > reg[i]) reg[i] = o.reg[i];
        }
    }
    double estimate() const {
        const double m = (double)reg.size();
        double sum = 0;
        size_t zeros = 0;
        for (uint8_t r : reg) {
            sum += pow(2.0, -(double)r);
            zeros += (r == 0);
        }
        double e = 0.7213 / (1.0 + 1.079 / m) * m * m / sum;
        if (e <= 2.5 * m && zeros > 0) {
            e = m * log(m / (double)zeros); // small-range correction
        }
        return e;
    }
};

size_t lebBytes(uint64_t v) {
    int bits = 0;
    while (v) { ++bits; v >>= 1; }
    const int b = (bits + 6) / 7;
    return b < 1 ? 1 : (size_t)b;
}

// Aggregates filled per group of IDs sharing one k-mer value.
struct GroupAgg {
    uint64_t groups = 0;
    uint64_t multiGroups = 0;
    uint64_t multiIds = 0;
    uint64_t deltaBitsHist[65] = {0};
    uint64_t groupVarintBytes = 0; // LEB(firstId) + sum LEB(sorted deltas)
    Hll hllAll;
    Hll hllMulti;
    vector<uint32_t> scratch;

    void processGroup(const vector<uint32_t> &ids) {
        ++groups;
        scratch = ids;
        sort(scratch.begin(), scratch.end());
        uint64_t h = 1469598103934665603ULL;
        for (uint32_t id : scratch) {
            h = mix64(h ^ id);
        }
        hllAll.add(h);
        groupVarintBytes += lebBytes(scratch[0]);
        if (scratch.size() > 1) {
            ++multiGroups;
            multiIds += scratch.size();
            hllMulti.add(h);
            for (size_t i = 1; i < scratch.size(); ++i) {
                const uint64_t d = scratch[i] - scratch[i - 1];
                int bits = 0;
                uint64_t x = d;
                while (x) { ++bits; x >>= 1; }
                deltaBitsHist[bits]++;
                groupVarintBytes += lebBytes(d);
            }
        }
    }

    void merge(const GroupAgg &o) {
        groups += o.groups;
        multiGroups += o.multiGroups;
        multiIds += o.multiIds;
        for (int i = 0; i < 65; ++i) deltaBitsHist[i] += o.deltaBitsHist[i];
        groupVarintBytes += o.groupVarintBytes;
        hllAll.merge(o.hllAll);
        hllMulti.merge(o.hllMulti);
    }
};

// Streams the deltas of one diffIdx byte range (format 0 or 1).
class DeltaRangeWalker {
    FILE *fp = nullptr;
    vector<uint8_t> buf;
    uint8_t *p = nullptr;
    uint8_t *end = nullptr;
    uint64_t remaining; // range bytes not yet read into buf
    const int fmt;

    bool refill() {
        const size_t carry = (size_t)(end - p);
        memmove(buf.data(), p, carry);
        const size_t want = (size_t)min<uint64_t>(remaining, buf.size() - carry);
        const size_t got = fread(buf.data() + carry, 1, want, fp);
        remaining -= got;
        p = buf.data();
        end = buf.data() + carry + got;
        return carry + got > 0;
    }

public:
    DeltaRangeWalker(const string &fileName, uint64_t byteStart, uint64_t byteEnd, int fmt)
        : buf(8u << 20), remaining(byteEnd - byteStart), fmt(fmt) {
        fp = fopen(fileName.c_str(), "rb");
        if (fp) {
            fseek(fp, (long)byteStart, SEEK_SET);
        }
        p = end = buf.data();
    }
    ~DeltaRangeWalker() {
        if (fp) fclose(fp);
    }
    bool ok() const { return fp != nullptr; }

    bool next(uint64_t &delta) {
        if (end - p < 10 && remaining > 0) {
            refill();
        }
        if (p >= end) {
            return false;
        }
        delta = 0;
        if (fmt == 1) { // LEB128: 7-bit fragments, end flag on the last byte
            uint8_t b;
            do {
                b = *p++;
                delta = (delta << 7) | (b & 0x7F);
            } while ((b & 0x80) == 0 && p < end);
        } else { // 15-bit fragments in uint16 words
            uint16_t w;
            do {
                memcpy(&w, p, 2);
                p += 2;
                delta = (delta << 15) | (w & 0x7FFF);
            } while ((w & 0x8000) == 0 && p < end);
        }
        return true;
    }
};

// Per-range result for sequential stitching of boundary-crossing runs.
struct RangeOut {
    vector<uint32_t> extendIds; // IDs of the leading zero-delta run
    vector<uint32_t> openIds;   // IDs of the group still open at range end
    bool sawNonzero = false;
    uint64_t idCount = 0;
};

} // namespace

int analyzeInfo(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string dbDir = par.filenames[0];
    const string diffFile = dbDir + "/diffIdx";
    const string infoFile = dbDir + "/info";
    const string splitFile = dbDir + "/split";

    const InfoIndexMetadata meta = InfoIndex::loadMetadata(infoFile);
    if (!meta.isPacked()) {
        fprintf(stderr, "Error: analyze-info currently supports packed info only "
                        "(this DB: %s).\n", meta.format.c_str());
        return 1;
    }
    const int fmt = DeltaIdxReader::loadDiffidxFormat(diffFile);
    const uint64_t diffSize = FileUtil::getFileSize(diffFile);
    const size_t idCount = meta.idCount;
    const size_t idSpace = (size_t)1 << meta.idBits;
    const double infoGB = InfoIndex::packedWordCount(idCount, meta.idBits) * 8 / 1e9;

    // Range boundaries: (diffIdx byte offset, info ordinal) pairs from split.
    vector<pair<uint64_t, uint64_t>> bounds;
    bounds.push_back({0, 0});
    if (FileUtil::fileExists(splitFile.c_str())) {
        MmapedData<DiffIdxSplit> splits = mmapData<DiffIdxSplit>(splitFile.c_str(), 3);
        const size_t numSplits = splits.fileSize / sizeof(DiffIdxSplit);
        for (size_t i = 1; i < numSplits; ++i) {
            if (splits.data[i].ADkmer == 0 || splits.data[i].ADkmer == UINT64_MAX) {
                continue;
            }
            const uint64_t off = (fmt == 1) ? splits.data[i].diffIdxOffset
                                            : splits.data[i].diffIdxOffset * 2;
            if (off > bounds.back().first && off < diffSize) {
                bounds.push_back({off, splits.data[i].infoIdxOffset});
            }
        }
        munmap(splits.data, splits.fileSize);
    }
    bounds.push_back({diffSize, idCount});
    const size_t nRanges = bounds.size() - 1;

    vector<RangeOut> outs(nRanges);
    vector<uint64_t> freq(idSpace, 0); // merged after the parallel region
    GroupAgg agg;                      // merged after the parallel region
    bool ioError = false;

    #pragma omp parallel
    {
        vector<uint64_t> freqLocal(idSpace, 0);
        GroupAgg aggLocal;
        vector<uint32_t> groupIds;
        groupIds.reserve(1024);

        #pragma omp for schedule(dynamic)
        for (size_t r = 0; r < nRanges; ++r) {
            DeltaRangeWalker walker(diffFile, bounds[r].first, bounds[r + 1].first, fmt);
            InfoIndexReader infoReader(infoFile, 1 << 20);
            if (!walker.ok() || !infoReader.loadAt(bounds[r].second)) {
                #pragma omp critical
                ioError = true;
                continue;
            }
            RangeOut &out = outs[r];
            const uint64_t nIds = bounds[r + 1].second - bounds[r].second;
            groupIds.clear();
            uint64_t delta;
            for (uint64_t k = 0; k < nIds; ++k) {
                uint32_t id;
                if (!walker.next(delta) || !infoReader.next(id)) {
                    #pragma omp critical
                    ioError = true;
                    break;
                }
                freqLocal[id]++;
                ++out.idCount;
                if (delta == 0) {
                    if (!out.sawNonzero) {
                        out.extendIds.push_back(id);
                    } else {
                        groupIds.push_back(id);
                    }
                } else {
                    if (out.sawNonzero) {
                        aggLocal.processGroup(groupIds);
                    }
                    out.sawNonzero = true;
                    groupIds.clear();
                    groupIds.push_back(id);
                }
            }
            if (out.sawNonzero) {
                out.openIds = groupIds;
            }
            groupIds.clear();
        }

        #pragma omp critical
        {
            for (size_t i = 0; i < idSpace; ++i) {
                freq[i] += freqLocal[i];
            }
            agg.merge(aggLocal);
        }
    }
    if (ioError) {
        fprintf(stderr, "Error: diffIdx/info walk failed or streams disagree.\n");
        return 1;
    }

    // Stitch identical-value runs that cross range starts.
    {
        vector<uint32_t> carry;
        for (const RangeOut &o : outs) {
            carry.insert(carry.end(), o.extendIds.begin(), o.extendIds.end());
            if (o.sawNonzero) {
                if (!carry.empty()) {
                    agg.processGroup(carry);
                }
                carry = o.openIds;
            }
        }
        if (!carry.empty()) {
            agg.processGroup(carry);
        }
    }

    uint64_t totalIds = 0;
    for (const RangeOut &o : outs) totalIds += o.idCount;

    // --- 1. frequency / entropy ---------------------------------------------
    uint64_t distinctIds = 0;
    double entropyBits = 0;
    for (size_t i = 0; i < idSpace; ++i) {
        if (freq[i] == 0) continue;
        ++distinctIds;
        const double pI = (double)freq[i] / (double)totalIds;
        entropyBits -= pI * log2(pI);
    }
    vector<uint64_t> sortedFreq;
    sortedFreq.reserve(distinctIds);
    for (size_t i = 0; i < idSpace; ++i) {
        if (freq[i]) sortedFreq.push_back(freq[i]);
    }
    sort(sortedFreq.rbegin(), sortedFreq.rend());
    // Frequency-ranked class coding: 2-bit selector + 7/14/21-bit payload.
    uint64_t classBits = 0, top128 = 0, top16k = 0;
    for (size_t rank = 0; rank < sortedFreq.size(); ++rank) {
        const uint64_t f = sortedFreq[rank];
        classBits += f * (rank < 128 ? 9 : rank < 16384 ? 16 : 23);
        if (rank < 128) top128 += f;
        if (rank < 16384) top16k += f;
    }

    const double curBitsPerId = 64.0 / InfoIndex::idsPerWord(meta.idBits);
    printf("=== info analysis (%zu ranges) ===\n", nRanges);
    printf("ids                  : %llu  (packed %u bits, %.1f effective, %.2f GB)\n",
           (unsigned long long)totalIds, meta.idBits, curBitsPerId, infoGB);
    printf("distinct ids         : %llu\n", (unsigned long long)distinctIds);
    printf("\n--- lever 1: frequency coding (needs renumbering by k-mer count) ---\n");
    printf("entropy              : %.2f bits/id  -> %.2f GB (theoretical bound)\n",
           entropyBits, entropyBits * totalIds / 8 / 1e9);
    printf("ranked classes 7/14/21b: %.2f bits/id -> %.2f GB (%.1f%% of current)\n",
           (double)classBits / totalIds, classBits / 8.0 / 1e9,
           100.0 * (classBits / 8.0) / (infoGB * 1e9));
    printf("top-128 ids cover    : %.2f%% of ids; top-16384: %.2f%%\n",
           100.0 * top128 / totalIds, 100.0 * top16k / totalIds);

    // --- 2. within-group deltas ---------------------------------------------
    printf("\n--- lever 2: within-group id deltas (needs dedup format) ---\n");
    printf("groups               : %llu  (multi: %llu, ids in multi: %llu = %.1f%%)\n",
           (unsigned long long)agg.groups, (unsigned long long)agg.multiGroups,
           (unsigned long long)agg.multiIds, 100.0 * agg.multiIds / totalIds);
    printf("sorted in-group delta bit-width histogram:\n");
    struct { int lo, hi; } gs[] = {{0,7},{8,10},{11,14},{15,17},{18,21},{22,32}};
    const uint64_t nDeltas = agg.multiIds - agg.multiGroups;
    for (auto g : gs) {
        uint64_t s = 0;
        for (int b = g.lo; b <= g.hi; ++b) s += agg.deltaBitsHist[b];
        if (s) printf("  %2d-%2d bits: %12llu  (%.2f%%)\n", g.lo, g.hi,
                      (unsigned long long)s, nDeltas ? 100.0 * s / nDeltas : 0);
    }
    printf("group varint projection (LEB first id + LEB sorted deltas, current numbering):\n");
    printf("  %.2f GB (%.1f%% of current info)\n",
           agg.groupVarintBytes / 1e9, 100.0 * agg.groupVarintBytes / (infoGB * 1e9));

    // --- 3. distinct id sets -------------------------------------------------
    const double setsAll = agg.hllAll.estimate();
    const double setsMulti = agg.hllMulti.estimate();
    printf("\n--- lever 3: id-set dictionary (HLL estimates, ~1%% error) ---\n");
    printf("distinct id sets     : ~%.3g of %llu groups (%.1f%% repeated)\n",
           setsAll, (unsigned long long)agg.groups,
           agg.groups ? 100.0 * (1.0 - setsAll / agg.groups) : 0);
    printf("distinct multi sets  : ~%.3g of %llu multi groups (%.1f%% repeated)\n",
           setsMulti, (unsigned long long)agg.multiGroups,
           agg.multiGroups ? 100.0 * (1.0 - setsMulti / agg.multiGroups) : 0);
    printf("set-reference cost   : %.1f bits/group -> %.2f GB refs (+ dictionary)\n",
           setsAll > 1 ? log2(setsAll) : 0,
           agg.groups * (setsAll > 1 ? log2(setsAll) : 0) / 8 / 1e9);

    return 0;
}
