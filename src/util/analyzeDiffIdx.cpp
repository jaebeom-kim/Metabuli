#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "Kmer.h"
#include "Mmap.h"
#include "DeltaIdxReader.h"
#include <cstdint>
#include <cstdio>
#include <vector>
#include <zlib.h>
#ifdef OPENMP
#include <omp.h>
#endif

using namespace std;

// Hidden developer tool: streams the 15-bit varint diffIdx, reconstructs each
// delta, and reports the delta distribution, identical-k-mer (delta = 0) group
// sizes, and projected sizes under alternative encodings (LEB128, flagless)
// plus a generic-compression estimate. Used to decide how to shrink diffIdx
// for cache-constrained machines.
//
// The scan is parallelized over the sibling split file's anchors: every anchor
// byte offset is delta-aligned, and anchors sit at AA-group starts, so each
// range decodes complete deltas independently. Identical-k-mer runs CAN span a
// range start (the anchor k-mer's duplicates follow it), so per-range runs are
// stitched sequentially at merge time.
//
// Usage: metabuli analyze-diffidx <diffIdx> [maxBytes]

namespace {

// Exact sizes 1..8, then power-of-two buckets (9-16, 17-32, ...).
constexpr int GROUP_BUCKETS = 40;

int groupBucket(uint64_t s) {
    if (s <= 8) {
        return (int)(s - 1);
    }
    const int lg = 63 - __builtin_clzll(s - 1); // s=9..16 -> 3, 17..32 -> 4, ...
    const int idx = 8 + (lg - 3);
    return idx < GROUP_BUCKETS ? idx : GROUP_BUCKETS - 1;
}

void groupBucketLabel(int idx, char *out, size_t n) {
    if (idx < 8) {
        snprintf(out, n, "%d", idx + 1);
    } else {
        const uint64_t lo = (1ULL << (idx - 8 + 3)) + 1;
        const uint64_t hi = 1ULL << (idx - 8 + 4);
        snprintf(out, n, "%llu-%llu", (unsigned long long)lo, (unsigned long long)hi);
    }
}

struct ScanStats {
    uint64_t nKmer = 0;
    uint64_t nFragments = 0;
    uint64_t nZeroDelta = 0;
    uint64_t fragHist[6] = {0};
    uint64_t bitHist[65] = {0};
    uint64_t curBytes = 0;
    uint64_t leb128Bytes = 0;
    uint64_t flaglessBytes = 0;
    // Identical-value run bookkeeping. A group = one k-mer value's entries:
    // a nonzero delta starts a group, following zero deltas extend it.
    uint64_t extendPrev = 0;   // zero deltas before this range's first nonzero
    bool sawNonzero = false;
    uint64_t openGroup = 0;    // size of the group still open at range end
    uint64_t groupCnt = 0;     // groups CLOSED inside this range
    uint64_t groupHist[GROUP_BUCKETS] = {0};       // closed groups by size
    uint64_t groupKmerHist[GROUP_BUCKETS] = {0};   // k-mers in those groups
    uint64_t maxGroup = 0;

    void closeGroup(uint64_t s) {
        ++groupCnt;
        const int b = groupBucket(s);
        groupHist[b]++;
        groupKmerHist[b] += s;
        if (s > maxGroup) {
            maxGroup = s;
        }
    }
};

// Decode complete deltas in [byteStart, byteEnd). byteStart must be
// delta-aligned; a partial delta at byteEnd (capped sample) is dropped.
bool scanRange(const string &fileName, uint64_t byteStart, uint64_t byteEnd,
               ScanStats &st) {
    FILE *fp = fopen(fileName.c_str(), "rb");
    if (!fp) {
        return false;
    }
    fseek(fp, (long)byteStart, SEEK_SET);
    const size_t BUF = 1u << 22; // 4M uint16 = 8 MB
    vector<uint16_t> buf(BUF);
    uint64_t remaining = (byteEnd - byteStart) / 2; // uint16 fragments
    uint64_t delta = 0;
    uint32_t fragsThisDelta = 0;
    while (remaining > 0) {
        const size_t want = remaining < BUF ? (size_t)remaining : BUF;
        const size_t n = fread(buf.data(), sizeof(uint16_t), want, fp);
        if (n == 0) {
            break;
        }
        remaining -= n;
        for (size_t i = 0; i < n; ++i) {
            const uint16_t w = buf[i];
            ++st.nFragments;
            ++fragsThisDelta;
            delta = (delta << 15) | (uint16_t)(w & 0x7FFF);
            if (w & 0x8000) {
                ++st.nKmer;
                if (delta == 0) {
                    ++st.nZeroDelta;
                    if (!st.sawNonzero) {
                        ++st.extendPrev;
                    } else {
                        ++st.openGroup;
                    }
                } else {
                    if (st.sawNonzero) {
                        st.closeGroup(st.openGroup);
                    }
                    st.sawNonzero = true;
                    st.openGroup = 1;
                }
                if (fragsThisDelta <= 5) {
                    ++st.fragHist[fragsThisDelta];
                }
                int bits = 0;
                uint64_t d = delta;
                while (d) { ++bits; d >>= 1; }
                st.bitHist[bits]++;
                st.curBytes += (uint64_t)fragsThisDelta * 2;
                int leb = (bits + 6) / 7;
                if (leb < 1) leb = 1;
                st.leb128Bytes += leb;
                int groups = (bits + 14) / 15;
                if (groups < 1) groups = 1;
                st.flaglessBytes += ((groups * 15 + 2 + 7) / 8); // 15-bit groups + 2-bit length tag
                delta = 0;
                fragsThisDelta = 0;
            }
        }
    }
    fclose(fp);
    return true;
}

// Stitch per-range stats in file order: a range's leading zero run extends the
// group left open by the previous range.
void mergeInOrder(const vector<ScanStats> &parts, ScanStats &total) {
    uint64_t carry = 0;
    for (const ScanStats &p : parts) {
        total.nKmer += p.nKmer;
        total.nFragments += p.nFragments;
        total.nZeroDelta += p.nZeroDelta;
        for (int i = 0; i < 6; ++i) total.fragHist[i] += p.fragHist[i];
        for (int i = 0; i < 65; ++i) total.bitHist[i] += p.bitHist[i];
        total.curBytes += p.curBytes;
        total.leb128Bytes += p.leb128Bytes;
        total.flaglessBytes += p.flaglessBytes;
        for (int i = 0; i < GROUP_BUCKETS; ++i) {
            total.groupHist[i] += p.groupHist[i];
            total.groupKmerHist[i] += p.groupKmerHist[i];
        }
        if (p.maxGroup > total.maxGroup) total.maxGroup = p.maxGroup;
        carry += p.extendPrev;
        if (p.sawNonzero) {
            if (carry > 0) {
                total.closeGroup(carry); // group left open by the previous range
            }
            carry = p.openGroup;
        } // else: the whole range was zero deltas; keep extending
    }
    if (carry > 0) {
        total.closeGroup(carry); // the file's final group
    }
    // Total group count = everything in the histogram (in-range closures came
    // in via groupHist above, stitched ones via closeGroup).
    total.groupCnt = 0;
    for (int i = 0; i < GROUP_BUCKETS; ++i) total.groupCnt += total.groupHist[i];
}

} // namespace

int analyzeDiffIdx(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string diffIdxFileName = par.filenames[0];
    // Optional second positional: cap bytes scanned for a quick serial sample.
    uint64_t maxBytes = (par.filenames.size() >= 2)
                            ? strtoull(par.filenames[1].c_str(), nullptr, 10) : 0;

    if (DeltaIdxReader::loadDiffidxFormat(diffIdxFileName) == 1) {
        fprintf(stderr, "Error: this diffIdx is LEB128 (Diffidx_format 1); "
                        "analyze-diffidx reads the 15-bit varint format only.\n");
        return 1;
    }
    const uint64_t fileSize = FileUtil::getFileSize(diffIdxFileName);
    if (fileSize == 0) {
        fprintf(stderr, "Error: %s is empty or missing.\n", diffIdxFileName.c_str());
        return 1;
    }

    // Range boundaries from the sibling split file (delta-aligned anchors).
    vector<uint64_t> bounds;
    bounds.push_back(0);
    const bool capped = maxBytes > 0 && maxBytes < fileSize;
    if (!capped) {
        const size_t slash = diffIdxFileName.find_last_of('/');
        const string splitFileName = (slash == string::npos)
            ? "split" : diffIdxFileName.substr(0, slash + 1) + "split";
        if (FileUtil::fileExists(splitFileName.c_str())) {
            MmapedData<DiffIdxSplit> splits = mmapData<DiffIdxSplit>(splitFileName.c_str(), 3);
            const size_t numSplits = splits.fileSize / sizeof(DiffIdxSplit);
            for (size_t i = 1; i < numSplits; ++i) {
                if (splits.data[i].ADkmer == 0 || splits.data[i].ADkmer == UINT64_MAX) {
                    continue;
                }
                const uint64_t off = splits.data[i].diffIdxOffset * 2; // uint16 -> bytes
                if (off > bounds.back() && off < fileSize) {
                    bounds.push_back(off);
                }
            }
            munmap(splits.data, splits.fileSize);
        }
    }
    bounds.push_back(capped ? (maxBytes & ~1ULL) : fileSize);
    const size_t nRanges = bounds.size() - 1;

    vector<ScanStats> parts(nRanges);
    bool ioError = false;
    #pragma omp parallel for schedule(dynamic)
    for (size_t r = 0; r < nRanges; ++r) {
        if (!scanRange(diffIdxFileName, bounds[r], bounds[r + 1], parts[r])) {
            #pragma omp critical
            ioError = true;
        }
    }
    if (ioError) {
        perror(diffIdxFileName.c_str());
        return 1;
    }
    ScanStats st;
    mergeInOrder(parts, st);
    const uint64_t bytesScanned = bounds.back();

    printf("=== diffIdx delta analysis%s (%zu ranges) ===\n",
           capped ? " (CAPPED sample)" : "", nRanges);
    printf("k-mers (deltas)      : %llu\n", (unsigned long long)st.nKmer);
    printf("identical k-mers (0) : %llu  (%.2f%%)\n", (unsigned long long)st.nZeroDelta,
           st.nKmer ? 100.0*st.nZeroDelta/st.nKmer : 0);
    printf("distinct k-mers      : %llu  (avg copies per value: %.3f)\n",
           (unsigned long long)st.groupCnt,
           st.groupCnt ? (double)st.nKmer/st.groupCnt : 0);
    printf("total fragments      : %llu\n", (unsigned long long)st.nFragments);
    printf("avg fragments/kmer   : %.3f\n", st.nKmer ? (double)st.nFragments/st.nKmer : 0);
    printf("bytes scanned        : %.2f GB\n", bytesScanned/1e9);

    printf("\nidentical-value group-size histogram (groups / k-mers covered):\n");
    for (int i = 0; i < GROUP_BUCKETS; ++i) {
        if (st.groupHist[i] == 0) {
            continue;
        }
        char label[32];
        groupBucketLabel(i, label, sizeof(label));
        printf("  size %-9s: %12llu groups (%.2f%%)  %12llu k-mers (%.2f%%)\n",
               label,
               (unsigned long long)st.groupHist[i],
               st.groupCnt ? 100.0*st.groupHist[i]/st.groupCnt : 0,
               (unsigned long long)st.groupKmerHist[i],
               st.nKmer ? 100.0*st.groupKmerHist[i]/st.nKmer : 0);
    }
    printf("  max group size: %llu\n", (unsigned long long)st.maxGroup);

    printf("\nfragment-count histogram (deltas needing N fragments):\n");
    for (int f = 1; f <= 5; ++f)
        if (st.fragHist[f]) printf("  %d frag: %12llu  (%.2f%%)\n", f,
            (unsigned long long)st.fragHist[f], 100.0*st.fragHist[f]/st.nKmer);
    printf("\ndelta bit-width histogram (grouped):\n");
    struct { int lo, hi; } groups[] = {{0,7},{8,15},{16,23},{24,31},{32,47},{48,64}};
    for (auto g : groups) {
        uint64_t s = 0; for (int b = g.lo; b <= g.hi; ++b) s += st.bitHist[b];
        if (s) printf("  %2d-%2d bits: %12llu  (%.2f%%)\n", g.lo, g.hi,
            (unsigned long long)s, 100.0*s/st.nKmer);
    }
    printf("\n=== projected diffIdx size ===\n");
    printf("current (15-bit varint) : %.2f GB  (%.3f B/kmer)\n",
           st.curBytes/1e9, st.nKmer ? (double)st.curBytes/st.nKmer : 0);
    if (st.curBytes) {
        printf("LEB128 (byte varint)    : %.2f GB  (%.1f%% of current, %.3f B/kmer)\n",
               st.leb128Bytes/1e9, 100.0*st.leb128Bytes/st.curBytes,
               st.nKmer ? (double)st.leb128Bytes/st.nKmer : 0);
        printf("flagless 15b + len tag  : %.2f GB  (%.1f%% of current)\n",
               st.flaglessBytes/1e9, 100.0*st.flaglessBytes/st.curBytes);
    }

    // Generic-compression estimate on the first 256 MB.
    const size_t SAMPLE = 256u << 20;
    const size_t sampleSize = (size_t)(fileSize < SAMPLE ? fileSize : SAMPLE);
    FILE *fp = fopen(diffIdxFileName.c_str(), "rb");
    if (fp) {
        vector<unsigned char> sample(sampleSize);
        const size_t got = fread(sample.data(), 1, sampleSize, fp);
        fclose(fp);
        if (got > 0) {
            uLongf bound = compressBound(got);
            vector<unsigned char> out(bound);
            if (compress2(out.data(), &bound, sample.data(), got, 9) == Z_OK) {
                printf("\nzlib -9 on %.0f MB raw sample: %.1f%% of original\n"
                       "  (generic recompression upper bound; zstd-19 usually beats this by 10-20%%)\n",
                       got/1e6, 100.0*bound/got);
            }
        }
    }
    return 0;
}
