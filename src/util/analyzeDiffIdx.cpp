#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include <cstdint>
#include <cstdio>
#include <vector>
#include <zlib.h>

using namespace std;

// Hidden developer tool: streams the 15-bit varint diffIdx, reconstructs each
// delta, and reports the delta distribution plus projected sizes under
// alternative encodings (LEB128, flagless) and a generic-compression estimate.
// Used to decide how to shrink diffIdx for cache-constrained machines.
//
// Usage: metabuli analyze-diffidx <diffIdx> [maxBytes]
int analyzeDiffIdx(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string diffIdxFileName = par.filenames[0];
    // Optional second positional: cap bytes scanned for a quick sample.
    uint64_t maxBytes = (par.filenames.size() >= 2)
                            ? strtoull(par.filenames[1].c_str(), nullptr, 10) : 0;

    FILE *fp = fopen(diffIdxFileName.c_str(), "rb");
    if (!fp) { perror(diffIdxFileName.c_str()); return 1; }

    const size_t BUF = 1u << 22; // 4M uint16 = 8 MB
    vector<uint16_t> buf(BUF);

    uint64_t nKmer = 0, nFragments = 0;
    uint64_t fragHist[6] = {0};
    uint64_t bitHist[65] = {0};
    uint64_t curBytes = 0, leb128Bytes = 0, flaglessBytes = 0;

    const size_t SAMPLE = 256u << 20; // 256 MB raw-byte sample for the zlib estimate
    vector<unsigned char> sample;
    sample.reserve(SAMPLE);

    uint64_t delta = 0;
    uint32_t fragsThisDelta = 0;
    uint64_t bytesRead = 0;
    bool capped = false;
    size_t n;
    while ((n = fread(buf.data(), sizeof(uint16_t), BUF, fp)) > 0) {
        for (size_t i = 0; i < n; ++i) {
            uint16_t w = buf[i];
            if (sample.size() < SAMPLE) {
                sample.push_back((unsigned char)(w & 0xFF));
                sample.push_back((unsigned char)(w >> 8));
            }
            ++nFragments;
            ++fragsThisDelta;
            delta = (delta << 15) | (uint16_t)(w & 0x7FFF);
            if (w & 0x8000) {
                ++nKmer;
                if (fragsThisDelta <= 5) ++fragHist[fragsThisDelta];
                int bits = 0; uint64_t d = delta; while (d) { ++bits; d >>= 1; }
                bitHist[bits]++;
                curBytes += (uint64_t)fragsThisDelta * 2;
                int leb = (bits + 6) / 7; if (leb < 1) leb = 1;
                leb128Bytes += leb;
                int groups = (bits + 14) / 15; if (groups < 1) groups = 1;
                flaglessBytes += ((groups * 15 + 2 + 7) / 8); // 15-bit groups + 2-bit length tag
                delta = 0; fragsThisDelta = 0;
            }
        }
        bytesRead += n * sizeof(uint16_t);
        if (maxBytes && bytesRead >= maxBytes) { capped = true; break; }
    }
    fclose(fp);

    printf("=== diffIdx delta analysis%s ===\n", capped ? " (CAPPED sample)" : "");
    printf("k-mers (deltas)      : %llu\n", (unsigned long long)nKmer);
    printf("total fragments      : %llu\n", (unsigned long long)nFragments);
    printf("avg fragments/kmer   : %.3f\n", nKmer ? (double)nFragments/nKmer : 0);
    printf("bytes scanned        : %.2f GB\n", bytesRead/1e9);
    printf("\nfragment-count histogram (deltas needing N fragments):\n");
    for (int f = 1; f <= 5; ++f)
        if (fragHist[f]) printf("  %d frag: %12llu  (%.2f%%)\n", f,
            (unsigned long long)fragHist[f], 100.0*fragHist[f]/nKmer);
    printf("\ndelta bit-width histogram (grouped):\n");
    struct { int lo, hi; } groups[] = {{0,7},{8,15},{16,23},{24,31},{32,47},{48,64}};
    for (auto g : groups) {
        uint64_t s = 0; for (int b = g.lo; b <= g.hi; ++b) s += bitHist[b];
        if (s) printf("  %2d-%2d bits: %12llu  (%.2f%%)\n", g.lo, g.hi,
            (unsigned long long)s, 100.0*s/nKmer);
    }
    printf("\n=== projected diffIdx size ===\n");
    printf("current (15-bit varint) : %.2f GB  (%.3f B/kmer)\n",
           curBytes/1e9, nKmer ? (double)curBytes/nKmer : 0);
    if (curBytes) {
        printf("LEB128 (byte varint)    : %.2f GB  (%.1f%% of current, %.3f B/kmer)\n",
               leb128Bytes/1e9, 100.0*leb128Bytes/curBytes, nKmer ? (double)leb128Bytes/nKmer : 0);
        printf("flagless 15b + len tag  : %.2f GB  (%.1f%% of current)\n",
               flaglessBytes/1e9, 100.0*flaglessBytes/curBytes);
    }
    if (!sample.empty()) {
        uLongf bound = compressBound(sample.size());
        vector<unsigned char> out(bound);
        if (compress2(out.data(), &bound, sample.data(), sample.size(), 9) == Z_OK) {
            printf("\nzlib -9 on %.0f MB raw sample: %.1f%% of original\n"
                   "  (generic recompression upper bound; zstd-19 usually beats this by 10-20%%)\n",
                   sample.size()/1e6, 100.0*bound/sample.size());
        }
    }
    return 0;
}
