#include "LocalParameters.h"
#include <Command.h>
#include "FileUtil.h"
#include "DeltaIdxReader.h"
#include "Kmer.h"
#include "common.h"
#include "BitManipulateMacros.h"
#include <cstdint>
#include <cstdio>
#include <vector>
#include <string>
#include <algorithm>
#include <dirent.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>

using namespace std;

// Re-encode a diffIdx sub-DB from the 15-bit fragment varint (format 0) to
// LEB128 byte varint (format 1), without rebuilding. The k-mer DELTAS are
// transcoded (no running sum / seed needed); info/taxID_list/taxonomyDB are
// untouched; the split file's diffIdxOffset is recomputed as a byte position in
// the new stream. Auto-select: if LEB128 would not be smaller (sparse DBs), the
// k-mer list is left in format 0 and nothing is rewritten.
//
// Parallel: the split file's anchors are k-mer boundaries, so the stream is
// partitioned into byte ranges at anchors and each range is transcoded
// independently (pread/pwrite at explicit offsets). Byte-identical to a serial
// run because transcoding is per-delta and chunks concatenate in order.
//
// Usage: metabuli reformat-db <db dir> [--sub-dbs <name>] [--out-sub-db <name>]

// Append one LEB128-encoded delta (7-bit groups, high first, 0x80 on last).
static void appendLEB(uint64_t delta, vector<uint8_t> & out) {
    uint8_t tmp[10];
    int idx = 9;
    tmp[9] = 0x80U | (uint8_t)(delta & 0x7FU); // last byte, end flag
    delta >>= 7;
    while (delta) {
        --idx;
        tmp[idx] = (uint8_t)(delta & 0x7FU);
        delta >>= 7;
    }
    out.insert(out.end(), tmp + idx, tmp + 10);
}

// A k-mer boundary in the 15-bit stream: byte position, k-mer ordinal.
struct Anchor { size_t byteStart; uint64_t ordinal; };

// Dedup-mode (format 2) bitmap output of one chunk. Interior full words are
// pwritten directly; words shared with a neighbor chunk (the chunk's bit range
// does not start/end on a 64-bit boundary) are handed back for an OR-merge.
struct ChunkBits {
    int fd = -1; // < 0: bitmap disabled (sizing pass or LEB mode)
    bool hasHead = false, hasTail = false;
    uint64_t headIdx = 0, headWord = 0;
    uint64_t tailIdx = 0, tailWord = 0;
};

// Transcode input bytes [byteStart, byteEnd) (which begin at ordinal ordStart)
// from 15-bit to LEB128. In dedup mode, zero deltas (identical k-mers) emit no
// LEB bytes; instead every entry contributes one bit (1 = same value as the
// previous entry) to the dupBits bitmap at its global ordinal. For each anchor
// ordinal in chunkOrds (sorted, within this range) record the LEB byte offset
// relative to the chunk start in relOut. fdOut < 0: measure sizes only.
// Otherwise pwrite the LEB bytes at baseOffset.
static void processChunk(int fdIn, size_t byteStart, size_t byteEnd,
                         uint64_t ordStart, const vector<uint64_t> & chunkOrds,
                         vector<uint64_t> & relOut,
                         uint64_t & lebSizeOut, uint64_t & kmerCntOut,
                         int fdOut, size_t baseOffset,
                         bool dedup, ChunkBits & bits) {
    const size_t RBUF = 1u << 22;
    vector<uint8_t> rbuf(RBUF);
    vector<uint8_t> wbuf;
    if (fdOut >= 0) wbuf.reserve(1u << 23);
    uint64_t ordinal = ordStart, lebSize = 0, kmers = 0;
    size_t ti = 0, wflushed = 0, pos = byteStart;
    uint64_t delta = 0; uint32_t frags = 0;
    // Bitmap state (dedup write pass only).
    const bool emitBits = dedup && bits.fd >= 0;
    const uint64_t headIdx = ordStart / 64;
    const bool headShared = (ordStart % 64) != 0;
    uint64_t bitWord = 0;
    uint64_t wordIdx = headIdx;
    vector<uint64_t> bitBuf;
    uint64_t bitBufStartWord = 0;
    bool bitBufActive = false;
    if (emitBits) bitBuf.reserve(1u << 19);
    auto flushBitBuf = [&]() {
        if (!bitBuf.empty()) {
            pwrite(bits.fd, bitBuf.data(), bitBuf.size() * 8, (off_t)(bitBufStartWord * 8));
            bitBuf.clear();
        }
        bitBufActive = false;
    };
    auto emitWord = [&](uint64_t idx, uint64_t w) {
        if (idx == headIdx && headShared) {
            bits.hasHead = true; bits.headIdx = idx; bits.headWord = w;
            return;
        }
        if (!bitBufActive) { bitBufStartWord = idx; bitBufActive = true; }
        bitBuf.push_back(w);
        if (bitBuf.size() >= (1u << 19)) flushBitBuf();
    };
    while (pos < byteEnd) {
        size_t want = min(RBUF, byteEnd - pos);
        ssize_t got = pread(fdIn, rbuf.data(), want, pos);
        if (got <= 0) break;
        if (got & 1) --got;              // never split a 2-byte fragment
        pos += got;
        size_t nfrag = (size_t)got / 2;
        const uint16_t * w16 = reinterpret_cast<const uint16_t *>(rbuf.data());
        for (size_t i = 0; i < nfrag; ++i) {
            if (frags == 0 && ti < chunkOrds.size() && chunkOrds[ti] == ordinal) {
                relOut[ti] = lebSize;
                ++ti;
            }
            uint16_t w = w16[i];
            ++frags;
            delta = (delta << 15) | (uint16_t)(w & 0x7FFF);
            if (w & 0x8000) {
                const bool skipDelta = dedup && delta == 0;
                if (emitBits) {
                    if (skipDelta) bitWord |= 1ULL << (ordinal % 64);
                    if ((ordinal % 64) == 63) {
                        emitWord(wordIdx, bitWord);
                        bitWord = 0; ++wordIdx;
                    }
                }
                if (!skipDelta) {
                    if (fdOut >= 0) {
                        size_t before = wbuf.size();
                        appendLEB(delta, wbuf);
                        lebSize += wbuf.size() - before;
                        if (wbuf.size() >= (1u << 22)) {
                            pwrite(fdOut, wbuf.data(), wbuf.size(), baseOffset + wflushed);
                            wflushed += wbuf.size();
                            wbuf.clear();
                        }
                    } else {
                        int bits2 = 0; uint64_t d = delta; while (d) { ++bits2; d >>= 1; }
                        int leb = (bits2 + 6) / 7; if (leb < 1) leb = 1;
                        lebSize += leb;
                    }
                }
                ++ordinal; ++kmers;
                delta = 0; frags = 0;
            }
        }
    }
    if (fdOut >= 0 && !wbuf.empty()) {
        pwrite(fdOut, wbuf.data(), wbuf.size(), baseOffset + wflushed);
    }
    if (emitBits) {
        flushBitBuf();
        if ((ordinal % 64) != 0) { // partial last word, shared with the next chunk
            bits.hasTail = true; bits.tailIdx = wordIdx; bits.tailWord = bitWord;
        }
    }
    lebSizeOut = lebSize; kmerCntOut = kmers;
}

// Copy every regular file in srcDir to dstDir except those in 'skip'.
static void copyDirExcept(const string & srcDir, const string & dstDir,
                          const vector<string> & skip) {
    DIR * dp = opendir(srcDir.c_str());
    if (!dp) return;
    struct dirent * e;
    while ((e = readdir(dp)) != nullptr) {
        const string name = e->d_name;
        if (name == "." || name == "..") continue;
        bool skipIt = false;
        for (const string & s : skip) if (name == s) { skipIt = true; break; }
        if (skipIt) continue;
        const string src = srcDir + "/" + name;
        if (!FileUtil::fileExists(src.c_str()) || FileUtil::directoryExists(src.c_str())) continue;
        FileUtil::copyFile(src, dstDir + "/" + name);
    }
    closedir(dp);
}

// Re-encode the diffIdx at srcDir into LEB128 (format 1) or, with dedup, into
// distinct-value LEB128 plus a same-value bitmap (format 2) at dstDir,
// recomputing split byte offsets. srcDir == dstDir converts in place (temp +
// rename). threads > 1 transcodes anchor-aligned byte ranges in parallel.
// Returns 0 on success (including "kept format 0"), 1 on error.
static int reformatOne(const string & srcDir, const string & dstDir, int threads,
                       bool dedup) {
    const bool inPlace = (srcDir == dstDir);
    const string srcDiffIdx = srcDir + "/diffIdx";
    const string srcSplit   = srcDir + "/split";
    if (!FileUtil::fileExists(srcDiffIdx.c_str())) {
        cout << "  skip " << srcDir << " (no diffIdx)" << endl;
        return 0;
    }
    const uint8_t srcFmt = DeltaIdxReader::loadDiffidxFormat(srcDiffIdx);
    if (srcFmt != 0) {
        if (dedup) {
            cout << "  " << srcDir << " is format " << (int)srcFmt
                 << "; --dedup needs the format 0 (15-bit) source. Skipping." << endl;
        } else {
            cout << "  " << srcDir << " already format " << (int)srcFmt << "; skipping" << endl;
        }
        return 0;
    }

    // Load split records (kept verbatim; only diffIdxOffset is rewritten).
    // In format 0, diffIdxOffset is a uint16 fragment count -> byte = x2.
    vector<DiffIdxSplit> splits;
    {
        FILE * sf = fopen(srcSplit.c_str(), "rb");
        if (!sf) { cout << "    cannot open split" << endl; return 1; }
        DiffIdxSplit rec;
        while (fread(&rec, sizeof(DiffIdxSplit), 1, sf) == 1) splits.push_back(rec);
        fclose(sf);
    }

    int fdIn = open(srcDiffIdx.c_str(), O_RDONLY);
    if (fdIn < 0) { cout << "    cannot open diffIdx" << endl; return 1; }
    struct stat st;
    fstat(fdIn, &st);
    const size_t fileSize = (size_t)st.st_size;

    // Anchors = k-mer boundaries at which a thread can start decoding: the
    // implicit stream start plus every distinct split boundary. A split's
    // diffIdxOffset addresses the byte where the k-mer at ordinal infoIdxOffset
    // begins; ordinal 0 maps to byte 0 (covers the zero-padding split records).
    vector<Anchor> anchors;
    anchors.push_back({0, 0});
    for (const DiffIdxSplit & s : splits) {
        if (s.infoIdxOffset != 0) {
            anchors.push_back({s.diffIdxOffset * 2, s.infoIdxOffset});
        }
    }
    sort(anchors.begin(), anchors.end(),
         [](const Anchor & a, const Anchor & b) { return a.ordinal < b.ordinal; });
    anchors.erase(unique(anchors.begin(), anchors.end(),
                         [](const Anchor & a, const Anchor & b) { return a.ordinal == b.ordinal; }),
                  anchors.end());
    const size_t nAnchor = anchors.size();

    // Partition anchors into chunk boundaries spread evenly by byte position.
    size_t T = (size_t)max(1, threads);
    if (T > nAnchor) T = nAnchor;
    vector<size_t> bnd; bnd.push_back(0);
    for (size_t t = 1; t < T; ++t) {
        size_t targetByte = (size_t)((double)fileSize * t / (double)T);
        size_t idx = bnd.back() + 1;
        while (idx < nAnchor && anchors[idx].byteStart < targetByte) ++idx;
        if (idx <= bnd.back()) idx = bnd.back() + 1;
        if (idx >= nAnchor) break;
        bnd.push_back(idx);
    }
    const size_t nChunk = bnd.size();

    // Per-chunk inputs and outputs.
    vector<size_t>   chByteStart(nChunk), chByteEnd(nChunk);
    vector<uint64_t> chOrdStart(nChunk);
    vector<vector<uint64_t>> chOrds(nChunk), chRel(nChunk);
    vector<uint64_t> chLeb(nChunk, 0), chKmer(nChunk, 0);
    for (size_t c = 0; c < nChunk; ++c) {
        const size_t a0 = bnd[c];
        const size_t a1 = (c + 1 < nChunk) ? bnd[c + 1] : nAnchor;
        chByteStart[c] = anchors[a0].byteStart;
        chByteEnd[c]   = (a1 < nAnchor) ? anchors[a1].byteStart : fileSize;
        chOrdStart[c]  = anchors[a0].ordinal;
        for (size_t a = a0; a < a1; ++a) chOrds[c].push_back(anchors[a].ordinal);
        chRel[c].assign(chOrds[c].size(), 0);
    }

    // Pass A: size every chunk in parallel (no output).
#pragma omp parallel for schedule(dynamic) num_threads((int)T)
    for (size_t c = 0; c < nChunk; ++c) {
        ChunkBits noBits;
        processChunk(fdIn, chByteStart[c], chByteEnd[c], chOrdStart[c],
                     chOrds[c], chRel[c], chLeb[c], chKmer[c], -1, 0, dedup, noBits);
    }
    uint64_t lebBytes = 0, kmerCnt = 0;
    for (size_t c = 0; c < nChunk; ++c) { lebBytes += chLeb[c]; kmerCnt += chKmer[c]; }
    const uint64_t bitWords = dedup ? (kmerCnt + 63) / 64 : 0;
    const uint64_t newBytes = lebBytes + bitWords * 8;

    if (dedup) {
        cout << "  " << srcDir << ": " << kmerCnt << " k-mers, 15-bit "
             << (fileSize/1e9) << " GB -> dedup LEB128 " << (lebBytes/1e9)
             << " GB + bitmap " << (bitWords*8/1e9) << " GB ("
             << (100.0*newBytes/fileSize) << "%)" << endl;
    } else {
        cout << "  " << srcDir << ": " << kmerCnt << " k-mers, 15-bit "
             << (fileSize/1e9) << " GB -> LEB128 " << (lebBytes/1e9) << " GB ("
             << (100.0*lebBytes/fileSize) << "%)" << endl;
    }

    if (!inPlace && !FileUtil::directoryExists(dstDir.c_str())) {
        FileUtil::makeDir(dstDir.c_str());
    }
    if (newBytes >= fileSize) {
        close(fdIn);
        if (inPlace) {
            cout << "    not smaller; keeping format 0." << endl;
            return 0;
        }
        copyDirExcept(srcDir, dstDir, {});
        cout << "    not smaller; copied verbatim as format 0." << endl;
        return 0;
    }

    // Absolute byte offset of each chunk in the output, then of each anchor.
    vector<size_t> chBase(nChunk, 0);
    for (size_t c = 1; c < nChunk; ++c) chBase[c] = chBase[c - 1] + chLeb[c - 1];
    vector<uint64_t> anchorOffset(nAnchor, 0);
    for (size_t c = 0; c < nChunk; ++c) {
        const size_t a0 = bnd[c];
        for (size_t j = 0; j < chOrds[c].size(); ++j) {
            anchorOffset[a0 + j] = chBase[c] + chRel[c][j];
        }
    }

    // New offset for each split record: look up its ordinal among the anchors.
    for (DiffIdxSplit & s : splits) {
        size_t lo = 0, hi = nAnchor; // lower_bound on ordinal
        while (lo < hi) { size_t m = (lo + hi) / 2; if (anchors[m].ordinal < s.infoIdxOffset) lo = m + 1; else hi = m; }
        s.diffIdxOffset = (lo < nAnchor && anchors[lo].ordinal == s.infoIdxOffset)
                              ? anchorOffset[lo] : 0;
    }

    // Pass B: transcode + write every chunk in parallel at its base offset.
    const string newDiffIdx = inPlace ? (srcDiffIdx + ".leb.tmp") : (dstDir + "/diffIdx");
    const string dstBits = dstDir + "/dupBits";
    const string newBits = inPlace ? (dstBits + ".tmp") : dstBits;
    int fdOut = open(newDiffIdx.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
    if (fdOut < 0) { close(fdIn); cout << "    cannot open output" << endl; return 1; }
    if (ftruncate(fdOut, (off_t)lebBytes) != 0) { /* best effort */ }
    int fdBits = -1;
    if (dedup) {
        fdBits = open(newBits.c_str(), O_RDWR | O_CREAT | O_TRUNC, 0644);
        if (fdBits < 0) { close(fdIn); close(fdOut); cout << "    cannot open dupBits" << endl; return 1; }
        if (ftruncate(fdBits, (off_t)(bitWords * 8)) != 0) { /* best effort */ }
    }
    vector<ChunkBits> chBits(nChunk);
#pragma omp parallel for schedule(dynamic) num_threads((int)T)
    for (size_t c = 0; c < nChunk; ++c) {
        uint64_t leb = 0, km = 0;
        vector<uint64_t> dummy(chOrds[c].size(), 0);
        chBits[c].fd = fdBits;
        processChunk(fdIn, chByteStart[c], chByteEnd[c], chOrdStart[c],
                     chOrds[c], dummy, leb, km, fdOut, chBase[c], dedup, chBits[c]);
    }
    if (dedup) {
        // OR-merge bitmap words shared between neighboring chunks, then write.
        vector<pair<uint64_t, uint64_t>> shared; // (wordIdx, bits)
        for (const ChunkBits & b : chBits) {
            if (b.hasHead) shared.push_back({b.headIdx, b.headWord});
            if (b.hasTail) shared.push_back({b.tailIdx, b.tailWord});
        }
        sort(shared.begin(), shared.end());
        for (size_t i = 0; i < shared.size();) {
            uint64_t w = shared[i].second;
            size_t j = i + 1;
            while (j < shared.size() && shared[j].first == shared[i].first) {
                w |= shared[j].second;
                ++j;
            }
            pwrite(fdBits, &w, 8, (off_t)(shared[i].first * 8));
            i = j;
        }
        close(fdBits);
    }
    close(fdOut);
    close(fdIn);

    if (inPlace) {
        if (rename(newDiffIdx.c_str(), srcDiffIdx.c_str()) != 0) {
            cout << "    rename failed" << endl; return 1;
        }
        if (dedup && rename(newBits.c_str(), dstBits.c_str()) != 0) {
            cout << "    dupBits rename failed" << endl; return 1;
        }
    } else {
        copyDirExcept(srcDir, dstDir, {"diffIdx", "split"});
    }
    {
        FILE * sf = fopen((dstDir + "/split").c_str(), "wb");
        fwrite(splits.data(), sizeof(DiffIdxSplit), splits.size(), sf);
        fclose(sf);
    }
    appendDbParameter(dstDir, "Diffidx_format", dedup ? "2" : "1");
    if (dedup) {
        cout << "    reformatted to dedup format 2 (" << (lebBytes/1e9)
             << " GB deltas + " << (bitWords*8/1e9) << " GB bitmap)"
             << (inPlace ? "" : " at " + dstDir) << "." << endl;
    } else {
        cout << "    reformatted to LEB128 (" << (lebBytes/1e9) << " GB)"
             << (inPlace ? "" : " at " + dstDir) << "." << endl;
    }
    return 0;
}

int reformatDb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string dbDir = par.filenames[0];
    // --sub-dbs selects the source sub-DB(s); for --out-sub-db it must name one.
    const string subName = par.subDbList;
    const string outName = par.outSubDb;
    const int threads = par.threads;
    const string listsDir = dbDir + "/" + SUBDB_DIR_NAME;

    // --out-sub-db: write the result to a new sub-DB, keep the source intact.
    if (!outName.empty()) {
        if (outName == SUBDB_BASE_NAME) {
            cerr << "Error: --out-sub-db cannot be \"" << SUBDB_BASE_NAME
                 << "\" (reserved for the top-level list)." << endl;
            return 1;
        }
        if (subName.find(',') != string::npos) {
            cerr << "Error: --out-sub-db needs a single source sub-DB; "
                 << "--sub-dbs \"" << subName << "\" names more than one." << endl;
            return 1;
        }
        vector<SubDbEntry> srcs = resolveSubDbs(dbDir, subName);
        if (srcs.size() != 1) {
            cerr << "Error: --out-sub-db needs exactly one source sub-DB; found "
                 << srcs.size() << ". Select one with --sub-dbs." << endl;
            return 1;
        }
        const string & srcDir = srcs[0].dir;
        if (!FileUtil::directoryExists(listsDir.c_str())) {
            FileUtil::makeDir(listsDir.c_str());
        }
        const string dstDir = listsDir + "/" + outName;
        if (FileUtil::directoryExists(dstDir.c_str())) {
            cerr << "Error: sub-DB \"" << outName << "\" already exists at "
                 << dstDir << "; remove it or choose another name." << endl;
            return 1;
        }
        cout << "Reformatting " << srcDir << " -> new sub-DB \"" << outName << "\"" << endl;
        int rc = reformatOne(srcDir, dstDir, threads, par.dedupKmers != 0);
        if (rc == 0) {
            cout << "Test with:  classify ... " << dbDir << " --sub-dbs " << outName
                 << "   (baseline: --sub-dbs " << srcs[0].name << ")" << endl;
        }
        return rc;
    }

    // In-place conversion of the selected sub-DB(s) (default: all; flat -> base).
    vector<SubDbEntry> subs = resolveSubDbs(dbDir, subName);
    int rc = 0;
    for (const SubDbEntry & e : subs) {
        if (reformatOne(e.dir, e.dir, threads, par.dedupKmers != 0) != 0) rc = 1;
    }
    return rc;
}
