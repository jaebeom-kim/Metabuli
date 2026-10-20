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
#include <dirent.h>

using namespace std;

// Re-encode a diffIdx sub-DB from the 15-bit fragment varint (format 0) to
// LEB128 byte varint (format 1), without rebuilding. The k-mer sequence is
// decoded and re-encoded; info/taxID_list/taxonomyDB are untouched; the split
// file's diffIdxOffset is recomputed as a byte position in the new stream.
// Auto-select: if LEB128 would not be smaller (sparse DBs), the k-mer list is
// left in format 0 and nothing is rewritten.
//
// Usage: metabuli reformat-db <db dir> [--sub-dbs <name>] [--out-sub-db <name>]
//   Operates on <db dir>/kmer-lists/<name>/ (default: every sub-DB), or on a
//   legacy flat <db dir>/diffIdx. --sub-dbs selects which sub-DB(s) to convert.
//   With --out-sub-db, writes the result to kmer-lists/<name>/ as a new sub-DB
//   and leaves the source untouched (for A/B testing via classify --sub-dbs).

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

// One pass over the diffIdx: measure both encodings' sizes (bytes).
struct SizeProbe { uint64_t curBytes = 0; uint64_t lebBytes = 0; uint64_t kmerCnt = 0; };

static SizeProbe probeSizes(const string & diffIdxFile) {
    SizeProbe s;
    FILE* fp = fopen(diffIdxFile.c_str(), "rb");
    const size_t BUF = 1u << 22;
    vector<uint16_t> buf(BUF);
    uint64_t delta = 0; uint32_t frags = 0; size_t n;
    while ((n = fread(buf.data(), sizeof(uint16_t), BUF, fp)) > 0) {
        for (size_t i = 0; i < n; ++i) {
            uint16_t w = buf[i];
            ++frags;
            delta = (delta << 15) | (uint16_t)(w & 0x7FFF);
            if (w & 0x8000) {
                s.curBytes += (uint64_t)frags * 2;
                int bits = 0; uint64_t d = delta; while (d) { ++bits; d >>= 1; }
                int leb = (bits + 6) / 7; if (leb < 1) leb = 1;
                s.lebBytes += leb;
                ++s.kmerCnt;
                delta = 0; frags = 0;
            }
        }
    }
    fclose(fp);
    return s;
}

// Copy every regular file in srcDir to dstDir except those in 'skip'.
static void copyDirExcept(const string & srcDir, const string & dstDir,
                          const vector<string> & skip) {
    DIR* dp = opendir(srcDir.c_str());
    if (!dp) return;
    struct dirent* e;
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

// Re-encode the diffIdx at srcDir into LEB128 at dstDir, recomputing split
// byte offsets. srcDir == dstDir performs an in-place conversion (temp + rename).
// Returns 0 on success (including "kept format 0"), 1 on error.
static int reformatOne(const string & srcDir, const string & dstDir) {
    const bool inPlace = (srcDir == dstDir);
    const string srcDiffIdx = srcDir + "/diffIdx";
    const string srcSplit   = srcDir + "/split";
    if (!FileUtil::fileExists(srcDiffIdx.c_str())) {
        cout << "  skip " << srcDir << " (no diffIdx)" << endl;
        return 0;
    }
    if (DeltaIdxReader::loadDiffidxFormat(srcDiffIdx) == 1) {
        cout << "  " << srcDir << " already format 1 (LEB128); skipping" << endl;
        return 0;
    }

    SizeProbe s = probeSizes(srcDiffIdx);
    cout << "  " << srcDir << ": " << s.kmerCnt << " k-mers, 15-bit "
         << (s.curBytes/1e9) << " GB -> LEB128 " << (s.lebBytes/1e9) << " GB ("
         << (100.0*s.lebBytes/s.curBytes) << "%)" << endl;

    if (!inPlace && !FileUtil::directoryExists(dstDir.c_str())) {
        FileUtil::makeDir(dstDir.c_str());
    }

    if (s.lebBytes >= s.curBytes) {
        if (inPlace) {
            cout << "    LEB128 not smaller; keeping format 0." << endl;
            return 0;
        }
        // Out-of-place but not worth converting: copy the sub-DB verbatim so the
        // requested new sub-DB still exists (stays format 0).
        copyDirExcept(srcDir, dstDir, {});
        cout << "    LEB128 not smaller; copied verbatim as format 0." << endl;
        return 0;
    }

    // Load split anchors (ADkmer, old diffIdxOffset [uint16 count], infoIdxOffset).
    vector<DiffIdxSplit> splits;
    {
        FILE* sf = fopen(srcSplit.c_str(), "rb");
        if (!sf) { cout << "    cannot open split" << endl; return 1; }
        DiffIdxSplit rec;
        while (fread(&rec, sizeof(DiffIdxSplit), 1, sf) == 1) splits.push_back(rec);
        fclose(sf);
    }

    // A split's diffIdxOffset is the byte where decoding RESUMES: setReadPosition
    // seeds lastValue = ADkmer (the k-mer at ordinal infoIdxOffset-1) and then
    // reads the delta at diffIdxOffset, producing the k-mer at ordinal
    // infoIdxOffset. So the new byte offset is captured just before encoding the
    // delta of the k-mer at ordinal == infoIdxOffset.
    vector<size_t> splitResumeOrdinal(splits.size());
    for (size_t k = 0; k < splits.size(); ++k) {
        splitResumeOrdinal[k] = splits[k].infoIdxOffset;
    }

    // Re-encode: decode 15-bit, emit LEB128, tracking byte position and k-mer
    // ordinal to rewrite each split's diffIdxOffset.
    const string newDiffIdx = inPlace ? (srcDiffIdx + ".leb.tmp") : (dstDir + "/diffIdx");
    FILE* out = fopen(newDiffIdx.c_str(), "wb");
    if (!out) { cout << "    cannot open output" << endl; return 1; }
    FILE* in = fopen(srcDiffIdx.c_str(), "rb");
    const size_t BUF = 1u << 22;
    vector<uint16_t> rbuf(BUF);
    vector<uint8_t> wbuf; wbuf.reserve(1u << 23);

    uint64_t delta = 0; uint32_t frags = 0;
    uint64_t ordinal = 0;      // index of the k-mer whose delta begins next
    uint64_t outBytes = 0;     // bytes written so far (byte offset of next delta)
    size_t nread;
    while ((nread = fread(rbuf.data(), sizeof(uint16_t), BUF, in)) > 0) {
        for (size_t i = 0; i < nread; ++i) {
            if (frags == 0) {
                for (size_t k = 0; k < splits.size(); ++k) {
                    if (splitResumeOrdinal[k] == ordinal) {
                        splits[k].diffIdxOffset = outBytes;
                    }
                }
            }
            uint16_t w = rbuf[i];
            ++frags;
            delta = (delta << 15) | (uint16_t)(w & 0x7FFF);
            if (w & 0x8000) {
                size_t before = wbuf.size();
                appendLEB(delta, wbuf);
                outBytes += wbuf.size() - before;
                ++ordinal;
                delta = 0; frags = 0;
                if (wbuf.size() >= (1u << 22)) {
                    fwrite(wbuf.data(), 1, wbuf.size(), out);
                    wbuf.clear();
                }
            }
        }
    }
    if (!wbuf.empty()) fwrite(wbuf.data(), 1, wbuf.size(), out);
    fclose(in);
    fclose(out);

    if (inPlace) {
        if (rename(newDiffIdx.c_str(), srcDiffIdx.c_str()) != 0) {
            cout << "    rename failed" << endl; return 1;
        }
    } else {
        // Copy the rest of the sub-DB (info, taxID_list, db.parameters, ...);
        // diffIdx was written above, split is written below.
        copyDirExcept(srcDir, dstDir, {"diffIdx", "split"});
    }

    // Write split (new byte offsets) and stamp the format.
    {
        FILE* sf = fopen((dstDir + "/split").c_str(), "wb");
        fwrite(splits.data(), sizeof(DiffIdxSplit), splits.size(), sf);
        fclose(sf);
    }
    appendDbParameter(dstDir, "Diffidx_format", "1");
    cout << "    reformatted to LEB128 (" << (outBytes/1e9) << " GB)"
         << (inPlace ? "" : " at " + dstDir) << "." << endl;
    return 0;
}

int reformatDb(int argc, const char **argv, const Command &command) {
    LocalParameters &par = LocalParameters::getLocalInstance();
    par.parseParameters(argc, argv, command, false, Parameters::PARSE_ALLOW_EMPTY, 0);
    const string dbDir = par.filenames[0];
    // --sub-dbs selects the source sub-DB(s); for --out-sub-db it must name one.
    const string subName = par.subDbList;
    const string outName = par.outSubDb;
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
        // Resolve the single source sub-DB directory (flat DB -> "base" at dbDir).
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
        int rc = reformatOne(srcDir, dstDir);
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
        if (reformatOne(e.dir, e.dir) != 0) rc = 1;
    }
    return rc;
}
