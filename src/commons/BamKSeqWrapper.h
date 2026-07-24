#ifndef METABULI_BAM_KSEQ_WRAPPER_H
#define METABULI_BAM_KSEQ_WRAPPER_H

#include "KSeqWrapper.h"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <string>
#include <vector>

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

// Presents the reads stored in a BAM file as a single-end KSeqWrapper stream,
// so the existing query read loops can consume a BAM exactly like a FASTA/FASTQ
// file. BAM is BGZF (concatenated gzip blocks), which zlib's gz* API inflates
// transparently; we then parse the BAM binary records ourselves. Only the read
// name, SEQ and QUAL are needed (CIGAR and aux fields are skipped by advancing
// over them). Records that are secondary (FLAG 0x100), supplementary (0x800),
// QC-fail (0x200) or carry no stored sequence are skipped, so each read is
// emitted once; reverse-strand reads (0x10) are reverse-complemented back to
// their original sequenced orientation.
//
// Assumes a little-endian host (BAM is little-endian); Metabuli's targets are.
class BamKSeqWrapper : public KSeqWrapper {
public:
    explicit BamKSeqWrapper(const char *fileName) {
        type = KSeqWrapper::KSEQ_GZIP;
        entry.name = {0, 0, nullptr};
        entry.comment = {0, 0, nullptr};
        entry.sequence = {0, 0, nullptr};
        entry.qual = {0, 0, nullptr};
#ifdef HAVE_ZLIB
        file = gzopen(fileName, "rb");
        ok = (file != nullptr) && readHeader();
#else
        (void) fileName;
        ok = false;
#endif
    }

    bool ReadEntry() override {
#ifdef HAVE_ZLIB
        if (!ok) {
            return false;
        }
        while (true) {
            int32_t blockSize;
            if (!readFull(&blockSize, sizeof(blockSize))) {
                return false; // clean EOF (or truncated trailing record)
            }
            if (blockSize < 32) {
                return false; // malformed
            }
            rec.resize(static_cast<size_t>(blockSize));
            if (!readFull(rec.data(), static_cast<size_t>(blockSize))) {
                return false;
            }

            const uint8_t *p = reinterpret_cast<const uint8_t *>(rec.data());
            const uint8_t  l_read_name = p[8];
            const uint16_t n_cigar_op  = readU16(p + 12);
            const uint16_t flag        = readU16(p + 14);
            const uint32_t l_seq       = readU32(p + 16);

            // Emit each read once, and only when it carries a sequence.
            if ((flag & 0x100) || (flag & 0x800) || (flag & 0x200) || l_seq == 0 || l_read_name == 0) {
                continue;
            }
            // Guard against a corrupt record whose fields overrun the block.
            const size_t need = 32 + static_cast<size_t>(l_read_name)
                              + 4 * static_cast<size_t>(n_cigar_op)
                              + (static_cast<size_t>(l_seq) + 1) / 2
                              + static_cast<size_t>(l_seq);
            if (need > rec.size()) {
                return false;
            }

            const uint8_t *name = p + 32;
            const uint8_t *seqPacked = name + l_read_name + 4 * static_cast<size_t>(n_cigar_op);
            const uint8_t *qualRaw = seqPacked + (static_cast<size_t>(l_seq) + 1) / 2;

            // name: l_read_name includes the trailing NUL.
            nameBuf.assign(reinterpret_cast<const char *>(name), static_cast<size_t>(l_read_name) - 1);

            // sequence: 4-bit packed, high nibble first.
            static const char code2base[16] =
                {'=', 'A', 'C', 'M', 'G', 'R', 'S', 'V', 'T', 'W', 'Y', 'H', 'K', 'D', 'B', 'N'};
            seqBuf.resize(l_seq);
            for (uint32_t i = 0; i < l_seq; ++i) {
                const uint8_t byte = seqPacked[i >> 1];
                const uint8_t nib = (i & 1u) ? (byte & 0x0Fu) : (byte >> 4);
                seqBuf[i] = code2base[nib];
            }

            // quality: 0xFF in the first byte means "no quality stored".
            const bool hasQual = (qualRaw[0] != 0xFF);
            if (hasQual) {
                qualBuf.resize(l_seq);
                for (uint32_t i = 0; i < l_seq; ++i) {
                    qualBuf[i] = static_cast<char>(qualRaw[i] + 33); // Phred+33
                }
            } else {
                qualBuf.clear();
            }

            if (flag & 0x10) { // reverse strand -> restore original orientation
                reverseComplement(seqBuf);
                if (hasQual) {
                    std::reverse(qualBuf.begin(), qualBuf.end());
                }
            }

            setEntry(hasQual);
            return true;
        }
#else
        return false;
#endif
    }

    ~BamKSeqWrapper() override {
#ifdef HAVE_ZLIB
        if (file != nullptr) {
            gzclose(file);
        }
#endif
    }

private:
#ifdef HAVE_ZLIB
    gzFile file = nullptr;
#endif
    bool ok = false;
    std::vector<char> rec;
    std::string nameBuf, seqBuf, qualBuf;
    std::string emptyBuf; // stable, NUL-terminated backing for the (unused) comment field

    static uint16_t readU16(const uint8_t *p) { uint16_t v; std::memcpy(&v, p, 2); return v; }
    static uint32_t readU32(const uint8_t *p) { uint32_t v; std::memcpy(&v, p, 4); return v; }

#ifdef HAVE_ZLIB
    bool readFull(void *buf, size_t n) {
        size_t got = 0;
        char *out = static_cast<char *>(buf);
        while (got < n) {
            const int r = gzread(file, out + got, static_cast<unsigned>(n - got));
            if (r <= 0) {
                return false;
            }
            got += static_cast<size_t>(r);
        }
        return true;
    }

    bool skip(size_t n) {
        char tmp[4096];
        while (n > 0) {
            const size_t chunk = n < sizeof(tmp) ? n : sizeof(tmp);
            if (!readFull(tmp, chunk)) {
                return false;
            }
            n -= chunk;
        }
        return true;
    }

    bool readHeader() {
        char magic[4];
        if (!readFull(magic, 4) || std::memcmp(magic, "BAM\1", 4) != 0) {
            return false;
        }
        int32_t l_text;
        if (!readFull(&l_text, 4) || l_text < 0 || !skip(static_cast<size_t>(l_text))) {
            return false;
        }
        int32_t n_ref;
        if (!readFull(&n_ref, 4) || n_ref < 0) {
            return false;
        }
        for (int32_t i = 0; i < n_ref; ++i) {
            int32_t l_name;
            if (!readFull(&l_name, 4) || l_name < 0 || !skip(static_cast<size_t>(l_name))) {
                return false;
            }
            int32_t l_ref;
            if (!readFull(&l_ref, 4)) { // reference length, unused
                return false;
            }
        }
        return true;
    }
#endif

    static char complement(char c) {
        switch (c) {
            case 'A': return 'T'; case 'T': return 'A';
            case 'C': return 'G'; case 'G': return 'C';
            case 'a': return 't'; case 't': return 'a';
            case 'c': return 'g'; case 'g': return 'c';
            default:  return 'N';
        }
    }

    static void reverseComplement(std::string &s) {
        const size_t n = s.size();
        for (size_t i = 0; i < n / 2; ++i) {
            const char a = complement(s[i]);
            s[i] = complement(s[n - 1 - i]);
            s[n - 1 - i] = a;
        }
        if (n & 1u) {
            s[n / 2] = complement(s[n / 2]);
        }
    }

    void setEntry(bool hasQual) {
        entry.name.s = &nameBuf[0];        entry.name.l = nameBuf.size();      entry.name.m = nameBuf.size();
        entry.sequence.s = &seqBuf[0];     entry.sequence.l = seqBuf.size();   entry.sequence.m = seqBuf.size();
        entry.comment.s = &emptyBuf[0];    entry.comment.l = 0;                entry.comment.m = 0;
        if (hasQual) {
            entry.qual.s = &qualBuf[0];    entry.qual.l = qualBuf.size();      entry.qual.m = qualBuf.size();
        } else {
            entry.qual.s = &emptyBuf[0];   entry.qual.l = 0;                   entry.qual.m = 0;
        }
    }
};

#endif // METABULI_BAM_KSEQ_WRAPPER_H
