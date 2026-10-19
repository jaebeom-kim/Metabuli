#ifndef METABULI_OFFSET_KSEQ_WRAPPER_H
#define METABULI_OFFSET_KSEQ_WRAPPER_H

#include "KSeqWrapper.h"
#include "FileUtil.h"
#include "Util.h"
#include "Debug.h"

#ifdef HAVE_ZLIB
#include <zlib.h>
#endif

// Query readers with byte-based progress support. The vendored KSeqWrapper
// subclasses keep their file handles private, so the consumed-bytes position
// needs Metabuli-owned readers; createQueryKseqWrapper() returns only
// MetabuliKSeqWrapper instances.
class MetabuliKSeqWrapper : public KSeqWrapper {
public:
    // Bytes consumed from the underlying file so far (compressed bytes for
    // gzip input); 0 when the source has no meaningful offset (e.g., stdin).
    virtual size_t getSourceOffset() = 0;
};

// Adapter over the vendored readers. For a plain file, the entry's
// headerOffset is the byte position of the latest record, which is exact
// enough for progress; other sources (stdin, bz2) report no offset.
class OffsetAdapterKSeq : public MetabuliKSeqWrapper {
public:
    explicit OffsetAdapterKSeq(KSeqWrapper *inner) : inner(inner) {
        seq = nullptr;
        if (inner != nullptr) {
            type = inner->type;
        }
    }

    bool ReadEntry() override {
        if (inner == nullptr || !inner->ReadEntry()) {
            return false;
        }
        // Shallow copy is safe: inner is not advanced again before the caller
        // consumes the entry.
        this->entry = inner->entry;
        return true;
    }

    size_t getSourceOffset() override {
        return (inner != nullptr && inner->type == KSEQ_FILE) ? inner->entry.headerOffset : 0;
    }

    ~OffsetAdapterKSeq() override {
        delete inner;
    }

private:
    KSeqWrapper *inner;
};

#ifdef HAVE_ZLIB
namespace METABULI_KSEQGZ {
    KSEQ_INIT(gzFile, gzread)
}

// Mirror of the vendored KSeqGzip that owns its gzFile, so the compressed
// read position (gzoffset) is available for progress reporting.
class GzipOffsetKSeq : public MetabuliKSeqWrapper {
public:
    explicit GzipOffsetKSeq(const char *fileName) {
        seq = nullptr;
        if (FileUtil::fileExists(fileName) == false) {
            errno = ENOENT;
            perror(fileName);
            EXIT(EXIT_FAILURE);
        }
        file = gzopen(fileName, "r");
        if (file == NULL) {
            perror(fileName);
            EXIT(EXIT_FAILURE);
        }
        kseq = METABULI_KSEQGZ::kseq_init(file);
        type = KSEQ_GZIP;
    }

    bool ReadEntry() override {
        if (METABULI_KSEQGZ::kseq_read(kseq) < 0) {
            return false;
        }
        entry.name = kseq->name;
        entry.comment = kseq->comment;
        entry.sequence = kseq->seq;
        entry.qual = kseq->qual;
        entry.headerOffset = 0;
        entry.sequenceOffset = 0;
        entry.newlineCount = kseq->newlineCount;
        return true;
    }

    size_t getSourceOffset() override {
        z_off_t pos = gzoffset(file); // compressed bytes consumed
        return pos > 0 ? (size_t) pos : 0;
    }

    ~GzipOffsetKSeq() override {
        METABULI_KSEQGZ::kseq_destroy(kseq);
        gzclose(file);
    }

private:
    gzFile file = nullptr;
    METABULI_KSEQGZ::kseq_t *kseq = nullptr;
};
#endif

// Offset-aware replacement for KSeqFactory on query files.
inline MetabuliKSeqWrapper *makeOffsetKseq(const char *file) {
#ifdef HAVE_ZLIB
    if (Util::endsWith(".gz", file)) {
        return new GzipOffsetKSeq(file);
    }
#endif
    return new OffsetAdapterKSeq(KSeqFactory(file));
}

#endif // METABULI_OFFSET_KSEQ_WRAPPER_H
