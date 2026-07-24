#ifndef METABULI_INTERLEAVED_KSEQ_WRAPPER_H
#define METABULI_INTERLEAVED_KSEQ_WRAPPER_H

#include "KSeqWrapper.h"
#include "BamKSeqWrapper.h"
#include <string>

// Presents one interleaved paired-end file (records R1,R2,R1,R2,...) as a
// single-mate stream, so the existing two-file paired read loops (kseq1/kseq2)
// work unchanged:
//   - mate 0 exposes records 0, 2, 4, ...  (the R1 of each pair)
//   - mate 1 exposes records 1, 3, 5, ...  (the R2 of each pair)
// One ReadEntry() advances by one *pair* (it consumes both physical records),
// so the pair index lines up between the two instances. Each instance keeps
// its own underlying handle on the same file, so the file is parsed twice; the
// exposed entry is always the last physical record read, so no deep copy is
// needed (kseq reuses its buffers only on the next read).
class InterleavedKSeqWrapper : public KSeqWrapper {
public:
    InterleavedKSeqWrapper(const char *file, int mate)
        : inner(KSeqFactory(file)), mate(mate), started(false) {
        if (inner != nullptr) {
            type = inner->type;
        }
    }

    bool ReadEntry() override {
        if (inner == nullptr) {
            return false;
        }
        if (mate == 0) {
            if (!started) {
                started = true;
                if (!inner->ReadEntry()) { return false; } // record 0
            } else {
                if (!inner->ReadEntry()) { return false; } // skip previous pair's R2
                if (!inner->ReadEntry()) { return false; } // this pair's R1
            }
        } else { // mate == 1
            if (!inner->ReadEntry()) { return false; }     // skip this pair's R1
            if (!inner->ReadEntry()) { return false; }     // this pair's R2
        }
        // Expose the record just read; inner is not advanced again before the
        // caller consumes it, so a shallow copy of the entry is safe.
        this->entry = inner->entry;
        return true;
    }

    ~InterleavedKSeqWrapper() override {
        delete inner;
    }

private:
    KSeqWrapper *inner;
    int mate;      // 0 -> R1 stream, 1 -> R2 stream
    bool started;  // mate-0 only: first ReadEntry consumes a single record
};

// Factory for query readers. When interleaved, both mates are drawn from
// file0 (file1 is ignored); otherwise the requested mate's own file is used.
inline bool hasBamExtension(const std::string &path) {
    const std::string suffix = ".bam";
    return path.size() >= suffix.size()
        && path.compare(path.size() - suffix.size(), suffix.size(), suffix) == 0;
}

inline KSeqWrapper *createQueryKseqWrapper(const std::string &file0,
                                           const std::string &file1,
                                           bool interleaved,
                                           int mate) {
    const std::string &path = (mate == 0) ? file0 : file1;
    // BAM input is single-end: read sequences straight out of the alignment
    // records (see BamKSeqWrapper). Interleaved is not combined with BAM.
    if (hasBamExtension(path)) {
        return new BamKSeqWrapper(path.c_str());
    }
    if (interleaved) {
        return new InterleavedKSeqWrapper(file0.c_str(), mate);
    }
    return KSeqFactory(path.c_str());
}

#endif // METABULI_INTERLEAVED_KSEQ_WRAPPER_H
