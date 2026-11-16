#ifndef METABULI_REFORMATDB_H
#define METABULI_REFORMATDB_H

#include <string>

// In-place re-encode of one k-mer list directory (diffIdx/split/db.parameters):
// format 0 (15-bit fragment varint) -> LEB128 (dedup=false, format 1) or
// distinct-value LEB128 plus a same-value bitmap (dedup=true, format 2).
// Auto-select keeps format 0 when the re-encoding would not be smaller.
// Used by reformat-db and by build-time re-encoding (--diffidx-format).
// Returns 0 on success (including "kept format 0").
int reformatKmerList(const std::string & dir, int threads, bool dedup);

#endif // METABULI_REFORMATDB_H
