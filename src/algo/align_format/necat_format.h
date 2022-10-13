#ifndef __NECAT_FORMATH_H
#define __NECAT_FORMATH_H

#include "paf.h"
#include "sam.h"
#include "../../corelib/seqdb.h"
#include "../../corelib/raw_reads_reader.h"

#ifdef __cplusplus
extern "C" {
#endif

void 
dump_sam_hdr(const char* rg_id,
    const char* rg_sample,
    const CSeqDB* db,
    RawReadReader* rdb,
    int argc,
    char* argv[],
    FILE* out);

#ifdef __cplusplus
}
#endif

#endif // __NECAT_FORMATH_H