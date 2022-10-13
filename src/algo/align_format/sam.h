#ifndef __SAM_H
#define __SAM_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SAM_VERSION "1.6"

void sam_dump_hdr_hd(FILE* out);

void sam_dump_hdr_sq(const char* name, const int size, FILE* out);

void sam_dump_hdr_pg(const char* pg_id, const char* pg_name, const char* pg_version, int argc, char* argv[], FILE* out);

void sam_dump_hdr_rg(const char* rg_id, const char* sample_name, FILE* out);

void
dump_one_sam_result(
    const u8* qseq,
    const int qid,
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const int sid,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_quality,
    const int map_score,
    const double perc_identity,
    const double eff_perc_identity,
    const int is_complete_map,
    const BOOL dump_md,
    const char* rg_sample,
    kstring_t* out);

#ifdef __cplusplus
}
#endif

#endif // __SAM_H