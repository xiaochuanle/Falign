#ifndef __PAF_H
#define __PAF_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void
print_paf_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, kstring_t* out);

void
dump_one_paf_result(
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_score,
    const BOOL dump_cigar,
    const BOOL dump_md,
    kstring_t* out);

#ifdef __cplusplus
}
#endif

#endif // __PAF_H