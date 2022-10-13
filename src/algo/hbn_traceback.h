#ifndef __HBN_TRACE_BACK_H
#define __HBN_TRACE_BACK_H

#include "dalign.h"
#include "edlib_wrapper.h"
#include "ksw2.h"
#include "ksw2_wrapper.h"
#include "small_edlib_align.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qoff, qend, qsize;
    int soff, send, ssize;
    int score;
    int dist;
    double ident_perc;
    kstring_t qabuf;
    kstring_t sabuf;
    kstring_t ext_qabuf;
    kstring_t ext_sabuf;
    vec_u8 qfrag;
    vec_u8 sfrag;
    char* qas;
    char* qae;
    char* sas;
    char* sae;
    EdlibAlignData* edlib;
    Ksw2Data* ksw;
    DalignData* dalign;
    small_edlib_align_struct* small_edlib;
} HbnTracebackData;

#define HbnTracebackDataDump(output_func, out, data) \
    output_func(out, "[%d, %d, %d] x [%d, %d, %d], %g\n", \
        (data)->qoff, (data)->qend, (data)->qsize, \
        (data)->soff, (data)->send, (data)->ssize, \
        (data)->ident_perc)

HbnTracebackData*
HbnTracebackDataNew();

HbnTracebackData*
HbnTracebackDataFree(HbnTracebackData* data);

BOOL
truncate_align_bad_ends(const char* qaln,
    const char* saln,
    const int aln_size,
    int* qoff,
    int* qend,
    int* soff,
    int* send,
    const char** qas_,
    const char** qae_,
    const char** sas_,
    const char** sae_);

int 
porec_compute_traceback(HbnTracebackData* data,
    int qb,
    int qe,
    int sb,
    int se,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    const int min_align_size,
    const double min_ident_perc);

int
porec_align_overhang(HbnTracebackData* data,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    const int min_align_size,
    const double min_ident_perc);

#ifdef __cplusplus
}
#endif

#endif // __HBN_TRACE_BACK_H
