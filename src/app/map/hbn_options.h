#ifndef __HBN_OPTIONS_H
#define __HBN_OPTIONS_H

#include <stdlib.h>
#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef enum {
    eOutputFmt_SAM = 0,
    eOutputFmt_PAF,
    eOutputFmt_Invalid
} EOutputFmt;

const char* output_format_name(EOutputFmt fmt);

EOutputFmt name_to_output_format(const char* name);

typedef struct {
    /// input sequence options
    size_t          query_batch_size;
    size_t          query_upto;

    /// ddf scoring options on detecting candidate subject subsequences
    int             kmer_size;
    int             kmer_window;
    double          rep_frac;
    int             max_kmer_occ;
    double          ddf;
    int             kmer_dist;
    int             chain_score;

    /// restrict search or results
    double          perc_identity;
    int             max_hsps_per_subject;
    int             hitlist_size;

    int             num_threads;
    EOutputFmt      outfmt;
    int             dump_by_file;

    const char*     enzyme;
    const char*     output;
    const char*     output_dir;
} HbnProgramOptions;

void dump_usage_simple(const char* pn);

void dump_usage_full(const char* pn);

int parse_arguments(int argc, char* argv[], HbnProgramOptions* opts);

#ifdef __cplusplus
}
#endif

#endif // __HBN_OPTIONS_H