#ifndef __TRIM_OVERLAP_SUBSEQ_H
#define __TRIM_OVERLAP_SUBSEQ_H

#include "hbn_options.h"
#include "chain_align_list.hpp"
#include "../../algo/hbn_traceback.h"
#include "../../algo/seq_loader.h"
#include "../../corelib/restrict_enzyme_loci_list.h"

void
trim_overlap_subseqs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    SeqReader* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const char* query_qv,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    PoreCAlign* pca_a,
    int pca_c,
    const EChainType chain_type,
    int extended_subject,
    EOutputFmt outfmt,
    kstring_t* out);

#endif // __TRIM_OVERLAP_SUBSEQ_H