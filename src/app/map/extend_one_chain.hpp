#ifndef __EXTEND_ONE_CHAIN_HPP
#define __EXTEND_ONE_CHAIN_HPP

#include "../../algo/hbn_traceback.h"
#include "../../algo/hbn_traceback_aux.h"
#include "chain_align_list.hpp"

#include <string>
#include <vector>

BOOL 
refine_pi_for_one_pca(HbnTracebackData* tbck_data,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_size,
    const u8* subject,
    PoreCAlign* pca);

double
compute_pi_for_one_align(const char* qas1,
    const char* sas1,
    const int as_size1,
    const char* qas2 = nullptr,
    const char* sas2 = nullptr,
    const int as_size2 = 0);

void
extend_one_chain(const int query_id,
    const int query_dir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    DDFS_Seed* sa,
    int sc,
    const double perc_identity,
    HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    int* cov_stats,
    std::vector<PoreCAlign>& align_list);

#endif // __EXTEND_ONE_CHAIN_HPP