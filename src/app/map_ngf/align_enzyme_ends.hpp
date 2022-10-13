#ifndef __ALIGN_ENZYME_ENDS_HPP
#define __ALIGN_ENZYME_ENDS_HPP

#include "../../corelib/restrict_enzyme_loci_list.h"
#include "chain_align_list.hpp"

#include <vector>

int 
align_enzyme_ends(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const int qid,
    const u8* fwd_query,
    const u8* rev_query,
    int qdir,
    int qsize,
    int sid,
    const u8* fwd_subject,
    int ssize,
    int qb,
    int qe,
    int sb,
    int se,
    double perc_identity,
    int* cov_stats,
    std::vector<PoreCAlign>& pca_list);

#endif // __ALIGN_ENZYME_ENDS_HPP