#ifndef __FIX_PCA_OFFSET_WITH_ENZYME_INFO_H
#define __FIX_PCA_OFFSET_WITH_ENZYME_INFO_H

#include "../../corelib/restrict_enzyme_loci_list.h"
#include "chain_align_list.hpp"

int 
fix_pca_with_enzyme_info(RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    int qdir,
    int qsize,
    int sid,
    int ssize,
    int* qb_,
    int* qe_,
    int* sb_,
    int* se_,
    int* enzyme_qb_,
    int* enzyme_qe_,
    int* enzyme_sb_,
    int* enzyme_se_,
    int* l_enzyme_match_,
    int* r_enzyme_match_);

void
fix_imperfect_ends(HbnTracebackData* tbck_data,
    const u8* fwd_query,
    const u8* rev_query,
    const u8* subject,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    PoreCAlign* pca_a,
    int pca_c,
    const double min_perc_identity,
    std::vector<PoreCAlign>& new_pca_list);

#endif