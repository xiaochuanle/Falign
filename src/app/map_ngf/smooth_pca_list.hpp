#ifndef __SMOOTH_PCA_LIST_HPP
#define __SMOOTH_PCA_LIST_HPP

#include "../../algo/hbn_traceback.h"
#include "../../algo/hbn_traceback_aux.h"
#include "chain_align_list.hpp"
#include "hbn_options.h"

#include <string>
#include <vector>

void
smooth_pca_list(std::vector<PoreCAlign>& pca_list,
    EChainType& chain_type,
    std::vector<PoreCAlign>& all_pca_list,
    SeqReader* ref,
    const char* read_name,
    const u8* fwd_read,
    const u8* rev_read,
    const int* vdfa, 
    const int vdfc,
    const int enzyme_size,
    HbnTracebackData* tbck_data);

#endif // __SMOOTH_PCA_LIST_HPP