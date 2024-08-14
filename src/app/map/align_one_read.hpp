#ifndef __ALIGN_ONE_READ_HPP
#define __ALIGN_ONE_READ_HPP

#include "map_one_volume.hpp"
#include "trim_overlap_subseq.hpp"

#include <vector>

void
align_one_read(MapThreadData* data, 
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    std::vector<PoreCAlign>& all_pca_list,
    TrimPcaList& trim_pca_list);

#endif // __ALIGN_ONE_READ_HPP