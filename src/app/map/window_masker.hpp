#ifndef __WINDOW_MASKER_HPP
#define __WINDOW_MASKER_HPP

#include "../../corelib/lookup_table.hpp"
#include "../../corelib/seq_name2id_map.hpp"
#include "../../corelib/unpacked_seqdb.hpp"

#include <vector>

void compute_repeat_regions(HbnUnpackedDatabase* reference, 
    const int min_masked_intv_size,
    std::vector<ChrIntv>& repeat_intv_list);

void load_bed_list(const char* bed_path, 
    SeqName2IdMap& refname2id, 
    std::vector<ChrIntv>& intv_list);

void compute_complement_intv_list(ChrIntv* a, size_t c, 
    HbnUnpackedDatabase& reference, 
    std::vector<ChrIntv>& complement_intv_list);

void repeat_intv_stats(ChrIntv* a, size_t c, const HbnUnpackedDatabase* reference);

#endif // __WINDOW_MASKER_HPP