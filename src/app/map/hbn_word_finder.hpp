#ifndef __HBN_WORD_FINDER_HPP
#define __HBN_WORD_FINDER_HPP

#include "../../corelib/lookup_table.hpp"
#include "hbn_options.hpp"

#include <random>
#include <vector>

class HbnWordFinder
{
public:
    HbnWordFinder(HbnUnpackedDatabase& db,
        const char* repeat_bed_path,
        const bool skip_repeat_regions,
        int num_threads,
        int kmer_size,
        int kmer_window,
        int max_kmer_occ,
        double repeat_rep_frac,
        double non_repeat_rep_frac);
    
    ~HbnWordFinder()
    {
        delete M_lktbl;
    }

    void extract_kmer_matches(const int query_id,
        const u8* fwd_query,
        const u8* rev_query,
        const int query_size,
        std::vector<KmerMatch>& km_list,
        std::mt19937* gen,
        std::uniform_int_distribution<>* dist);
    
private:

    int                                 M_kmer_size;
    HbnLookupTable*                     M_lktbl;
};

#endif // __HBN_WORD_FINDER_HPP