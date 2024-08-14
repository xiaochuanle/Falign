#include "../../corelib/lookup_table.hpp"

extern void 
kmer_cov_stats(const char* referenece_path, 
    const char* repeat_bed_path,
    int repeat_kmer_size,
    int repeat_kmer_window,
    int repeat_max_kmer_occ,
    double repeat_rep_frac,
    int non_repeat_kmer_size,
    int non_repeat_kmer_window,
    int non_repeat_max_kmer_occ,
    double non_repeat_rep_frac,
    int num_threads);

int main(int argc, char* argv[])
{
    const char* reference_path = argv[1];
    const char* repeat_bed_path = argv[2];

    kmer_cov_stats(reference_path, repeat_bed_path, 15, 2, 100, 0.1, 15, 2, 100, 0.0001, 48);
    return 0;
}
