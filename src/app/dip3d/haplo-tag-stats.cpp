#include "frag-hap.hpp"
#include "../../corelib/pdqsort.h"

using namespace std;

int haplo_tag_stats_main(int argc, char* argv[])
{
    if (argc <= 2) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s frag-1 [...frag-n]\n", argv[0], argv[1]);
        return 1;
    }

    FragAndContactStats stats;
    vector<FragHapInfo> fhi_list;
    for (int i = 2; i < argc; ++i) {
        fprintf(stderr, "Add %s\n", argv[i]);
        load_frag_hap_info_list(argv[i], fhi_list);
        FragHapInfo* fhia = fhi_list.data();
        size_t fhic = fhi_list.size();

        pdqsort(fhia, fhia + fhic, [](const FragHapInfo& x, const FragHapInfo& y) { 
            return (x.read_id < y.read_id) || (x.read_id == y.read_id && x.frag_id < y.frag_id); });
    
        size_t ii = 0;
        while (ii < fhic) {
            size_t jj = ii + 1;
            while (jj < fhic && fhia[ii].read_id == fhia[jj].read_id) ++jj;
            stats.add_one_list(fhia + ii, jj - ii);
            ii = jj;
        }
    }

    stats.dump_stats();

    return 0;
}