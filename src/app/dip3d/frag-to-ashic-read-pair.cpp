#include "../../corelib/fasta.hpp"
#include "../../corelib/frag_id.hpp"
#include "../../corelib/pdqsort.h"
#include "frag-hap.hpp"

#include <string>
#include <vector>

using namespace std;

static constexpr const char* kFragAlleleNameList[] = {
    "both-ref",
    "ref",
    "alt"
};

int frag_to_ashic_read_pair_main(int argc, char* argv[])
{
	fprintf(stderr, "==== argc = %d\n", argc);
    if (argc != 6) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s hap-list chr min-mapQ output\n", argv[0], argv[1]);
        exit (EXIT_FAILURE);
    }
    const char* hap_list_path = argv[2];
    const char* chr = argv[3];
    int min_mapQ = atoi(argv[4]);
    const char* output = argv[5];

    fprintf(stderr, "Frag-hap: %s\n", hap_list_path);
    fprintf(stderr, "Chromosome: %s\n", chr);
    fprintf(stderr, "min-mapQ: %d\n", min_mapQ);
    fprintf(stderr, "output: %s\n", output);
    fprintf(stderr, "\n");

    vector<FragHapInfo> frag_list;
    load_frag_hap_info_list(hap_list_path, frag_list);
    pdqsort(frag_list.begin(), frag_list.end(), [](const FragHapInfo& x, const FragHapInfo& y) {
        return x.read_id < y.read_id;
    });
    FragHapInfo* a = frag_list.data();
    size_t c = frag_list.size();

    hbn_dfopen(out, output, "w");
    size_t i = 0;
    while (i < c) {
        const char* hp1 = kFragAlleleNameList[a[i].hp];
        size_t j = i + 1;
        while (j < c && a[i].read_id == a[j].read_id) {
            const char* hp2 = kFragAlleleNameList[a[j].hp];
            if (a[i].mapQ >= min_mapQ && a[j].mapQ >= min_mapQ) fprintf(out, "%s\t%d\t%s\t%s\t%d\t%s\n", chr, a[i].subject_offset, hp1, chr, a[j].subject_offset, hp2);
            ++j;
        }
        i = j;
    }
    hbn_fclose(out);

    return 0;
}
