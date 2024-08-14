#include "../../corelib/pdqsort.h"
#include "vcf-reader.hpp"

using namespace std;

static int
get_most_variant_phase_set(pair<PhasedHetSNP, bcf1_t*>* a, size_t c)
{
    pdqsort(a, a + c, [](const pair<PhasedHetSNP, bcf1_t*>& x, const pair<PhasedHetSNP, bcf1_t*>& y) {
        return x.first.phase_set < y.first.phase_set;
    });

    size_t max_cnt = 0;
    int max_phase_set = -1;
    size_t i = 0;
    while (i < c) {
        size_t j = i + 1;
        while (j < c && a[i].first.phase_set == a[j].first.phase_set) ++j;
        
        size_t cnt = j - i;
        if (cnt > max_cnt) {
            max_cnt = cnt;
            max_phase_set = a[i].first.phase_set;
        }

        i = j;
    }

    pdqsort(a, a + c, [](const pair<PhasedHetSNP, bcf1_t*>& x, const pair<PhasedHetSNP, bcf1_t*>& y) {
        return x.first.soff < y.first.soff;
    });

    return max_phase_set;
}

int extract_mvp_het_snp_main(int argc, char* argv[])
{
    const char* input = argv[2];
    const char* output = argv[3];
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s input-vcf output-vcf\n", argv[0], argv[1]);
        exit (1);
    }

    map<string, int>* phase_set_name2id = nullptr;
    vector<string>* phase_set_name_list = nullptr;
    vector<pair<PhasedHetSNP, bcf1_t*>> phased_snps;
    bcf_hdr_t* bcf_hdr = nullptr;
    load_phased_het_snps_from_vcf_file(input, true, &phase_set_name2id, &phase_set_name_list, &bcf_hdr, phased_snps);
    
    vcfFile* out = vcf_open(output, "w");
    hts_set_threads(out, 8);
    if (vcf_hdr_write(out, bcf_hdr)) HBN_ERR("FAIL at writing VCF header to %s", output);

    pair<PhasedHetSNP, bcf1_t*>* ppsa = phased_snps.data();
    size_t ppsc = phased_snps.size();
    size_t i = 0;
    while (i < ppsc) {
        size_t j = i + 1;
        while (j < ppsc && ppsa[i].first.sid == ppsa[j].first.sid) ++j;

        int max_phase_set = get_most_variant_phase_set(ppsa + i, j - i);
        for (size_t k = i; k < j; ++k) {
            if (ppsa[k].first.phase_set == max_phase_set)
            if (vcf_write(out, bcf_hdr, ppsa[k].second)) HBN_ERR("FAIL at writing VCF record");
        }

        i = j;
    }
    vcf_close(out);

    fprintf(stderr, "free name2id\n");
    if (phase_set_name2id) delete phase_set_name2id;
    fprintf(stderr, "free names\n");
    if (phase_set_name_list) delete phase_set_name_list;
    fprintf(stderr, "free vcf hdr\n");
    if (bcf_hdr) bcf_hdr_destroy(bcf_hdr);

    fprintf(stderr, "Free bcf records\n");
    for (auto& x : phased_snps) if (x.second) bcf_destroy1(x.second);

    return 0;
}