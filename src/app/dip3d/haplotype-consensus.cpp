#include "../../corelib/frag_id.hpp"
#include "../../corelib/pdqsort.h"
#include "../../htslib/sam.h"
#include "bam-writer.hpp"
#include "frag-hap.hpp"
#include "x-tag-bam.hpp"

using namespace std;

static void
s_load_frag_haps_from_bams(const char* bam_path, vector<FragHapInfo>& frag_list)
{
    BamChunkReader in(bam_path, 8);
    FragHapInfo frag;
    char hptagname[2] = { 'H', 'P' };
    char pitagname[2] = { 'p', 'i' };
    while (1) {
        bam1_t* bam = in.get_next_bam();
        if (!bam) break;
        in.advanve();

        const char* qn = bam_get_qname(bam);
        const int qnl = strlen(qn);
        extract_frag_id_from_name(qn, qnl, &frag.read_id, &frag.frag_id, &frag.subject_id, &frag.subject_offset);

        frag.mapQ = bam->core.qual;
        frag.frag_size = bam->core.l_qseq;

        frag.hp = 0;
        uint8_t* hptag = bam_aux_get(bam, hptagname);
        if (hptag) frag.hp = bam_aux2i(hptag);

        uint8_t* pitag = bam_aux_get(bam, pitagname);
        if (!pitag) HBN_ERR("Could not retrieve alignment identity via pi tag from BAM record");
        frag.identity = bam_aux2f(pitag);

        frag_list.push_back(frag);
    }
}

static void
s_impute_one_list_pore_c(FragHapInfo* fhia, int fhic)
{
    int hp[3]; hp[0] = hp[1] = hp[2] = 0;
    for (int i = 0; i < fhic; ++i) ++hp[fhia[i].hp];
    bool r = (hp[0] > 0) && (hp[1] > 0 || hp[2] > 0);
    if (!r) return;

    if (hp[1] > 0 && hp[2] == 0) {
        for (int i = 0; i < fhic; ++i) if (fhia[i].mapQ >= 1) fhia[i].hp = 1;
        return;
    }
    if (hp[1] == 0 && hp[2] > 0) {
        for (int i = 0; i < fhic; ++i) if (fhia[i].mapQ >= 1) fhia[i].hp = 2;
        return;
    }

    double p1 = 1.0 * hp[1] / (hp[1] + hp[2]);
    double p2 = 1.0 * hp[2] / (hp[1] + hp[2]);
    if (p1 > 0.8) {
        for (int i = 0; i < fhic; ++i) if (fhia[i].mapQ >= 1) fhia[i].hp = 1;
        return;
    }
    if (p2 > 0.8) {
        for (int i = 0; i < fhic; ++i) if (fhia[i].mapQ >= 1) fhia[i].hp = 2;
        return;
    }
}

int haplotype_consensus_main(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s whatshap-tagged-bam wrk-dir\n", argv[0], argv[1]);
        exit (EXIT_FAILURE);
    }
    const char* whatshap_tagged_bam_path = argv[2];
    const char* wrk_dir = argv[3];
    create_directory(wrk_dir);

    vector<FragHapInfo> frag_list;
    s_load_frag_haps_from_bams(whatshap_tagged_bam_path, frag_list);

    pdqsort(frag_list.begin(), frag_list.end(), [](const FragHapInfo& x, const FragHapInfo& y) { return x.read_id < y.read_id; });

    HBN_LOG("whatshap-tagged frag and contact stats for %s", whatshap_tagged_bam_path);
    frag_and_contact_stats_hap_list(frag_list.data(), frag_list.size());

    vector<FragHapInfo> cns_frag_list = frag_list;
    FragHapInfo* fhia = cns_frag_list.data();
    int fhic = cns_frag_list.size();
    int i = 0;
    while (i < fhic) {
        int j = i + 1;
        while (j < fhic && fhia[i].read_id == fhia[j].read_id) ++j;
        bool has_tag = false, has_untag = false;
        for (int k = i; k < j; ++k) {
            if (fhia[k].hp == 0) {
                has_untag = true;
            } else {
                has_tag = true;
            }
        }
        if (has_tag && has_untag) {
            s_impute_one_list_pore_c(fhia + i, j - i);
        }
        i = j;
    }

    HBN_LOG("frag and contact stats with haplotype cosensus for %s", whatshap_tagged_bam_path);
    frag_and_contact_stats_hap_list(cns_frag_list.data(), cns_frag_list.size());

    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s/snp-tagged-frag-hap-list", wrk_dir);
    save_frag_hap_info_list(frag_list.data(), frag_list.size(), path);

    snprintf(path, HBN_MAX_PATH_LEN, "%s/imputed-frag-hap-list", wrk_dir);
    save_frag_hap_info_list(cns_frag_list.data(), cns_frag_list.size(), path);

    snprintf(path, HBN_MAX_PATH_LEN, "%s/tagged.bam", wrk_dir);
    tag_bam_st(cns_frag_list.data(), cns_frag_list.size(), 8, whatshap_tagged_bam_path, path);

    return 0;
}