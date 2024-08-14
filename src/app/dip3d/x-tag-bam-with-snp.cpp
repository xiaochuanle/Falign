#include "x-tag-bam-with-snp.hpp"

#include "../../corelib/frag_id.hpp"
#include "../../corelib/pdqsort.h"
#include "sam_map_info.hpp"

#include <mutex>
#include <set>

using namespace std;

class SnpTagBamThreadWorkData
{
public:
    SnpTagBamThreadWorkData(BamTagOptions* options,
        HbnUnpackedDatabase* reference,
        SeqName2IdMap* ref_name2id,
        VcfPhasedHetSnpList* snp_list,
        BamChunkReader* bam,
        vector<FragHapInfo>* frag_hap_list)
    :   
        M_options(options),
        M_reference(reference),
        M_ref_name2id(ref_name2id),
        M_bam(bam),
        M_snps(snp_list),
        M_frag_hap_list(frag_hap_list),
        M_hap_amb_frag(0)
    {

    }

    ~SnpTagBamThreadWorkData() {

    }



    HbnUnpackedDatabase* reference() {
        return M_reference;
    }

    SeqName2IdMap* ref_name2id() {
        return M_ref_name2id;
    }

    bam1_t* get_next_sam() {
        lock_guard<mutex> _(M_bam_mutex);
        bam1_t* bam = M_bam->get_next_bam();
        if (bam) M_bam->advanve();
        return bam;
    }

    VcfPhasedHetSnpList* snps() {
        return M_snps;
    }

    void add_one_frag_hap(FragHapInfo& hap) {
        lock_guard<mutex> _(M_frag_hap_list_mutex);
        M_frag_hap_list->push_back(hap);
    }

    sam_hdr_t* sam_hdr() {
        return M_bam->sam_hdr();
    }

    BamTagOptions* options() {
        return M_options;
    }

    void add_one_hap_amb_frag() {
        lock_guard<mutex> _(M_frag_hap_list_mutex);
        ++M_hap_amb_frag;
    }

    size_t hap_amb_frag() {
        return M_hap_amb_frag;
    }

private:
    BamTagOptions*          M_options;
    HbnUnpackedDatabase*    M_reference;
    SeqName2IdMap*          M_ref_name2id;
    BamChunkReader*         M_bam;
    mutex                   M_bam_mutex;
    VcfPhasedHetSnpList*    M_snps;
    vector<FragHapInfo>*    M_frag_hap_list;
    mutex                   M_frag_hap_list_mutex;
    size_t                  M_hap_amb_frag;
};

static void
s_resolve_haplotype_at_one_overlapped_snp(const char* qas, const char* sas, const int as_size,
    const int* qas_pos_list, const int* sas_pos_list,
    int qb, int qe, int sb, int se,
    const int match_base,
    const PhasedHetSNP& snp,
    set<int>& snp_sites,
    int* hp_cnt)
{
    int as_i = 0;
    while (as_i < as_size) {
        if (sas[as_i] != GAP_CHAR && sas_pos_list[as_i] == snp.soff) break;
        ++as_i;
    }
    hbn_assert(as_i < as_size);
    hbn_assert(sas[as_i] != GAP_CHAR);
    hbn_assert(sas_pos_list[as_i] == snp.soff);
    hbn_assert(sas[as_i] == snp.REF);
    if (qas[as_i] == GAP_CHAR) { ++hp_cnt[0]; return; }

    bool l_perfect = true;
    for (int i = 1; i <= match_base && as_i - i >= 0; ++i) {
        if (sas[as_i - i] != GAP_CHAR) {
            int soff = sas_pos_list[as_i - i];
            hbn_assert(soff >= 0);
            if (snp_sites.find(soff) != snp_sites.end()) continue;
        }
        if (qas[as_i - i] != sas[as_i - i]) l_perfect = false;
    }

    bool r_perfect = true;
    for (int i = 1; i <= match_base && as_i + i < as_size; ++i) {
        if (sas[as_i + i] != GAP_CHAR) {
            int soff = sas_pos_list[as_i + i];
            if (snp_sites.find(soff) != snp_sites.end()) continue;
        }
        if (qas[as_i + i] != sas[as_i + i]) r_perfect = false;
    }

#if 0
    if (l_perfect == false || r_perfect == false) {
        fprintf(stderr, "=============\n");
        for (int i = 1; i <= 3 && as_i - i >= 0; ++i) fprintf(stderr, "%d %c - %c\n", as_i - i, qas[as_i - i], sas[as_i - i]);
        fprintf(stderr, "*%d %c - %c\n", as_i,  qas[as_i], sas[as_i]);
        for (int i = 1; i <= 3 && as_i + i < as_size; ++i) fprintf(stderr, "%d %c - %c\n", as_i + i, qas[as_i + i], sas[as_i + i]);
        fprintf(stderr, "%d\t%c\t%c\t%d\t%d\n", snp.soff, snp.REF, snp.ALT, snp.gt[0], snp.gt[1]);
    }
#endif
    if (l_perfect == false || r_perfect == false) { ++hp_cnt[0]; return; }

    int hp = 0;
    if (qas[as_i] == snp.REF) {
        hp = snp.gt[0] + 1;
    } else if (qas[as_i] == snp.ALT) {
        hp = snp.gt[1] + 1;
    }
    hbn_assert(hp == 0 || hp == 1 || hp == 2);
    ++hp_cnt[hp];
}

static bool
snp_tag_one_bam(BamTagOptions* options,
        HbnUnpackedDatabase* reference,
        SeqName2IdMap* ref_name2id,
        VcfPhasedHetSnpList* snp_list,
        sam_hdr_t* sam_hdr,
        bam1_t* bam,
        SAM_MapInfo* sam_maps,
        FragHapInfo* hap)
{
    if (!sam_maps->parse(sam_hdr, bam, reference, ref_name2id)) return false;
    //sam_maps->dump();
    const char* qn = sam_maps->query_name();
    const int qnl = strlen(qn);
    extract_frag_id_from_name(qn, qnl, &hap->read_id, &hap->frag_id, &hap->subject_id, &hap->subject_offset);
    hap->frag_size = sam_maps->query_size();
    hap->identity = sam_maps->identity();
    hap->mapQ = sam_maps->mapQ();
    hap->hp = 0;
    hap->amb_hp = 0;
    bool r = (sam_maps->mapQ() >= options->snp_tag_min_mapQ)
             &&
             (sam_maps->identity() >= options->snp_tag_identity);
    if (!r) return true;

    int qb = sam_maps->qb();
    int qe = sam_maps->qe();
    int sb = sam_maps->sb();
    int se = sam_maps->se();
    const char* qas = sam_maps->qas();
    const char* sas = sam_maps->sas();
    const int as_size = sam_maps->as_size();
    const int* qas_pos_list = sam_maps->qas_pos_list();
    const int* sas_pos_list = sam_maps->sas_pos_list();

    PhasedHetSNP* chr_snp_list = nullptr;
    int chr_snp_list_idx = -1;
    int chr_snp_list_size = 0;
    if (!snp_list->get_snp_idx(sam_maps->sname(), sb, &chr_snp_list, &chr_snp_list_idx, &chr_snp_list_size)) {
        return true;
    }
    if (chr_snp_list_idx) hbn_assert(sb > chr_snp_list[chr_snp_list_idx-1].soff);
    hbn_assert(sb <= chr_snp_list[chr_snp_list_idx].soff);
    set<int> ovlp_snps;
    for (int i = chr_snp_list_idx; i < chr_snp_list_size; ++i) {
        if (chr_snp_list[i].soff >= se) break;
        ovlp_snps.insert(chr_snp_list[i].soff);
    }

    int hp_cnt[3] = { 0, 0, 0 };
    for (int i = chr_snp_list_idx; i < chr_snp_list_size; ++i) {
        if (chr_snp_list[i].soff >= se) break;
        s_resolve_haplotype_at_one_overlapped_snp(qas, sas, as_size, qas_pos_list, sas_pos_list,
            qb, qe, sb, se, options->snp_tag_match_base, 
            chr_snp_list[i], ovlp_snps, hp_cnt);
    }

    //fprintf(stderr, "haps: %d, %d, %d\n", hp_cnt[0], hp_cnt[1], hp_cnt[2]);
    if (hp_cnt[1] > hp_cnt[2]) {
        hap->hp = 1;
    } else if (hp_cnt[1] < hp_cnt[2]) {
        hap->hp = 2;
    } 
    return true;
}

////////////////

static void
s_resolve_haplotype_at_one_overlapped_snp_amb(const char* qas, const char* sas, const int as_size,
    const int* qas_pos_list, const int* sas_pos_list,
    int qb, int qe, int sb, int se,
    const int match_base,
    const PhasedHetSNP& snp,
    set<int>& snp_sites,
    int* hp_cnt)
{
    int as_i = 0;
    while (as_i < as_size) {
        if (sas[as_i] != GAP_CHAR && sas_pos_list[as_i] == snp.soff) break;
        ++as_i;
    }
    hbn_assert(as_i < as_size);
    hbn_assert(sas[as_i] != GAP_CHAR);
    hbn_assert(sas_pos_list[as_i] == snp.soff);
    hbn_assert(sas[as_i] == snp.REF);
    if (qas[as_i] == GAP_CHAR) { ++hp_cnt[0]; return; }

    int hp = 0;
    if (qas[as_i] == snp.REF) {
        hp = snp.gt[0] + 1;
    } else if (qas[as_i] == snp.ALT) {
        hp = snp.gt[1] + 1;
    }
    hbn_assert(hp == 0 || hp == 1 || hp == 2);
    ++hp_cnt[hp];
}

static bool
snp_tag_one_bam_amb(BamTagOptions* options,
        HbnUnpackedDatabase* reference,
        SeqName2IdMap* ref_name2id,
        VcfPhasedHetSnpList* snp_list,
        sam_hdr_t* sam_hdr,
        bam1_t* bam,
        SAM_MapInfo* sam_maps,
        FragHapInfo* hap)
{
    if (!sam_maps->parse(sam_hdr, bam, reference, ref_name2id)) return false;
    //sam_maps->dump();
    const char* qn = sam_maps->query_name();
    const int qnl = strlen(qn);
    extract_frag_id_from_name(qn, qnl, &hap->read_id, &hap->frag_id, &hap->subject_id, &hap->subject_offset);
    hap->frag_size = sam_maps->query_size();
    hap->identity = sam_maps->identity();
    hap->mapQ = sam_maps->mapQ();
    hap->hp = 0;
    hap->amb_hp = 0;
    //bool r = (sam_maps->mapQ() >= 5) && (sam_maps->identity() >= 85.0);
    bool r = (sam_maps->mapQ() >= 5);
    if (!r) return true;

    int qb = sam_maps->qb();
    int qe = sam_maps->qe();
    int sb = sam_maps->sb();
    int se = sam_maps->se();
    const char* qas = sam_maps->qas();
    const char* sas = sam_maps->sas();
    const int as_size = sam_maps->as_size();
    const int* qas_pos_list = sam_maps->qas_pos_list();
    const int* sas_pos_list = sam_maps->sas_pos_list();

    PhasedHetSNP* chr_snp_list = nullptr;
    int chr_snp_list_idx = -1;
    int chr_snp_list_size = 0;
    if (!snp_list->get_snp_idx(sam_maps->sname(), sb, &chr_snp_list, &chr_snp_list_idx, &chr_snp_list_size)) {
        return true;
    }
    if (chr_snp_list_idx) hbn_assert(sb > chr_snp_list[chr_snp_list_idx-1].soff);
    hbn_assert(sb <= chr_snp_list[chr_snp_list_idx].soff);
    set<int> ovlp_snps;
    for (int i = chr_snp_list_idx; i < chr_snp_list_size; ++i) {
        if (chr_snp_list[i].soff >= se) break;
        ovlp_snps.insert(chr_snp_list[i].soff);
    }

    int hp_cnt[3] = { 0, 0, 0 };
    for (int i = chr_snp_list_idx; i < chr_snp_list_size; ++i) {
        if (chr_snp_list[i].soff >= se) break;
        s_resolve_haplotype_at_one_overlapped_snp_amb(qas, sas, as_size, qas_pos_list, sas_pos_list,
            qb, qe, sb, se, options->snp_tag_match_base, 
            chr_snp_list[i], ovlp_snps, hp_cnt);
    }

    //fprintf(stderr, "haps: %d, %d, %d\n", hp_cnt[0], hp_cnt[1], hp_cnt[2]);
    if (hp_cnt[1] > hp_cnt[2]) {
        hap->amb_hp = 1;
    } else if (hp_cnt[1] < hp_cnt[2]) {
        hap->amb_hp = 2;
    }
    return true;
}

//////////////////

static void
s_resolve_haplotypes_one_read_amb(BamTagOptions* options, FragHapInfo* a, int c)
{
    int read_hp = 0;
    int read_hp1 = 0, read_hp2 = 0;
    for (int i = 0; i < c; ++i) {
        if (a[i].hp == 1) ++read_hp1;
        if (a[i].hp == 2) ++read_hp2;
    }

#if 0
   if (read_hp1 > 0 && read_hp2 > 0) {
	double hp1f = 1.0 * read_hp1 / (read_hp1 + read_hp2);
	double hp2f = 1.0 * read_hp2 / (read_hp1 + read_hp2);
	if (hp1f >= 0.8) {
		for (int i = 0; i < c; ++i) if (a[i].hp == 2) a[i].hp = 1;
	}
	if (hp2f >= 0.8) {
		for (int i = 0; i < c; ++i) if (a[i].hp == 1) a[i].hp = 2;
	}
   }

    read_hp = 0;
    read_hp1 = 0, read_hp2 = 0;
    for (int i = 0; i < c; ++i) {
        if (a[i].hp == 1) ++read_hp1;
        if (a[i].hp == 2) ++read_hp2;
    }
#endif

    if (read_hp1 == 0 && read_hp2 > 0) read_hp = 2;
    if (read_hp1 > 0 && read_hp2 == 0) read_hp = 1;

    int read_amb_hp1 = 0, read_amb_hp2 = 0, read_amb_hp = 0;
    for (int i = 0; i < c; ++i) {
        if (a[i].hp > 0) continue;
        if (a[i].mapQ < 5) continue;
        //if (a[i].identity < options->snp_tag_identity) continue;
        if (a[i].amb_hp == 1) ++read_amb_hp1;
        if (a[i].amb_hp == 2) ++read_amb_hp2;
    }

    if (read_hp == 0) {
        if (read_amb_hp1 == 0 && read_amb_hp2 > 0) read_amb_hp = 2;
        if (read_amb_hp1 > 0 && read_amb_hp2 == 0) read_amb_hp = 1;
        if (read_amb_hp == 0) return;

        for (int i = 0; i < c; ++i) {
            if (a[i].hp > 0) continue;
            if (a[i].mapQ < 5) continue;
            //if (a[i].identity < options->snp_tag_identity) continue;
            if (a[i].amb_hp == 0) continue;
            a[i].hp = a[i].amb_hp;
        }

        return;
    }

    for (int i = 0; i < c; ++i) {
        if (a[i].hp > 0) continue;
        if (a[i].mapQ < 5) continue;
        //if (a[i].identity < options->snp_tag_identity) continue;
        if (a[i].amb_hp == read_hp) a[i].hp = read_hp;
    }
}

static void
s_resolve_haplotypes_amb(BamTagOptions* options, FragHapInfo* a, int c)
{
    pdqsort(a, a + c, [](const FragHapInfo& x, const FragHapInfo& y) { return x.read_id < y.read_id; });
    int i = 0;
    while (i < c) {
        int j = i + 1;
        while (j < c && a[i].read_id == a[j].read_id) ++j;
        s_resolve_haplotypes_one_read_amb(options, a + i, j - i);
        i = j;
    }
}

/////////////////

static void*
snp_tag_thread(void* params)
{
    SnpTagBamThreadWorkData* data = (SnpTagBamThreadWorkData*)(params);
    BamTagOptions* options = data->options();
    HbnUnpackedDatabase* reference = data->reference();
    SeqName2IdMap* ref_name2id = data->ref_name2id();
    sam_hdr_t* sam_hdr = data->sam_hdr();
    VcfPhasedHetSnpList* snp_list = data->snps();
    SAM_MapInfo sam_maps;

    int cnt = 0;
    FragHapInfo hap, hap_amb;
    while (1) {
        bam1_t* bam = data->get_next_sam();
        if (!bam) break;

        if (snp_tag_one_bam(options, reference, ref_name2id, snp_list, sam_hdr, bam, &sam_maps, &hap)) {
            snp_tag_one_bam_amb(options, reference, ref_name2id, snp_list, sam_hdr, bam, &sam_maps, &hap_amb);
            hap.amb_hp = hap_amb.amb_hp;
            data->add_one_frag_hap(hap);
        }

        bam_destroy1(bam);
        //if (++cnt == 100) break;
    }

    return nullptr;
}

void snp_tag_bam_mt(BamTagOptions* options, std::vector<FragHapInfo>& hap_list)
{
    HbnUnpackedDatabase referenece(options->reference_path);
    referenece.load_next_batch();
    SeqName2IdMap ref_name2id;
    int num_chr = referenece.NumSeqs();
    for (int i = 0; i < num_chr; ++i) ref_name2id.add_one_name(referenece.SeqName(i));
    ref_name2id.build_name2id_map();
    VcfPhasedHetSnpList snps(options->vcf_path);
    BamChunkReader bam(options->bam_path, 8);
    pthread_t jobs[options->num_threads];
    SnpTagBamThreadWorkData data(options,
        &referenece,
        &ref_name2id,
        &snps,
        &bam,
        &hap_list);

    for (int i = 0; i < options->num_threads; ++i) {
        pthread_create(jobs + i, nullptr, snp_tag_thread, &data);
    }
    for (int i = 0; i < options->num_threads; ++i) {
        pthread_join(jobs[i], nullptr);
    }

    s_resolve_haplotypes_amb(options, hap_list.data(), hap_list.size());
}
