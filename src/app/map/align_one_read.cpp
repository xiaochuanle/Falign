#include "align_one_read.hpp"

#include <algorithm>
#include <cmath>
#include <vector>

#include "../../algo/make_candidate_kmer_chain.h"
#include "trim_overlap_subseq.hpp"
#include "extend_hit_list.hpp"
#include "smooth_pca_list.hpp"

using namespace std;

static const double kHomProb = 0.75;

typedef struct {
    int sid;
    int max_score;
    int hit_offset;
    int hit_cnt;
} SubjectInitHitInfo;

////////////////// aux function

static double 
s_map_q_p1(EChainType chain_type)
{
    if (chain_type == ePerfectChain) {
        return 0.9999999;
    } else if (chain_type == eCompleteChain) {
        return 0.99999;
    } else{
        return 0.999;
    }
}

static void
s_find_hom_pca(const EChainType chain_type, PoreCAlign* aln, PoreCAlign* all_pca_a, int all_pca_c, vector<PoreCAlign>& hom_pca_list)
{
    const double chain_type_prob = s_map_q_p1(chain_type);
    PoreCAlign* pi = aln;
    for (int j = 0; j < all_pca_c; ++j) {
        PoreCAlign* pj = all_pca_a + j;
	    if (pi->qdir == pj->qdir && pi->qoff == pj->qoff && pi->qend == pj->qend && pi->sid == pj->sid && pi->soff == pj->soff && pi->send == pj->send) {
			hom_pca_list.push_back(*pj);
            continue;
	    }

        if (pi->chain_qoff == pj->chain_qoff && pi->chain_qend == pj->chain_qend) {
            if (pi->sid == pj->sid && pi->soff <= pj->soff && pi->send > pj->soff) {
                continue;
            }
            if (pi->sid == pj->sid && pj->soff <= pi->soff && pj->send > pi->soff) {
                continue;
            }

            hbn_assert(pi->pi > 0.0);
            double p2 = fabs(1.0 - pj->pi / pi->pi);
            int map_q =ceil(-10.0 * log(1.0 - chain_type_prob * p2));
            if (map_q < 5) hom_pca_list.push_back(*pj);
        }
    }
}

static void
s_pick_hom_pca_kmer(const EChainType chain_type, PoreCAlign* hom_pca_a, int hom_pca_c, PoreCAlign* pca_a, int pca_c, PoreCAlign* result)
{
    result->map_q = 0;
    const double chain_type_prob = s_map_q_p1(chain_type);
    sort(hom_pca_a, hom_pca_a + hom_pca_c, [](const PoreCAlign& x, const PoreCAlign& y) { return x.chain_score > y.chain_score; });
    int s0 = hom_pca_a[0].chain_score;
    int s1 = hom_pca_a[1].chain_score;
    if (s0 - s1 > 1) {
        double p2 = fabs(1.0 - 1.0 * s1 / s0);
        int map_q = ceil(-10.0 * log(1.0 - chain_type_prob * kHomProb * p2));
#if 0
            fprintf(stderr, "0 ============= select hom pca %d by chain score\n", supp_list[0].id);
            for (int x = 0; x < hom_pca_c; ++x) dump_chain_pca(fprintf, stderr, hom_pca_a[x], x);
            fprintf(stderr, "s0 = %d, s1 = %d, p1 = %g, p2 = %g, p3 = %g, map_q %d\n", s0, s1, chain_type_prob, kHomProb, p2, map_q);
#endif

        if (map_q > 60) map_q = 60;
        if (map_q < 0) map_q = 0;

        *result = hom_pca_a[0];
        result->map_q = map_q;     
    }
}

static void
s_pick_hom_pca_chr(const EChainType chain_type, PoreCAlign* hom_pca_a, int hom_pca_c, PoreCAlign* pca_a, int pca_c, PoreCAlign* result)
{
    result->map_q = 0;
    const double chain_type_prob = s_map_q_p1(chain_type);
    struct PcaSupportInfo { int id; int supp; };
    vector<PcaSupportInfo> supp_list(hom_pca_c);
    for (size_t i = 0; i < hom_pca_c; ++i) { supp_list[i].id = i; supp_list[i].supp = 0; }

    for (int i = 0; i < hom_pca_c; ++i) {
        PoreCAlign* pi = hom_pca_a + i;
        for (int j = 0; j < pca_c; ++j) {
            PoreCAlign* pj = pca_a + j;
            if (pi->chain_qoff == pj->chain_qoff && pi->chain_qend == pj->chain_qend) continue;
            if (pj->map_q < 5) continue;
            if (pi->sid == pj->sid) ++supp_list[i].supp;
        }
    }
    sort(supp_list.begin(), supp_list.end(), [](const PcaSupportInfo& x, const PcaSupportInfo& y) { return x.supp > y.supp; });
    int cnt0 = supp_list[0].supp;
    int cnt1 = supp_list[1].supp;
    if (cnt0 - cnt1 > 1) {
        double p2 = fabs(1.0 - 1.0 * cnt1 / cnt0);
        int map_q = ceil(-10.0 * log(1.0 - chain_type_prob * kHomProb * p2));
#if 0
        fprintf(stderr, "1 ============= select hom pca %d by chr support\n", supp_list[0].id);
        for (int x = 0; x < hom_pca_c; ++x) dump_chain_pca(fprintf, stderr, hom_pca_a[x], x);
        fprintf(stderr, "cnt0 = %d, cnt1 = %d, p1 = %g, p2 = %g, p3 = %g, map_q %d\n", cnt0, cnt1, chain_type_prob, kHomProb, p2, map_q);
#endif

        if (map_q > 60) map_q = 60;
        if (map_q < 0) map_q = 0;
        *result = hom_pca_a[supp_list[0].id];
        result->map_q = map_q;
    }
}

static bool 
s_pick_hom_list(const EChainType chain_type, PoreCAlign* hom_pca_a, int hom_pca_c, PoreCAlign* pca_a, int pca_c, PoreCAlign* result)
{
    PoreCAlign p;
    int max_map_q = 0;

    s_pick_hom_pca_kmer(chain_type, hom_pca_a, hom_pca_c, pca_a, pca_c, &p);
    if (p.map_q > max_map_q) {
        *result = p;
        max_map_q = p.map_q;
    }

    s_pick_hom_pca_chr(chain_type, hom_pca_a, hom_pca_c, pca_a, pca_c, &p);
    if (p.map_q > max_map_q) {
        *result = p;
        max_map_q = p.map_q;
    }

    return max_map_q > 0;
}

static void
s_set_map_q_for_hom_pca(PoreCAlign* pca_a, int pca_c, const EChainType chain_type, PoreCAlign* all_pca_a, int all_pca_c)
{
    vector<PoreCAlign> hom_pca_list;
    PoreCAlign pca;
    for (int i = 0; i < pca_c; ++i) {
        if (pca_a[i].map_q >= 5) continue;
        hom_pca_list.clear();
        s_find_hom_pca(chain_type, pca_a + i, all_pca_a, all_pca_c, hom_pca_list);

        //fprintf(stderr, "---- find %d hom pca\n", hom_pca_list.size()); dump_chain_pca(fprintf, stderr, pca_a[i], i);
        for (int x = 0; x < hom_pca_list.size(); ++x) {
            PoreCAlign& hpca = hom_pca_list[x];
            //fprintf(stderr, "\t"); dump_chain_pca(fprintf, stderr, hpca, x); fprintf(stderr, "chain score = %d\n", hpca.chain_score);
        }

        if (hom_pca_list.size() < 2) continue;
        bool x = s_pick_hom_list(chain_type, hom_pca_list.data(), hom_pca_list.size(), pca_a, pca_c, &pca);
        if (!x) continue;
        pca.map_q = -pca.map_q;
        pca_a[i] = pca;
    }
    for (int i = 0; i < pca_c; ++i) pca_a[i].map_q = abs(pca_a[i].map_q);
}

static void
s_set_hom_pca_idx(PoreCAlign* pca_a, int pca_c, vector<PoreCAlign>& all_pca_list)
{
    for (int i = 0; i < pca_c; ++i) pca_a[i].qc = -1;
    PoreCAlign* apa = all_pca_list.data();
    int apc = all_pca_list.size();
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pi = pca_a + i;
        double rep_pi = 0.0;
	    int rep_j = -1;
        for (int j = 0; j < apc; ++j) {
            PoreCAlign* pj = apa + j;
	        if (pi->qdir == pj->qdir && pi->qoff == pj->qoff && pi->qend == pj->qend && pi->sid == pj->sid && pi->soff == pj->soff && pi->send == pj->send) {
			    continue;
	        }
            if (pi->chain_qoff == pj->chain_qoff && pi->chain_qend == pj->chain_qend) {
                if (pi->sid == pj->sid && pi->soff <= pj->soff && pi->send > pj->soff) {
                    continue;
                }
                if (pi->sid == pj->sid && pj->soff <= pi->soff && pj->send > pi->soff) {
                    continue;
                }
                if (pj->pi > rep_pi) { rep_pi = pj->pi; rep_j = j; }
            }
        }
	    pi->qc = rep_j;
    }    
}

static void
s_set_map_q(EChainType chain_type, PoreCAlign* pca_a, int pca_c, vector<PoreCAlign>& all_pca_list)
{
    PoreCAlign* apa = all_pca_list.data();
    int apc = all_pca_list.size();
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pi = pca_a + i;
        double rep_pi = 0.0;
	    int rep_j = -1;
        for (int j = 0; j < apc; ++j) {
            PoreCAlign* pj = apa + j;
	        if (pi->qdir == pj->qdir && pi->qoff == pj->qoff && pi->qend == pj->qend && pi->sid == pj->sid && pi->soff == pj->soff && pi->send == pj->send) {
			    continue;
	        }
            if (pi->chain_qoff == pj->chain_qoff && pi->chain_qend == pj->chain_qend) {
                if (pi->sid == pj->sid && pi->soff <= pj->soff && pi->send > pj->soff) {
                    continue;
                }
                if (pi->sid == pj->sid && pj->soff <= pi->soff && pj->send > pi->soff) {
                    continue;
                }
                if (pj->pi > rep_pi) { rep_pi = pj->pi; rep_j = j; }
            }
        }
        double p1 = s_map_q_p1(chain_type);
        double p2 = fabs(1.0 - rep_pi / pi->pi);
        int map_q = ceil(-10 * log10(1.0 - p1 * p2));
	    if (map_q > 60) map_q = 60;
	    if (map_q < 0) map_q = 0;
        //dump_chain_pca(fprintf, stderr, *pi, i);
        //fprintf(stderr, "p1 = %g, pi2 = %g, mapq = %d\n", p1, rep_pi, map_q);
        pi->map_q = map_q;
	    pi->qc = rep_j;
    }

    s_set_map_q_for_hom_pca(pca_a, pca_c, chain_type, apa, apc);
    s_set_hom_pca_idx(pca_a, pca_c, all_pca_list);
}

///// end aux functions

static void
s_detect_hom_hits(HbnInitHit* hita, int hitc, int read_size)
{
    for (int i = 0; i < hitc; ++i) hita[i].is_hom = 0;
    for (int i = 0; i < hitc; ++i) {
        int iqb = hita[i].qbeg;
        int iqe = hita[i].qend;
        if (hita[i].qdir == REV) {
            iqb = read_size - hita[i].qend;
            iqe = read_size - hita[i].qbeg;
        }
        for (int j = i + 1; j < hitc; ++j) {
            int jqb = hita[j].qbeg;
            int jqe = hita[j].qend;
            if (hita[j].qdir == REV) {
                jqb = read_size - hita[j].qend;
                jqe = read_size - hita[j].qbeg;
            }     
            if (abs(iqb - jqb) <= 10 && abs(jqb - jqe) <= 10) {
                hita[i].is_hom = 1;
                hita[j].is_hom = 1;
            }       
        }
    }
}

static void
s_detect_candidates(MapThreadData* data, 
    const int qidx,
    vector<HbnInitHit>& hit_list,
    vector<SubjectInitHitInfo>& sbjct_hit_info_list)
{
    hit_list.clear();
    sbjct_hit_info_list.clear();

    const char* query_name = SeqReader_SeqName(data->queries, qidx);
    const u8* fwd_query = SeqReader_Seq(data->queries, qidx, FWD);
    const u8* rev_query = SeqReader_Seq(data->queries, qidx, REV);
    const char* query_qv = SeqReader_QV(data->queries, qidx);
    const int query_size = SeqReader_SeqSize(data->queries, qidx);
    WordFinderThreadData* word_data = data->word_data;
    hbn_ddfkm(word_data, qidx);
    kv_dinit(vec_init_hit, l_hit_list);
    SubjectInitHitInfo sbjct_hitinfo;
    vector<HbnKmerMatch> fwd_seeds, rev_seeds;

    for (size_t skmi_i = 0; skmi_i < kv_size(word_data->skmi_list); ++skmi_i) {
        SubjectKmInfo skmi = kv_A(word_data->skmi_list, skmi_i);
        const int sid = skmi.sid;
        DDFKmerMatch* kma = kv_data(word_data->ddfkm_list) + skmi.km_offset;
        const int kmc = skmi.km_cnt;
        const int subject_size = SeqReader_SeqSize(data->subjects, skmi.sid);
        fwd_seeds.clear();
        rev_seeds.clear();
        for (int k = 0; k < kmc; ++k) {
            HbnKmerMatch hkm;
            hkm.qoff = kma[k].qoff;
            hkm.soff = kma[k].soff;
            hkm.length = data->opts->kmer_size;
            if (kma[k].context & 1) {
                rev_seeds.push_back(hkm);
            } else {
                fwd_seeds.push_back(hkm);
            }
        }

        sbjct_hitinfo.sid = skmi.sid;
        sbjct_hitinfo.max_score = 0;
        sbjct_hitinfo.hit_offset = hit_list.size();
        sbjct_hitinfo.hit_cnt = 0;

        HbnKmerMatch* hkma = fwd_seeds.data();
        int hkmc = fwd_seeds.size();
        kv_clear(l_hit_list);
        make_candidate_kmer_chain(hkma,
            hkmc,
            qidx,
            FWD,
            query_size,
            sid,
            subject_size,
            data->opts->kmer_dist,
            data->opts->ddf,
            data->opts->chain_score,
            &l_hit_list);
        HbnInitHit* hit_array = kv_data(l_hit_list);
        int hit_count = kv_size(l_hit_list);
        for (int x = 0; x < hit_count; ++x) {
            HbnInitHit hit = kv_A(l_hit_list, x);
            hit_list.push_back(hit);
            sbjct_hitinfo.hit_cnt++;
            sbjct_hitinfo.max_score = hbn_max(sbjct_hitinfo.max_score, hit.score);
        }

        hkma = rev_seeds.data();
        hkmc = rev_seeds.size();
        kv_clear(l_hit_list);
        make_candidate_kmer_chain(hkma,
            hkmc,
            qidx,
            REV,
            query_size,
            sid,
            subject_size,
            data->opts->kmer_dist,
            data->opts->ddf,
            data->opts->chain_score,
            &l_hit_list);
        hit_array = kv_data(l_hit_list);
        hit_count = kv_size(l_hit_list);
        for (int x = 0; x < hit_count; ++x) {
            HbnInitHit hit = kv_A(l_hit_list, x);
            hit_list.push_back(hit);
            sbjct_hitinfo.hit_cnt++;
            sbjct_hitinfo.max_score = hbn_max(sbjct_hitinfo.max_score, hit.score);
        }

        if (sbjct_hitinfo.hit_cnt > 0) sbjct_hit_info_list.push_back(sbjct_hitinfo);
    }
    kv_destroy(l_hit_list);

    sort(sbjct_hit_info_list.begin(), sbjct_hit_info_list.end(),
        [](const SubjectInitHitInfo& a, const SubjectInitHitInfo& b)->bool { return a.max_score > b.max_score; });
    s_detect_hom_hits(hit_list.data(), hit_list.size(), query_size);
}

static void
s_extend_gapped_filling_hits(HbnTracebackData* tbck_data,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const int query_id,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    SeqReader* subjects,
    SubjectInitHitInfo* sihia,
    int sihic,
    HbnInitHit* hit_list,
    int* cov_stats,
    vector<PoreCAlign>& pca_list)
{
    fill(cov_stats, cov_stats + query_size, 0);
    for (auto& pca : pca_list) {
	bool r = (pca.lqm && pca.rqm) || (pca.lsm && pca.rsm);
	if (!r) continue;
        for (int x = pca.chain_qoff; x < pca.chain_qend; ++x) ++cov_stats[x];
    }

    vector<HbnInitHit> gf_hit_list;
    for (int i = 0; i < sihic; ++i) {
        HbnInitHit* hita = hit_list + sihia[i].hit_offset;
        int hitc = sihia[i].hit_cnt;
        int subject_id = hita[0].sid;
        const u8* subject = SeqReader_Seq(subjects, subject_id, FWD);
	    const int subject_size = SeqReader_SeqSize(subjects, subject_id);
        for (int p = 0; p < hitc; ++p) gf_hit_list.push_back(hita[p]);
    }
    sort(gf_hit_list.begin(), gf_hit_list.end(), [](const HbnInitHit& a, const HbnInitHit& b) { return a.score > b.score; });
    HbnInitHit* hita = gf_hit_list.data();
    int hitc = gf_hit_list.size();
    vector<PoreCAlign> l_pca_list;

    int extended_gf = 0;
    for (int p = 0; p < hitc; ++p) {
        if (hita[p].qdir == REV){
        	int qb = query_size - hita[p].qend;
        	int qe = query_size - hita[p].qbeg;
        	hita[p].qbeg = qb;
        	hita[p].qend = qe;
        }
        int cov = 0;
        for (int x = hita[p].qbeg; x < hita[p].qend; ++x) if (cov_stats[x]) ++cov;
        if (hita[p].is_hom) continue;
        if (extended_gf >= 20 && hita[p].qend - hita[p].qbeg - cov < 2) continue;
        ++extended_gf;
        l_pca_list.clear();
        int subject_id = hita[p].sid;
        const u8* subject = SeqReader_Seq(subjects, subject_id, FWD);
	    const int subject_size = SeqReader_SeqSize(subjects, subject_id);
        int qb = hita[p].qbeg;
        int qe = hita[p].qend;
        int sb = hita[p].sbeg;
        int se = hita[p].send;
        if (hita[p].qdir == REV) {
            qb = query_size - hita[p].qend;
            qe = query_size - hita[p].qbeg;
        }
        align_subseq(query_id, hita[p].qdir, fwd_query, rev_query, query_size,
                hita[p].sid, subject, subject_size, qb, qe, sb, se, hita[p].score,
                opts->perc_identity, tbck_data, reloci_list, qvep_list, l_pca_list);
        pca_list.insert(pca_list.end(), l_pca_list.begin(), l_pca_list.end());
    }
}

void
align_one_read(MapThreadData* data, const int qidx, kstring_t* out)
{
    const int verbose = 0;
    const char* query_name = SeqReader_SeqName(data->queries, qidx);
    const u8* fwd_query = SeqReader_Seq(data->queries, qidx, FWD);
    const u8* rev_query = SeqReader_Seq(data->queries, qidx, REV);
    const char* query_qv = SeqReader_QV(data->queries, qidx);
    const int query_size = SeqReader_SeqSize(data->queries, qidx);

    //HBN_LOG("mappig query %d:%s:%d", qidx, query_name, query_size);

    vector<HbnInitHit> hit_list;
    vector<SubjectInitHitInfo> sbjct_hit_info_list;
    s_detect_candidates(data, qidx, hit_list, sbjct_hit_info_list);
    SubjectInitHitInfo* shia = sbjct_hit_info_list.data();
    int shic = sbjct_hit_info_list.size();
    if (!shic) return;

    QueryVdfEndPointList_Setup(&data->qvep_list, &data->reloci_list->enzyme, fwd_query, rev_query, query_size);
    const int* vdfa = kv_data(data->qvep_list.fwd_vdf_endpoint_list);
    const int vdfc = kv_size(data->qvep_list.fwd_vdf_endpoint_list);
    const int enzyme_size = data->reloci_list->enzyme.enzyme_size;
    EChainType chain_type;
    vector<PoreCAlign> all_pca_list;
    vector<PoreCAlign> l_pca_list, l1_pca_list, chain, global_chain;
    int cov = 0;

    vector<int> _cov_stats(query_size); int* cov_stats = _cov_stats.data();
    for (int i = 0; i < shic && i < data->opts->hitlist_size; ++i) {
        HbnInitHit* hita = hit_list.data() + shia[i].hit_offset;
        int hitc = shia[i].hit_cnt;
        int sid = hita[0].sid;
        const char* sname = SeqReader_SeqName(data->subjects, sid);
        //HBN_LOG("mapping %d hits for subject %d:%s, score = %d", hitc, sid, sname, shia[i].max_score);
        l_pca_list.clear();
        extend_hit_list(data->tbck_data, 
            data->opts, 
            data->reloci_list, 
            &data->qvep_list, 
            data->subjects, 
            qidx,
            query_name,
            fwd_query,
            rev_query,
            query_qv,
            query_size,
            hita,
            hitc,
            data->opts->max_hsps_per_subject,
            cov_stats,
            l_pca_list);
	    l1_pca_list.assign(l_pca_list.begin(), l_pca_list.end());
        remove_contained_pca(l1_pca_list);
        all_pca_list.insert(all_pca_list.end(), l1_pca_list.begin(), l1_pca_list.end());
        if (select_pca_chain(data->pca_chain_data, vdfa, vdfc, l1_pca_list.data(), l1_pca_list.size(), ePerfectChain, chain, &cov)) {
            if (verbose) {HBN_LOG("1 find perfect chain, cov = %d", cov);for (auto& xp : chain) dump_chain_pca(fprintf, stderr, xp, -1);}
            chain_type = ePerfectChain;
            smooth_pca_list(chain, chain_type, all_pca_list, data->subjects, query_name, fwd_query, rev_query, vdfa, vdfc, enzyme_size, data->tbck_data);
            s_set_map_q(chain_type, chain.data(), chain.size(), all_pca_list);
            trim_overlap_subseqs(data->tbck_data, data->reloci_list, &data->qvep_list, data->subjects, query_name, qidx, fwd_query, rev_query, query_qv, query_size,
                all_pca_list.data(), all_pca_list.size(), chain.data(), chain.size(), chain_type, 0, data->opts->outfmt, out);
            return;
        }
    }

    if (shic > data->opts->hitlist_size) {
        s_extend_gapped_filling_hits(data->tbck_data,
            data->opts,
            data->reloci_list,
            &data->qvep_list,
            qidx,
            query_name,
            fwd_query,
            rev_query,
            query_size,
            data->subjects,
            shia + data->opts->hitlist_size,
            shic - data->opts->hitlist_size,
            hit_list.data(),
            cov_stats,
            all_pca_list);
    }
    l1_pca_list.assign(all_pca_list.begin(), all_pca_list.end());
    if (select_pca_chain(data->pca_chain_data, vdfa, vdfc, l1_pca_list.data(), l1_pca_list.size(), ePerfectChain, chain, &cov)) {
        if (verbose) { HBN_LOG("3 find perfect chain, cov = %d", cov); for (auto& xp : chain) dump_chain_pca(fprintf, stderr, xp, -1); }
        chain_type = ePerfectChain;
        smooth_pca_list(chain, chain_type, all_pca_list, data->subjects, query_name, fwd_query, rev_query, vdfa, vdfc, enzyme_size, data->tbck_data);
        s_set_map_q(chain_type, chain.data(), chain.size(), all_pca_list);
        trim_overlap_subseqs(data->tbck_data, data->reloci_list, &data->qvep_list, data->subjects, query_name, qidx, fwd_query, rev_query, query_qv, query_size,
            all_pca_list.data(), all_pca_list.size(), chain.data(), chain.size(), chain_type, 0, data->opts->outfmt, out);
        return;
    }

    //// find pseudo perfect chain on single chromosome
    PoreCAlign* pca_a = nullptr;
    int pca_c = 0;
    pca_a = all_pca_list.data();
    pca_c = all_pca_list.size();
    sort(pca_a, pca_a + pca_c, [](const PoreCAlign& a, const PoreCAlign& b) { return a.sid < b.sid; });
    vector<PoreCAlign> best_chain;
    int max_cov = 0;
    int i = 0;
    while (i < pca_c) {
        int j = i + 1;
        while (j < pca_c && pca_a[i].sid == pca_a[j].sid) ++j;
        if (!select_pca_chain(data->pca_chain_data, vdfa, vdfc, pca_a + i, j - i, ePseudoPerfectChain, chain, &cov)) {
            i = j;
            continue;
        }
        if (verbose) { HBN_LOG("1 find pseudo perfect chain, cov = %d", cov); for (auto& xp : chain) dump_chain_pca(fprintf, stderr, xp, -1); }
        if (cov > max_cov) {
            max_cov = cov;
            best_chain.assign(chain.begin(), chain.end());
        }
        i = j;
    }

    l1_pca_list.assign(all_pca_list.begin(), all_pca_list.end());
    if (!select_pca_chain(data->pca_chain_data, vdfa, vdfc, l1_pca_list.data(), l1_pca_list.size(), ePseudoPerfectChain, global_chain, &cov)) {
        global_chain.clear();
    }

    if (max_cov) {
        if (!global_chain.empty()) {
            double g_avg_pi = compute_pca_list_avg_pi(global_chain.data(), global_chain.size());
            double avg_pi = compute_pca_list_avg_pi(best_chain.data(), best_chain.size());
            if (g_avg_pi > avg_pi) {
                chain_type = eCompleteChain;
                smooth_pca_list(global_chain, chain_type, all_pca_list, data->subjects, query_name, fwd_query, rev_query, vdfa, vdfc, enzyme_size, data->tbck_data);
                s_set_map_q(chain_type, global_chain.data(), global_chain.size(), all_pca_list);
                trim_overlap_subseqs(data->tbck_data,
                    data->reloci_list,
                    &data->qvep_list,
                    data->subjects,
                    query_name,
                    qidx,
                    fwd_query,
                    rev_query,
                    query_qv,
                    query_size,
                    all_pca_list.data(),
                    all_pca_list.size(),
                    global_chain.data(),
                    global_chain.size(),
                    chain_type,
                    0,
                    data->opts->outfmt,
                    out); 
                return;               
            }
        }
        chain_type = eCompleteChain;
        smooth_pca_list(best_chain, chain_type, all_pca_list, data->subjects, query_name, fwd_query, rev_query, vdfa, vdfc, enzyme_size, data->tbck_data);
        s_set_map_q(chain_type, best_chain.data(), best_chain.size(), all_pca_list);
        trim_overlap_subseqs(data->tbck_data,
                data->reloci_list,
                &data->qvep_list,
                data->subjects,
                query_name,
                qidx,
                fwd_query,
                rev_query,
                query_qv,
                query_size,
                all_pca_list.data(),
                all_pca_list.size(),
                best_chain.data(),
                best_chain.size(),
                chain_type,
                0,
                data->opts->outfmt,
                out);
        return;
    }

    l1_pca_list.assign(all_pca_list.begin(), all_pca_list.end());
    if (select_pca_chain(data->pca_chain_data, vdfa, vdfc, l1_pca_list.data(), l1_pca_list.size(), eCompleteChain, chain, &cov)) {
        if (verbose) { HBN_LOG("2 find complete chain, cov = %d", cov); for (auto& xp : chain) dump_chain_pca(fprintf, stderr, xp, -1); }
        chain_type = eCompleteChain;
        smooth_pca_list(chain, chain_type, all_pca_list, data->subjects, query_name, fwd_query, rev_query, vdfa, vdfc, enzyme_size, data->tbck_data);
        s_set_map_q(chain_type, chain.data(), chain.size(), all_pca_list);
        trim_overlap_subseqs(data->tbck_data, data->reloci_list, &data->qvep_list, data->subjects, query_name, qidx, fwd_query, rev_query, query_qv, query_size,
            all_pca_list.data(), all_pca_list.size(), chain.data(), chain.size(), chain_type, 0, data->opts->outfmt, out);
        return;
    }

    l1_pca_list.assign(all_pca_list.begin(), all_pca_list.end());
    if (select_pca_chain(data->pca_chain_data, vdfa, vdfc, l1_pca_list.data(), l1_pca_list.size(), eMaxCovChain, chain, &cov)) {
        if (verbose) { HBN_LOG("2 find max cov chain, cov = %d", cov); for (auto& xp : chain) dump_chain_pca(fprintf, stderr, xp, -1); }   
        chain_type = eMaxCovChain;
        smooth_pca_list(chain, chain_type, all_pca_list, data->subjects, query_name, fwd_query, rev_query, vdfa, vdfc, enzyme_size, data->tbck_data);
        s_set_map_q(chain_type, chain.data(), chain.size(), all_pca_list);
        trim_overlap_subseqs(data->tbck_data, data->reloci_list, &data->qvep_list, data->subjects, query_name, qidx, fwd_query, rev_query, query_qv, query_size,
            all_pca_list.data(), all_pca_list.size(), chain.data(), chain.size(), chain_type, 0, data->opts->outfmt, out);
        return;
    }
}
