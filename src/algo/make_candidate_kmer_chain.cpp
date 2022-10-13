#include "make_candidate_kmer_chain.h"

#include <algorithm>
#include <vector>

#include <cmath>

using namespace std;

void
validate_hkm_list(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const HbnKmerMatch* hkma,
    const int hkmc)
{
#if 1
    for (int i = 0; i < hkmc; ++i) {
        int qi = hkma[i].qoff;
        int si = hkma[i].soff;
        for (int k = 0; k < hkma[k].length; ++k, ++qi, ++si) {
            hbn_assert(read[qi] == subject[si], "[%s, %s, %d] at (%d, %d, %d): qi = %d, si = %d, q = %d, s = %d", 
            HBN_LOG_ARGS_GENERIC, hkma[k].qoff, hkma[k].soff, hkma[k].length, qi, si, read[qi], subject[si]);
        }
    }
#endif
}

static BOOL
s_init_hit_is_contained(HbnInitHit* ha, int hc, HbnInitHit* hn)
{
    const int E = 10;
    for (int i = 0; i < hc; ++i) {
        int r = (hn->qbeg + E >= ha[i].qbeg)
                &&
                (hn->qend <= ha[i].qend + E)
                &&
                (hn->sbeg + E >= ha[i].sbeg)
                &&
                (hn->send <= ha[i].send + E);
        if (r) return TRUE;
    }
    return FALSE;
}

static BOOL 
s_cov_is_full(int* cov_stats, int qb, int qe)
{
    int cov = 0;
    for (int i = qb; i < qe; ++i) {
        if (cov_stats[i] < 10) ++cov;
    }
    return cov < 20;
}

extern "C"
void
make_candidate_kmer_chain(const HbnKmerMatch* hkma,
    const int hkmc,
    const int query_id,
    const int query_dir,
    const int query_size,
    const int subject_id,
    const int subject_size,
    const int max_kmer_dist,
    const double max_ddf,
    const int min_chain_score,
    vec_init_hit* init_hit_list)
{
    vector<int> _pred_list(hkmc); int* pred_list = _pred_list.data();
    vector<int> _score_list(hkmc); int* score_list = _score_list.data();
    fill(pred_list, pred_list + hkmc, -1);
    fill(score_list, score_list + hkmc, 1);

    for (int i = 0; i < hkmc; ++i) {
        int max_score = 1;
        int max_j = -1;
        int iqb = hkma[i].qoff;
        int isb = hkma[i].soff;
        for (int j = i - 1; j >= 0; --j) {
            int jqb = hkma[j].qoff;
            int jsb = hkma[j].soff;
            hbn_assert(jsb <= isb);
            if (isb - jsb > max_kmer_dist) break;
            if (jsb == isb || jqb >= iqb || iqb - jqb > max_kmer_dist) continue;
            int q_d = iqb - jqb;
            int s_d = isb - jsb;
            int d_d = abs(q_d - s_d); // distance difference
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > max_ddf) continue;
            if (d_d > 40) continue;
            int score = score_list[j] + 1;
            if (score > max_score) {
                max_score = score;
                max_j = j;
            } 
        }
        score_list[i] = max_score;
        pred_list[i] = max_j;
    }

    vector<int> _succ_list(hkmc); int* succ_list = _succ_list.data();
    fill(succ_list, succ_list + hkmc, 0);
    for (int i = 0; i < hkmc; ++i) {
        if (pred_list[i] >= 0) succ_list[pred_list[i]] = 1;
    }
    vector<pair<int, int>> chain_score_list;
    for (int i = 0; i < hkmc; ++i) {
        if (succ_list[i] == 0 && score_list[i] >= min_chain_score) {
            chain_score_list.emplace_back(i, score_list[i]);
        }
    }
    sort(chain_score_list.begin(), 
         chain_score_list.end(),
         [](const pair<int, int>& a, const pair<int, int>& b)->bool {
             return (a.second > b.second) || (a.second == b.second && a.first < b.first);
         });

    int* avail_list = succ_list;
    fill(avail_list, avail_list + hkmc, 1);
    vector<int> idx_list;
    HbnInitHit init_hit; memset(&init_hit, 0, sizeof(HbnInitHit));
    vector<HbnInitHit> l_hit_list;
    vector<int> cov_stats(query_size); fill(cov_stats.begin(), cov_stats.end(), 0);
    for (auto& ias : chain_score_list) {
        if (!avail_list[ias.first]) continue;
        int p = ias.first;
        idx_list.clear();
        while (p >= 0) {
            idx_list.push_back(p);
            p = pred_list[p];
        }
        reverse(idx_list.begin(), idx_list.end());
        const HbnKmerMatch* fhkm = hkma + idx_list.front();
        const HbnKmerMatch* lhkm = hkma + idx_list.back();
        init_hit.qbeg = fhkm->qoff;
        init_hit.qend = lhkm->qoff + lhkm->length;
        init_hit.sbeg = fhkm->soff;
        init_hit.send = lhkm->soff + lhkm->length;
        init_hit.score = ias.second;
        if (s_init_hit_is_contained(l_hit_list.data(), l_hit_list.size(), &init_hit)) continue;
        if (s_cov_is_full(cov_stats.data(), init_hit.qbeg, init_hit.qend)) continue;
        for (int i = init_hit.qbeg; i < init_hit.qend; ++i) ++cov_stats[i];
        l_hit_list.push_back(init_hit);
        for (auto& i : idx_list) avail_list[i] = 0;
        if (l_hit_list.size() == 200) break;
    }

    for (auto& hit : l_hit_list) {
        hit.qid = query_id;
        hit.qdir = query_dir;
        hit.qsize = query_size;
        hit.sid = subject_id;
        hit.sdir = FWD;
        hit.ssize = subject_size;
        kv_push(HbnInitHit, *init_hit_list, hit);
    }
}