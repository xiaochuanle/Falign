#include "chain_align_list.hpp"
#include "../../corelib/pdqsort.h"

using namespace std;

#include <array>
#include <vector>

using namespace std;

const char* get_chain_type_name(EChainType type)
{
    static const char* name[] = { "perfect", "complete", "complete", "incomplete" };
    return name[type];
}

BOOL set_pca_chain_offset(PoreCAlign* pca, const int enzyme_size)
{
    pca->lqm = (pca->qoff == 0) || (pca->qoff == pca->enzyme_qoff);
    pca->rqm = (pca->qend == pca->qsize) || (pca->qend == pca->enzyme_qend);
    pca->lsm = (pca->soff == 0) || (pca->soff == pca->enzyme_soff);
    pca->rsm = (pca->send == pca->ssize) || (pca->send == pca->enzyme_send);
    
    bool l_perfect = (pca->lqm && pca->lsm) || (pca->qoff == 0) || (pca->soff == 0);
    bool r_perfect = (pca->rqm && pca->rsm) || (pca->qend == pca->qsize) || (pca->send == pca->ssize);
    pca->is_perfect = l_perfect && r_perfect;
    pca->lp = l_perfect;
    pca->rp = r_perfect;

    bool l_pseudo_perfect = l_perfect || pca->lqm || pca->lsm;
    bool r_pseudo_perfect = r_perfect || pca->rqm || pca->rsm;
    pca->is_pseudo_perfect = l_pseudo_perfect && r_pseudo_perfect;
    pca->lpp = l_pseudo_perfect;
    pca->rpp = r_pseudo_perfect;

    pca->qc = pca->qend - pca->qoff;
    if (pca->qdir == FWD) {
        pca->chain_qoff = pca->qoff;
        pca->chain_qend = pca->qend;
        pca->lcp = pca->lp;
        pca->lcpp = pca->lpp;
        pca->rcp = pca->rp;
        pca->rcpp = pca->rpp;
    } else {
        pca->chain_qoff = pca->qsize - pca->qend;
        pca->chain_qend = pca->qsize - pca->qoff;
        pca->lcp = pca->rp;
        pca->lcpp = pca->rpp;
        pca->rcp = pca->lp;
        pca->rcpp = pca->lpp;

        if (pca->qend == pca->enzyme_qend && pca->chain_qoff >= enzyme_size) {
                pca->chain_qoff -= enzyme_size;
        }
        if (pca->chain_qend < pca->qsize && pca->qoff == pca->enzyme_qoff) {
                pca->chain_qend -= enzyme_size;
        }
    }
    if (pca->chain_qoff >= pca->chain_qend) return FALSE;
    if (pca->soff >= pca->send || pca->enzyme_soff > pca->enzyme_send) return FALSE;
    return TRUE;
}

static void
s_remove_contained_pca(PoreCAlign* pca_a, int pca_c, PoreCAlign* all_pca_a, int all_pca_c)
{
    if (pca_c == 1) return;

    for (int i = 0; i < pca_c; ++i) {
        if (pca_a[i].qid == -1) continue;
        PoreCAlign* pi = pca_a + i;
        for (int j = 0; j < all_pca_c; ++j) {
            PoreCAlign* pj = all_pca_a + j;
            if (pj->chain_qoff > pi->chain_qend) break;
            if (pi->qdir != pj->qdir || pi->chain_qend != pj->chain_qoff) continue;
            for (int k = i + 1; k < pca_c; ++k) {
                PoreCAlign* pk = pca_a + k;
                if (pj->qdir != pk->qdir || pj->chain_qend != pk->chain_qend) continue;
                HBN_LOG("find contained pca pair");
                dump_chain_pca(fprintf, stderr, *pi, -1);
                dump_chain_pca(fprintf, stderr, *pk, -1);
                dump_chain_pca(fprintf, stderr, *pj, -1);
                pi->qid = -1;
                pk->qid = -1;
                break;
            }
            if (pi->qid == -1) break;
        }
    }
}

void remove_contained_pca(PoreCAlign* pca_a, int* _pca_c)
{
    int pca_c = *_pca_c;
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pi = pca_a + i;
        if (pi->qid == -1) continue;
        for (int j = i + 1; j < pca_c; ++j) {
            PoreCAlign* pj = pca_a + j;
            if (pj->qid == -1) continue;
            if (pi->chain_qoff == pj->chain_qoff && pi->chain_qend == pj->chain_qend && pi->soff == pj->soff && pi->send == pj->send) {
                if (pi->pi >= pj->pi) {
                    pj->qid = -1;
                } else {
                    *pi = *pj;
                    pj->qid = -1;
                }
            }
        }
    }
    int n = 0;
    for (int i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
    pca_c = n;

    pdqsort(pca_a, pca_a + pca_c, 
        [](const PoreCAlign& a, const PoreCAlign& b) { return (a.chain_qoff < b.chain_qoff) || (a.chain_qoff == b.chain_qoff && a.chain_qend > b.chain_qend); });
    int i = 0;
    while (i < pca_c) {
        int j = i + 1;
        while (j < pca_c && pca_a[i].chain_qoff == pca_a[j].chain_qoff) ++j;
        s_remove_contained_pca(pca_a + i, j - i, pca_a + j, pca_c - j);
        i = j;
    }
    n = 0;
    for (i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
    pca_c = n;

    *_pca_c = n;    
}

void remove_contained_pca(std::vector<PoreCAlign>& pca_list)
{
    PoreCAlign* pca_a = pca_list.data();
    int pca_c = pca_list.size();
    remove_contained_pca(pca_a, &pca_c);
    pca_list.resize(pca_c);
}

void remove_duplicate_pca(PoreCAlign* pca_a, int* pca_c_)
{
    int pca_c = *pca_c_;
    pdqsort(pca_a, pca_a + pca_c,
        [](const PoreCAlign& a, const PoreCAlign& b)->bool {
            return a.qc > b.qc;
    });
    const int E = 10;

    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pi = pca_a + i;
        if (pi->qid == -1) continue;
        for (int j = i + 1; j < pca_c; ++j) {
            PoreCAlign* pj = pca_a + j;
            if (pj->qid == -1) continue;
            if (pj->qdir != pi->qdir) continue;
            int r = (abs(pi->qoff - pj->qoff) <= E)
                    &&
                    (abs(pi->qend - pj->qend) <= E)
                    &&
                    (abs(pi->soff - pj->soff) <= E)
                    &&
                    (abs(pi->send - pj->send) <= E);
            if (!r) continue;
            if (pi->is_perfect) {
                //HBN_LOG("1 remove contained pca");
                //dump_pca(fprintf, stderr, *pj, j);
                pj->qid = -1;
                continue;
            }
            if (pj->is_perfect) {
                //HBN_LOG("2 remove contained pca");
                //dump_pca(fprintf, stderr, *pi, i);
                *pi = *pj;
                pj->qid = -1;
                break;
            }
                //HBN_LOG("3 remove contained pca");
                //dump_pca(fprintf, stderr, *pj, j);
            pj->qid = -1;
        }
    }

    int n = 0;
    for (int i = 0; i < pca_c; ++i) {
        if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
    }
    *pca_c_ = n;
}

static BOOL 
s_is_perfect_chain_offset(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    if (al->is_perfect == false || ar->is_perfect == false) return FALSE;
    if (al->chain_qend != ar->chain_qoff) return FALSE;
    if (al->sid == ar->sid) {
        if (al->soff == ar->soff || al->send == ar->send) {
            bool r = (ar->soff >= al->soff && ar->send <= al->send)
                 ||
                 (al->soff >= ar->soff && al->send <= ar->send);
            if (r) return TRUE;
        }        
        if (al->soff <= ar->soff && al->send > ar->soff) return FALSE;
        if (ar->soff <= al->soff && ar->send > al->soff) return FALSE;
    }
    return TRUE;
}

static BOOL
s_is_validate_soff_relation_general(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    if (al->sid != ar->sid) return TRUE;

    const int E = 10;
    if (abs(al->soff - ar->soff) <= E || abs(al->send - ar->send) <= E) {
        bool r = (ar->soff + E >= al->soff && ar->send <= al->send + E)
                 ||
                 (al->soff + E >= ar->soff && al->send <= ar->send + E);
        if (r) return true;
    }

    if (al->soff <= ar->soff) {
        if (al->send > ar->soff) {
            if (al->send - ar->soff > max_align_ovlp) {
                return FALSE;
            }
        }
    } else {
        if (ar->send > al->soff) {
            if (ar->send - al->soff > max_align_ovlp) {
                return FALSE;
            }
        }
    }

    return TRUE;
}

static BOOL
s_is_pseudo_perfect_chain_offset(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    hbn_assert(al->chain_qoff <= ar->chain_qoff);
    if (ar->chain_qend - al->chain_qend <= 5 || ar->chain_qoff - al->chain_qoff <= 5) return FALSE;
    if (al->is_pseudo_perfect == false || ar->is_pseudo_perfect == false) return FALSE;

    if (al->rcp && ar->lcp && al->chain_qend != ar->chain_qoff) return FALSE;

    if (al->rcp && al->chain_qend > ar->chain_qoff) return FALSE;

    if (ar->lcp && al->chain_qend > ar->chain_qoff) return FALSE;

    if (al->chain_qend > ar->chain_qoff && al->chain_qend - ar->chain_qoff > max_align_ovlp) {
        return FALSE;
    }

    return s_is_validate_soff_relation_general(al, ar, max_align_ovlp);
}

BOOL 
is_validate_chain_offset_general(PoreCAlign* al, PoreCAlign* ar, int max_align_ovlp)
{
    hbn_assert(al->chain_qoff <= ar->chain_qoff);
    if (ar->chain_qend - al->chain_qend <= 5 || ar->chain_qoff - al->chain_qoff <= 5) return FALSE;
    
    if (al->rcp && ar->lcp) {
        if (al->chain_qend > ar->chain_qoff) return FALSE;
        //if (ar->chain_qoff > al->chain_qend && ar->chain_qoff - al->chain_qend > 50) return FALSE;
        if (al->sid == ar->sid) {
            if (al->soff <= ar->soff && al->send > ar->soff) return FALSE;
            if (ar->soff <= al->soff && ar->send > al->soff) return FALSE;
        }
        return TRUE;
    }
    
    if (al->chain_qend > ar->chain_qoff && al->chain_qend - ar->chain_qoff > max_align_ovlp) {
        return FALSE;
    }

    return s_is_validate_soff_relation_general(al, ar, max_align_ovlp);
}

static int
s_pca_list_relation(vector<int>& a_idx_list, vector<int>& b_idx_list)
{
    int* l_a = NULL;
    int l_c = 0;
    int* s_a = NULL;
    int s_c = 0;
    bool reverse_relation = false;
    if (a_idx_list.size() >= b_idx_list.size()) {
        l_a = a_idx_list.data();
        l_c = a_idx_list.size();
        s_a = b_idx_list.data();
        s_c = b_idx_list.size();
    } else {
        l_a = b_idx_list.data();
        l_c = b_idx_list.size();
        s_a = a_idx_list.data();
        s_c = a_idx_list.size();
        reverse_relation = true;
    }

    int l_i = 0, s_i = 0, n_eq = 0;
    while (l_i < l_c && s_i < s_c) {
        while (l_i < l_c && l_a[l_i] < s_a[s_i]) ++l_i;
        if (l_i >= l_c) break;
        while (s_i < s_c && s_a[s_i] < l_a[l_i]) ++s_i;
        if (s_i >= s_c) break;
        if (l_a[l_i] == s_a[s_i]) {
            ++n_eq;
            ++l_i;
            ++s_i;
        }
    }

    if (n_eq != s_c) return 0;
    if (reverse_relation) return -1;
    return 1;
}

static bool 
s_is_large_gap(const int* vdfa, const int vdfc, int from, int to)
{
    int i = 0;
    for (; i < vdfc; ++i) if (from >= vdfa[i]) break;
    if (i >= vdfc - 1) return false;
    return to <= vdfa[i+1];
}

bool s_is_complete_chain(const int* vdfa, const int vdfc, const PoreCAlign* pca_a, const int pca_c)
{
    if (vdfc <= 2) return true;

    if (pca_a[0].chain_qoff > 50) return false;
    /// left gap
    {
        int l_gap = 0;
        for (int i = 1; i < vdfc; ++i) {
            if (vdfa[i] > pca_a[0].chain_qoff) break;
            l_gap = vdfa[i];
        }
        //HBN_LOG("l_gap = %d", l_gap);
        if (l_gap > 50) return false;
    }

    if (pca_a[pca_c-1].qsize - pca_a[pca_c-1].chain_qend > 50) return false;
    /// right gap
    {
        int qsize = pca_a[0].qsize;
        int r_gap = 0;
        for (int i = vdfc - 3; i >= 0; --i) {
            if (vdfa[i] < pca_a[pca_c-1].chain_qend) break;
            r_gap = qsize - vdfa[i];
        }
        //HBN_LOG("r_gap = %d", r_gap);
        if (r_gap > 50) return false;
    }

    int gap_size = pca_a[0].chain_qoff;
    for (int i = 0; i < pca_c - 1; ++i) {
        int gap = pca_a[i+1].chain_qoff - pca_a[i].chain_qend;
        if (gap > 50) return false;
        if (gap > 0) gap_size += gap;
    }
    gap_size += (pca_a[pca_c-1].qsize - pca_a[pca_c-1].chain_qend);

    return (gap_size < kMaxMissedCovBase);
}

bool s_is_perfect_chain(PoreCAlign* pca_a, int pca_c)
{
    for (int i = 0; i < pca_c - 1; ++i) {
        if (!s_is_perfect_chain_offset(pca_a + i, pca_a + i + 1, PC_MAX_ALIGN_OVLP)) return FALSE;
    }
    return TRUE;
}

bool s_is_pseudo_perfect_chain(PoreCAlign* pca_a, int pca_c)
{
    for (int i = 0; i < pca_c - 1; ++i) {
        if (!s_is_pseudo_perfect_chain_offset(pca_a + i, pca_a + i + 1, PC_MAX_ALIGN_OVLP)) return FALSE;
    }
    return TRUE;
}

EChainType pca_chain_type(const int* vdfa, const int vdfc, PoreCAlign* pca_a, int pca_c)
{
    if (s_is_perfect_chain(pca_a, pca_c)) return ePerfectChain;
    if (s_is_pseudo_perfect_chain(pca_a, pca_c)) return ePseudoPerfectChain;
    if (s_is_complete_chain(vdfa, vdfc, pca_a, pca_c)) return eCompleteChain;
    return eMaxCovChain;
}

int
compute_pca_list_q_cov(PoreCAlign* pca_a, int pca_c)
{
    int q_cov = 0;
    for (int i = 0; i < pca_c; ++i) {
        q_cov += (pca_a[i].chain_qend - pca_a[i].chain_qoff);
    }
    int ovlp = 0;
    for (int i = 0; i < pca_c - 1; ++i) {
	hbn_assert(pca_a[i].chain_qoff < pca_a[i+1].chain_qoff);
        if (pca_a[i].chain_qend > pca_a[i+1].chain_qoff) {
            ovlp += (pca_a[i].chain_qend - pca_a[i+1].chain_qoff);
        }
    }
    hbn_assert(q_cov > ovlp);
    return q_cov - ovlp;
}

double 
compute_pca_list_avg_pi(PoreCAlign* pca_a, int pca_c)
{
    if (pca_c < 1) return 0.0;
    double sum = 0.0;
    for (int i = 0; i < pca_c; ++i) sum += pca_a[i].pi;
    return sum / pca_c;
}

static int
s_fill_left_most_gap(PoreCAlign* a0, PoreCAlign* pca_a, int pca_c)
{
  //  HBN_LOG("left most pca"); dump_chain_pca(fprintf, stderr, *a0, -1);
    int best_pi = 0.0;
    int best_i = -1;
    int max_cov = 0;
    for (int i = 0; i < pca_c; ++i) {
        //HBN_LOG("examine pca"); dump_chain_pca(fprintf, stderr, pca_a[i], i);
        if (pca_a[i].chain_qoff >= a0->chain_qend) break;
        if (pca_a[i].chain_qoff >= a0->chain_qoff) break;
        if (pca_a[i].chain_qoff == a0->chain_qoff && pca_a[i].chain_qend == a0->chain_qend) continue;
        if (is_validate_chain_offset_general(pca_a + i, a0, PC_MAX_ALIGN_OVLP)) {
	    int cov = pca_a[i].chain_qend - pca_a[i].chain_qoff;
	    if (pca_a[i].chain_qend > a0->chain_qoff) cov -= (pca_a[i].chain_qend - a0->chain_qoff);
	    if (cov == max_cov) {
                if (pca_a[i].pi > best_pi) {
                    best_pi = pca_a[i].pi;
                    best_i = i;
                }
            } else if (cov > max_cov) {
               // HBN_LOG("======== find filling pca"); dump_chain_pca(fprintf, stderr, pca_a[i], i);
		best_pi = pca_a[i].pi;
		best_i = i;
		max_cov = cov;
	    }
        }
    }
    return best_i;
}

static int
s_fill_right_most_gap(PoreCAlign* an, PoreCAlign* pca_a, int pca_c)
{
    int best_pi = 0.0;
    int best_i = -1;
    int max_cov = 0;
    for (int i = 0; i < pca_c; ++i) {
        if (an->chain_qoff >= pca_a[i].chain_qoff) continue;
        if (pca_a[i].chain_qoff == an->chain_qoff && pca_a[i].chain_qend == an->chain_qend) continue;
        if (is_validate_chain_offset_general(an, pca_a + i, PC_MAX_ALIGN_OVLP)) {
	    int cov = pca_a[i].chain_qend - pca_a[i].chain_qoff;
	    if (an->chain_qend > pca_a[i].chain_qoff) cov -= (an->chain_qend - pca_a[i].chain_qoff);
	    if (cov == max_cov) {
               	if (pca_a[i].pi > best_pi) {
                    best_pi = pca_a[i].pi;
                    best_i = i;
                }
	    } else if (cov > max_cov) {
		best_pi = pca_a[i].pi;
		best_i = i;
		max_cov = cov;
	    }
        }
    }
    return best_i;
}

static int
s_fill_one_gap(PoreCAlign* al, PoreCAlign* ar, PoreCAlign* pca_a, int pca_c)
{
    hbn_assert(al->chain_qend < ar->chain_qoff);
    const int verbose = 0;
    if (verbose) {
        HBN_LOG("filling gaps between two pca");
        dump_chain_pca(fprintf, stderr, *al, -1);
        dump_chain_pca(fprintf, stderr, *ar, -1);
    }

    int from = -1;
    for (int i = 0; i < pca_c; ++i) {
        if (pca_a[i].chain_qoff == al->chain_qoff && pca_a[i].chain_qend == al->chain_qend) continue;
        if (pca_a[i].chain_qoff == ar->chain_qoff && pca_a[i].chain_qend == ar->chain_qend) continue;
        if (pca_a[i].chain_qoff - al->chain_qend > PC_MAX_ALIGN_OVLP) break;
        if (al->chain_qend - pca_a[i].chain_qoff <= PC_MAX_ALIGN_OVLP) {
            from = i;
            break;
        }
    }
    if (from == -1) return -1;

    int to = from;
    for (; to < pca_c; ++to) {
        if (pca_a[to].chain_qoff - al->chain_qend > PC_MAX_ALIGN_OVLP) break;
    }

    if (verbose) {
        HBN_LOG("find candidate gapped filling pca");
        for (int i = from; i < to; ++i) {
            dump_chain_pca(fprintf, stderr, pca_a[i], i);
        }
    }

    int best_i = -1;
    double max_pi = 0.0;
    /// perfect pca
    for (int i = from; i < to; ++i) {
        if (pca_a[i].chain_qoff == al->chain_qend && pca_a[i].chain_qend == ar->chain_qoff) {
            if (pca_a[i].pi > max_pi) {
                best_i = i;
                max_pi = pca_a[i].pi;
            }
        }
    }
    if (best_i != -1) {
        if (verbose) {
            HBN_LOG("found perfect gapped filling pca");
            dump_chain_pca(fprintf, stderr, pca_a[best_i], best_i);
        }
        return best_i;
    }

    /// one-side perfect pca
    for (int i = from; i < to; ++i) {
        PoreCAlign* pi = pca_a + i;
        if (pi->chain_qoff != al->chain_qend && pi->chain_qend != ar->chain_qoff) continue;
        if (al->chain_qoff > pi->chain_qoff) continue;
        if (pi->chain_qoff > ar->chain_qoff) continue;
        if (al->chain_qend == pi->chain_qoff && pi->chain_qend < ar->chain_qend) {
            if (pi->chain_qoff <= ar->chain_qoff && is_validate_chain_offset_general(pi, ar, PC_MAX_ALIGN_OVLP)) {
                if (pi->pi > max_pi) {
                    best_i = i;
                    max_pi = pi->pi;
                }
            }
        } else if (pi->chain_qend == ar->chain_qoff && pi->chain_qoff > al->chain_qoff) {
            if (al->chain_qoff <= pi->chain_qoff && is_validate_chain_offset_general(al, pi, PC_MAX_ALIGN_OVLP)) {
                if (pi->pi > max_pi) {
                    best_i = i;
                    max_pi = pi->pi;
                }
            }
        }
    }
    if (best_i != -1) {
        if (verbose) {
            HBN_LOG("found one-side perfect gapped filling pca");
            dump_chain_pca(fprintf, stderr, pca_a[best_i], best_i);
        }
        return best_i;
    }

    /// imperfect pca
    for (int i = from; i < to; ++i) {
        PoreCAlign* pi = pca_a + i;
        if (!(pi->chain_qoff > al->chain_qoff && pi->chain_qend < ar->chain_qend)) continue;
        if (al->chain_qoff < pi->chain_qoff && is_validate_chain_offset_general(al, pi, PC_MAX_ALIGN_OVLP)
            &&
            pi->chain_qoff < ar->chain_qoff && is_validate_chain_offset_general(pi, ar, PC_MAX_ALIGN_OVLP)) {
            if (pi->pi > max_pi) {
                best_i = i;
                max_pi = pi->pi;
            }
        }
    }
    if (best_i != -1) {
        if (verbose) {
            HBN_LOG("found imperfect gapped filling pca");
            dump_chain_pca(fprintf, stderr, pca_a[best_i], best_i);
        }
        return best_i;
    }

    return -1;
}

static void
s_fill_pca_list_gaps(PoreCAlign* all_pca_a, 
    int all_pca_c, 
    const int* vdfa, const int vdfc, 
    vector<PoreCAlign>& pca_list, 
    PoreCAlignChain* chain_info)
{
    vector<PoreCAlign> new_pca_list;
    PoreCAlign* pca_a = pca_list.data();
    int pca_c = pca_list.size();
    int find_gapped_fill_pca = 0;
    if (pca_a[0].chain_qoff > 50) {
        int r = s_fill_left_most_gap(pca_a, all_pca_a, all_pca_c);
        if (r >= 0) {
            new_pca_list.push_back(all_pca_a[r]);
            find_gapped_fill_pca = 1;
        }
    }
    for (int i = 0; i < pca_c - 1; ++i) {
        new_pca_list.push_back(pca_a[i]);
        if (pca_a[i+1].chain_qoff > pca_a[i].chain_qend) {
            int r = s_fill_one_gap(pca_a + i, pca_a + i + 1, all_pca_a, all_pca_c);
            if (r >= 0) {
                new_pca_list.push_back(all_pca_a[r]);
                find_gapped_fill_pca = 1;
            }
        }
    }
    new_pca_list.push_back(pca_a[pca_c - 1]);
    if (pca_a[pca_c-1].qsize - pca_a[pca_c-1].chain_qend > 50) {
        int r = s_fill_right_most_gap(pca_a + pca_c - 1, all_pca_a, all_pca_c);
        if (r >= 0) {
            new_pca_list.push_back(all_pca_a[r]);
            find_gapped_fill_pca = 1;
        }
    }
    if (!find_gapped_fill_pca) return;

    if (chain_info) {
        chain_info->is_complete_chain = s_is_complete_chain(vdfa, vdfc, new_pca_list.data(), new_pca_list.size());
        chain_info->is_perfect_chain = s_is_perfect_chain(new_pca_list.data(), new_pca_list.size());
        chain_info->is_pseudo_perfect_chain = s_is_pseudo_perfect_chain(new_pca_list.data(), new_pca_list.size());
        chain_info->q_cov = compute_pca_list_q_cov(new_pca_list.data(), new_pca_list.size());
    }
    pca_list.clear();
    pca_list.insert(pca_list.end(), new_pca_list.begin(), new_pca_list.end());
}

static bool 
s_is_validate_pca_offset_relations(PoreCAlign* al, PoreCAlign* ar, const EChainType type)
{
    if (type == ePerfectChain) {
        return s_is_perfect_chain_offset(al, ar, PC_MAX_ALIGN_OVLP);
    } else if (type == ePseudoPerfectChain) {
        return s_is_pseudo_perfect_chain_offset(al, ar, PC_MAX_ALIGN_OVLP);
    }
    return is_validate_chain_offset_general(al, ar, PC_MAX_ALIGN_OVLP);
}

static bool
s_is_correct_chain_type(const int* vdfa, const int vdfc, PoreCAlign* pca_a, int pca_c, const EChainType type)
{
    if (type == ePerfectChain) {
        if (!s_is_complete_chain(vdfa, vdfc, pca_a, pca_c)) return false;
        if (!s_is_perfect_chain(pca_a, pca_c)) return false;
        return true;
    } else if (type == ePseudoPerfectChain) {
        if (!s_is_complete_chain(vdfa, vdfc, pca_a, pca_c)) return false;
        if (!s_is_pseudo_perfect_chain(pca_a, pca_c)) return false;
        return true;
    } else if (type == eCompleteChain) {
        if (!s_is_complete_chain(vdfa, vdfc, pca_a, pca_c)) return false;
        return true;
    }
    return true;
}

static void 
pca_supports_setup(pca_supports* support, PoreCAlign* pca_a, int pca_c, const EChainType chain_type)
{
    support->supports.resize(pca_c);
    for (int i = 0; i < pca_c; ++i) {
        pca_support* p = &support->supports[i];
        p->has_precessor = FALSE;
        p->has_successor = FALSE;
        p->succ_offset = 0;
        p->succ_cnt = 0;
        p->prec_offset = 0;
        p->prec_cnt = 0;
    }

    support->succ_idx_list.clear();
    for (int i = 0; i < pca_c; ++i) support->succ_idx_list.push_back(i);

    for (int i = 0; i < pca_c; ++i) {
        pca_support* pi = &support->supports[i];
        pi->succ_offset = support->succ_idx_list.size();
        for  (int j = i + 1; j < pca_c; ++j) {
            pca_support* pj = &support->supports[j];
            if (s_is_validate_pca_offset_relations(pca_a + i, pca_a + j, chain_type)) {
                support->succ_idx_list.push_back(j);
                ++pi->succ_cnt;
                pi->has_successor = TRUE;
            }
        }
    }

    for (int i = 0; i < pca_c; ++i) {
        pca_support* pi = &support->supports[i];
        pi->prec_offset = support->succ_idx_list.size();
        for (int j = i - 1; j >= 0; --j) {
            pca_support* pj = &support->supports[j];
            if (s_is_validate_pca_offset_relations(pca_a + j, pca_a + i, chain_type)) {
                support->succ_idx_list.push_back(j);
                ++pi->prec_cnt;
                pi->has_precessor = TRUE;
            }
        }
    }
}

struct FragChainPoint
{
    int pca_idx;
    int score;
    int qcov;
    int n_frag;
    int prev;
};

void
s_choose_a_solution(PoreCAlignChainData* pca_chain_data,
    const int* vdfa, 
    const int vdfc, 
    PoreCAlign* pca_a, 
    int pca_c, 
    const EChainType chain_type, 
    std::vector<PoreCAlign>& _chain, 
    int* cov)
{
    if (pca_c == 0) return;

    vector<pca_chain_dp_point> dp_p_list(pca_c);
    pca_chain_dp_point* dp_p_a = dp_p_list.data();

    /// try to find a complete and perfect chain
    for (int i = 0; i < pca_c; ++i) {
        pca_chain_dp_point* p = dp_p_a + i;
        p->pca_a_idx = i;
        p->parent = NULL;
        p->q_cov = pca_a[i].chain_qend - pca_a[i].chain_qoff;
        p->s_dist = 0;
    }
    for (int i = pca_c - 1; i >= 0; --i) {
        pca_chain_dp_point* dp_i = dp_p_a + i;
        PoreCAlign* pi = pca_a + i;
        for (int j = i - 1; j >= 0; -- j) {
            pca_chain_dp_point* dp_j = dp_p_a + j;
            PoreCAlign* pj = pca_a + j;
            if (!s_is_validate_pca_offset_relations(pj, pi, chain_type)) continue;
            int q_cov = pj->qc;
            if (pj->chain_qend > pi->chain_qoff) q_cov -= (pj->chain_qend - pi->chain_qoff);
            q_cov += dp_i->q_cov;
            if (q_cov > dp_j->q_cov) {
                dp_j->q_cov = q_cov;
                dp_j->parent = dp_i;
            }
        }
    }
    int max_q_cov = 0;
    int max_j = - 1;
    for (int i = 0; i < pca_c; ++i) {
        if (dp_p_a[i].q_cov > max_q_cov) {
            max_q_cov = dp_p_a[i].q_cov;
            max_j = i;
        }
    }
    vector<PoreCAlign> chain;
    pca_chain_dp_point* p = dp_p_a + max_j;
    PoreCAlignChain _best_chain, *best_chain = &_best_chain;
    while (p) {
        int pca_a_idx = p->pca_a_idx;
        chain.push_back(pca_a[pca_a_idx]);
        best_chain->pca_a_idx_list.push_back(pca_a_idx);
        best_chain->chain.push_back(pca_a[pca_a_idx]);
        p = (pca_chain_dp_point*)p->parent;
    }
    best_chain->is_valid = 1;
    best_chain->sid = -1;
    best_chain->q_cov = max_q_cov;
    best_chain->avg_pi = compute_pca_list_avg_pi(best_chain->chain.data(), best_chain->chain.size());

    //HBN_LOG("max_cov = %d, chain_type = %d", max_q_cov, chain_type);
    //for (auto& xpca : chain) dump_chain_pca(fprintf, stderr, xpca, -1);
    if (!s_is_correct_chain_type(vdfa, vdfc, chain.data(), chain.size(), chain_type)) return;

    const int kMaxQCov = compute_pca_list_q_cov(chain.data(), chain.size());
    const int query_size = pca_a[0].qsize;
    pca_supports_setup(&pca_chain_data->supports, pca_a, pca_c, chain_type);

    SmallObjectAlloc* soa = pca_chain_data->dp_point_soa;
    SmallObjectAllocClear(soa);
    vector<pca_chain_dp_point*> dp_stack;
    int n_chain = 0;

    const pca_chain_dp_point* dp_array = dp_p_list.data();
    for (int i = 0; i < pca_c; ++i) {
        //if (kv_A(pca_chain_data->supports.supports, i).has_precessor) continue;
        //if (pca_a[i].chain_qoff > kMaxGap) continue;
        if (dp_array[i].q_cov != kMaxQCov) continue;
        pca_chain_dp_point* point = (pca_chain_dp_point*)SmallObjectAllocAlloc(soa, 1);
        point->parent = NULL;
        point->pca_a_idx = i;
        point->q_cov = pca_a[i].qc;
        point->s_dist = 0;
        dp_stack.push_back(point);
    }

    PoreCAlignChain _dp_chain, *dp_chain = &_dp_chain;
    int found_chain = 0;
    while (!dp_stack.empty() && n_chain < 50) {
        pca_chain_dp_point* curr = dp_stack.back();
        dp_stack.pop_back();
        int n_succ = pca_chain_data->supports.supports[curr->pca_a_idx].succ_cnt;

        if (n_succ == 0) {
            dp_chain->clear();
            pca_chain_dp_point* tail = curr;
            while (tail) {
                int pca_a_idx = tail->pca_a_idx;
                dp_chain->pca_a_idx_list.push_back(pca_a_idx);
                dp_chain->chain.push_back(pca_a[pca_a_idx]);
                tail = (pca_chain_dp_point*)tail->parent;
            }
            reverse(dp_chain->pca_a_idx_list.begin(), dp_chain->pca_a_idx_list.end());
            reverse(dp_chain->chain.begin(), dp_chain->chain.end());
            int q_cov = compute_pca_list_q_cov(dp_chain->chain.data(), dp_chain->chain.size());
            hbn_assert(q_cov <= query_size);
            //if (q_cov + kMaxMissedCovBase < query_size) {
            if (q_cov < kMaxQCov) {
                memset(curr, 0, sizeof(pca_chain_dp_point));
                SmallocObjactAllocDeallocOne(soa, curr);
                continue;
            }
#if 0
            HBN_LOG("find one chain %d, q_cov = %d, target chain type = %d", n_chain, q_cov, chain_type);
            for (size_t i = 0; i < dp_chain->chain.size(); ++i) {
                PoreCAlign* pca = &dp_chain->chain[i];
                dump_chain_pca(fprintf, stderr, *pca, i);
            }
#endif      
	    ++found_chain;
            dp_chain->is_valid = 1;
            dp_chain->q_cov = q_cov;
            dp_chain->avg_pi = compute_pca_list_avg_pi(dp_chain->chain.data(), dp_chain->chain.size());
            
            if (s_is_correct_chain_type(vdfa, vdfc, dp_chain->chain.data(), dp_chain->chain.size(), chain_type)) {
                if (best_chain->chain.empty()) {
                    swap(best_chain, dp_chain);
                } else {
                    if (s_pca_list_relation(best_chain->pca_a_idx_list, dp_chain->pca_a_idx_list) == 1) {
                        swap(best_chain, dp_chain);
                    } else if (best_chain->avg_pi < dp_chain->avg_pi) {// || best_chain->chain.size() > dp_chain->chain.size()) {
                        swap(best_chain, dp_chain);
                    }
                }
                ++n_chain;
            }
            memset(curr, 0, sizeof(pca_chain_dp_point));
            SmallocObjactAllocDeallocOne(soa, curr);
	        if (found_chain == 10) break;
            continue;
        }

        PoreCAlign* al = pca_a + curr->pca_a_idx;
        int* pca_a_idx_list = pca_chain_data->supports.succ_idx_list.data()
                              +
                              pca_chain_data->supports.supports[curr->pca_a_idx].succ_offset;
        for (int i = 0; i < n_succ; ++i) {
            int pca_a_idx = pca_a_idx_list[i];
            PoreCAlign* ar = pca_a + pca_a_idx;
            int max_q_cov = curr->q_cov + dp_array[pca_a_idx].q_cov;
	    if (al->chain_qend > ar->chain_qoff) max_q_cov -= (al->chain_qend - ar->chain_qoff);
            if (max_q_cov < kMaxQCov) continue;
            if (!s_is_validate_pca_offset_relations(al, ar, chain_type)) continue;
            int q_cov = ar->chain_qend - ar->chain_qoff;
            if (al->chain_qend > ar->chain_qoff) q_cov -= (al->chain_qend - ar->chain_qoff);
            hbn_assert(q_cov >= 0);
            q_cov += curr->q_cov;
            hbn_assert(q_cov <= query_size);
            pca_chain_dp_point* son = (pca_chain_dp_point*)SmallObjectAllocAllocOne(soa);
            son->q_cov = q_cov;
            son->parent = curr;
            son->pca_a_idx = pca_a_idx;
            dp_stack.push_back(son);
        }
    }
    if (best_chain->chain.empty()) return;
    chain.assign(best_chain->chain.begin(), best_chain->chain.end());
    *cov = compute_pca_list_q_cov(chain.data(), chain.size());

    double pi1 = chain.empty() ? 0.0 : compute_pca_list_avg_pi(chain.data(), chain.size());
    double pi2 = _chain.empty() ? 0.0 : compute_pca_list_avg_pi(_chain.data(), _chain.size());
    if (pi1 > pi2 + 10.0) _chain.assign(chain.begin(), chain.end());
}

static bool 
x_select_pca_chain_by_cov(PoreCAlignChainData* pca_chain_data, PoreCAlign* pca_a, int pca_c, const int* vdfa, const int vdfc, EChainType chain_type, vector<PoreCAlign>& chain)
{
    vector<FragChainPoint> fcp_list(pca_c);
    FragChainPoint* fcp_a = fcp_list.data();
    for (int i = 0; i < pca_c; ++i) {
        FragChainPoint* p = fcp_a + i;
        p->pca_idx = i;
        p->qcov = pca_a[i].qc;
        p->score = pca_a[i].score;
        p->n_frag = 1;
        p->prev = -1;
    }
    for (int i = 0; i < pca_c; ++i) {
        int max_cov = pca_a[i].qc;
        int n_frag = 0;
        int prev = -1;
        PoreCAlign* pi = pca_a + i;
        for (int j = i - 1; j >= 0; --j) {
            PoreCAlign* pj = pca_a + j;
            if (!s_is_validate_pca_offset_relations(pj, pi, chain_type)) continue;
            int cov = pi->qc + fcp_a[j].qcov;
            if (pj->chain_qend > pi->chain_qoff) cov -= (pj->chain_qend - pi->chain_qoff);
            if (cov > max_cov) {
                max_cov = cov;
                n_frag = fcp_a[j].n_frag;
                prev = j;
            } else if (max_cov - cov <= 30 && fcp_a[j].n_frag < n_frag) {
                max_cov = cov;
                n_frag = fcp_a[j].n_frag;
                prev = j;
            }
        }
        if (prev >= 0) {
            fcp_a[i].n_frag = n_frag + 1;
            fcp_a[i].prev = prev;
            fcp_a[i].qcov = max_cov;
        }
    }

    int max_i = -1, max_cov = 0, n_frag = pca_c + 1;
    for (int i = 0; i < pca_c; ++i) {
        if (fcp_a[i].qcov > max_cov) {
            max_i = i;
            max_cov = fcp_a[i].qcov;
            n_frag = fcp_a[i].n_frag;
        } else if (fcp_a[i].qcov == max_cov && fcp_a[i].n_frag < n_frag) {
            max_i = i;
            n_frag = fcp_a[i].n_frag;
        }
    }

    chain.clear();
    int p = max_i;
    while (p >= 0) {
        chain.push_back(pca_a[p]);
        p = fcp_a[p].prev;
    }
    reverse(chain.begin(), chain.end());
    if (!s_is_correct_chain_type(vdfa, vdfc, chain.data(), chain.size(), chain_type)) return false;

    int cov = 0;
    s_choose_a_solution(pca_chain_data, vdfa, vdfc, pca_a, pca_c, chain_type, chain, &cov);
    return true;
}

static bool 
x_select_pca_chain_greedy(PoreCAlign* pca_a, int pca_c, vector<PoreCAlign>& chain)
{
    chain.clear();
    pdqsort(pca_a, pca_a + pca_c, [](const PoreCAlign& x, const PoreCAlign& y) { return x.qc > y.qc; });
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pi = pca_a + i;
        bool has_overlap = false;
        for (auto& xp : chain) {
            if (pi->soff <= xp.soff && pi->send > xp.soff) { has_overlap = true; break; }
            if (xp.soff <= pi->soff && xp.send > pi->soff) { has_overlap = true; break; }
            int ovlp = 0;
            if (pi->chain_qoff <= xp.chain_qoff && pi->chain_qend > xp.chain_qoff) {
                ovlp = pi->chain_qend - xp.chain_qoff;
            } else if (xp.chain_qoff <= pi->chain_qoff && xp.chain_qend > pi->chain_qoff) {
                ovlp = xp.chain_qend - pi->chain_qoff;
            }
            if (ovlp > pi->qc * 0.4) { has_overlap = true; break; }
        }
        if (!has_overlap) chain.push_back(*pi);
    }
    pdqsort(chain.begin(), chain.end(), [](const PoreCAlign& x, const PoreCAlign& y) { return x.chain_qoff < y.chain_qoff; });
    return true;
}

bool 
select_pca_chain(PoreCAlignChainData* pca_chain_data,
    const int* vdfa, 
    const int vdfc, 
    PoreCAlign* pca_array, 
    int pca_count, 
    const EChainType chain_type, 
    std::vector<PoreCAlign>& chain, 
    int* cov)
{
    if (pca_count == 0) return 0;
    pdqsort(pca_array, pca_array + pca_count, [](const PoreCAlign& a, const PoreCAlign& b)->bool { return a.chain_qoff < b.chain_qoff; });

#if 0
    HBN_LOG("org pca list");
    for (int i = 0; i < pca_count; ++i) {
        dump_chain_pca(fprintf, stderr, pca_array[i], i);
        //dump_pca(fprintf, stderr, pca_a[i], i);
    }
#endif

    vector<PoreCAlign> pca_list;
    pca_list.insert(pca_list.end(), pca_array, pca_array + pca_count);
    PoreCAlign* pca_a = pca_list.data();
    int pca_c = pca_list.size();
    if (chain_type == eMaxCovChain) {
        for (int i = 0; i < pca_c; ++i) {
	    //if (pca_a[i].qend - pca_a[i].qoff >= 500) continue;
            int ld = abs(pca_a[i].soff - pca_a[i].enzyme_soff) <= 20;
            int rd = abs(pca_a[i].send - pca_a[i].enzyme_send) <= 20;
            if (ld || rd) continue;
            pca_a[i].qid = -1;
        }
    }
    {
        int n = 0;
        for (int i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
        pca_c = n;
    }

#if 1
    pdqsort(pca_a, pca_a + pca_c, [](const PoreCAlign& x, const PoreCAlign& y) { return x.score > y.score; });
	const int E1 = 40, E2 = 40;
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pi = pca_a + i;
        if (pi->score < 50) { pi->qid = -1; continue; }
	if (pi->chain_qend - pi->chain_qoff < 50) { pi->qid = -1; continue; }
        if (pi->qid == -1) continue;
        for (int j = i + 1; j < pca_c; ++j) {
            PoreCAlign* pj = pca_a + j;
            if (pj->qid == -1) continue;
            bool r = abs(pi->chain_qoff - pj->chain_qoff) <= E1 && pj->chain_qend <= pi->chain_qend + E2;
            if (r) {
                pj->qid = -1;
                continue;
            }
            r = abs(pi->chain_qend - pj->chain_qend) <= E1 && pj->chain_qoff + E2 >= pi->chain_qoff;
            if (r) {
                pj->qid = -1;
                continue;
            }
        }
    }
    {
        int n = 0;
        for (int i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
        pca_list.resize(n);
    }
    if (pca_list.empty()) return 0;
    pca_a = pca_list.data();
    pca_c = pca_list.size();
    pdqsort(pca_a, pca_a + pca_c, [](const PoreCAlign& a, const PoreCAlign& b)->bool { return a.chain_qoff < b.chain_qoff; });
#endif
    for (int i = 0; i < pca_c; ++i) {
        pca_a[i].qc = pca_a[i].chain_qend - pca_a[i].chain_qoff;
    }

#if 0
    HBN_LOG("chaining pca list, type = %d", chain_type);
    for (int i = 0; i < pca_c; ++i) {
        dump_chain_pca(fprintf, stderr, pca_a[i], i);
        //dump_pca(fprintf, stderr, pca_a[i], i);
    }
#endif

    if (chain_type == eMaxCovChain) return x_select_pca_chain_greedy(pca_a, pca_c, chain);

    return x_select_pca_chain_by_cov(pca_chain_data, pca_a, pca_c, vdfa, vdfc, chain_type, chain);
}