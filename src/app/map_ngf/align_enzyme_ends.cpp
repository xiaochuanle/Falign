#include "align_enzyme_ends.hpp"

#include <algorithm>
#include <vector>
#include <cmath>

#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/edlib.h"

using namespace std;

static const int kEnzymeNbhd = 50;
static const int kMaxW = 50;
static const double kEpsilon = 0.2;

static void
s_add_enzyme_pos_left(const int* reloci_array,
    const int intv_cnt,
    const int intv_i,
    const int start_offset,
    const int end_offset,
    const int max_dist,
    vector<int>& offset_list)
{
    offset_list.clear();
    hbn_assert(start_offset >= reloci_array[intv_i] && start_offset <= reloci_array[intv_i + 1]);
    int i = intv_i;
    while (i >= 0) {
        if (start_offset - reloci_array[i] > max_dist) break;
        offset_list.push_back(reloci_array[i]);
        --i;
    }
    i = intv_i;
    while (i < intv_cnt) {
        if (reloci_array[i+1] - start_offset > max_dist) break;
        if (reloci_array[i+1] >= end_offset) break;
        offset_list.push_back(reloci_array[i+1]);
        ++i;
    }
    sort(offset_list.begin(), offset_list.end());
}

static void
s_add_enzyme_pos_right(const int* reloci_array,
    const int intv_cnt,
    const int intv_i,
    const int start_offset,
    const int end_offset,
    const int max_dist,
    vector<int>& offset_list)
{
    offset_list.clear();
    hbn_assert(end_offset >= reloci_array[intv_i] && end_offset <= reloci_array[intv_i + 1]);
    int i = intv_i;
    while (i >= 0) {
        if (reloci_array[i] <= start_offset) break;
        if (end_offset - reloci_array[i] > max_dist) break;
        offset_list.push_back(reloci_array[i]);
        --i;
    }    
    i = intv_i;
    while (i < intv_cnt) {
        if (reloci_array[i+1] - end_offset > max_dist) break;
        offset_list.push_back(reloci_array[i+1]);
        ++i;
    }
    sort(offset_list.begin(), offset_list.end());
}

static bool
s_realign_subseq_1(HbnTracebackData* tbck_data,
    int qdir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const u8* subject,
    int sid,
    int* _qb,
    int* _qe,
    int* _sb,
    int* _se,
    double* _pi)
{
    int qb = *_qb;
    int qe = *_qe;
    int sb = *_sb;
    int se = *_se;
    const u8* query = (qdir == FWD) ? fwd_query : rev_query;
    int q_l = qe - qb;
    int s_l = se - sb;
    int tolerance = hbn_max(q_l, s_l) * 0.36;
    int dist = small_edlib_nw_1(query + qb, q_l,
                subject + sb, s_l,
                tolerance,
                tbck_data->small_edlib);
    if (dist == -1) return false;


    int x_sb = 0;
    int x_se = s_l;
    hbn_assert(x_sb == 0);
    hbn_assert(x_se == s_l);
    *_qb = qb;
    *_qe = qe;
    *_sb = sb + x_sb;
    *_se = sb + x_se;
    *_pi = 100.0 * (1.0 - 2.0 * dist / (q_l + s_l));
    return true;
}

static bool
s_realign_subseq(HbnTracebackData* tbck_data,
    int qdir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const u8* subject,
    int sid,
    int* _qb,
    int* _qe,
    int* _sb,
    int* _se,
    double* _pi)
{
    int qb = *_qb;
    int qe = *_qe;
    int sb = *_sb;
    int se = *_se;
    if (qe - qb < SMALL_EDLIB_MAX_SEQ_SIZE
        &&
        se - sb < SMALL_EDLIB_MAX_SEQ_SIZE) {
        return s_realign_subseq_1(tbck_data, qdir, fwd_query, rev_query,
            query_size, subject, sid, _qb, _qe, _sb, _se, _pi);
    }
    //return false;

    const u8* query = (qdir == FWD) ? fwd_query : rev_query;
    int q_l = qe - qb;
    ks_set_size(&tbck_data->ext_qabuf, q_l);
    for (int i = qb, p = 0; i < qe; ++i, ++p) {
        ks_A(tbck_data->ext_qabuf, p) = DECODE_RESIDUE(query[i]);
    }
    int s_l = se - sb;
    ks_set_size(&tbck_data->ext_sabuf, s_l);
    for (int i = sb, p = 0; i < se; ++i, ++p) {
        ks_A(tbck_data->ext_sabuf, p) = DECODE_RESIDUE(subject[i]);
    }
    EdlibAlignTask task = EDLIB_TASK_LOC;
    int tolerance = hbn_max(q_l, s_l) * 0.36;
    EdlibAlignResult align = edlibAlign(ks_s(tbck_data->ext_qabuf),
                                q_l,
                                ks_s(tbck_data->ext_sabuf),
                                s_l,
                                edlibNewAlignConfig(tolerance, EDLIB_MODE_NW, task, NULL, 0));
    if (align.numLocations == 0) {
        edlibFreeAlignResult(align);
        return false;
    }

    int x_sb = align.startLocations[0];
    int x_se = align.endLocations[0] + 1;
    hbn_assert(x_sb == 0);
    hbn_assert(x_se == s_l);
    *_qb = qb;
    *_qe = qe;
    *_sb = sb + x_sb;
    *_se = sb + x_se;
    *_pi = 100.0 * (1.0 - 2.0 * align.editDistance / (q_l + s_l));
    edlibFreeAlignResult(align);
    return true;
}

static int
s_align_right_enzyme_ends(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const int qid,
    const u8* fwd_query,
    const u8* rev_query,
    int qdir,
    int qsize,
    int sid,
    const u8* fwd_subject,
    int ssize,
    int qb,
    int qe,
    int sb,
    int se,
    double perc_identity,
    int* cov_stats,
    vector<PoreCAlign>& pca_list)
{
//fprintf(stderr, "find right enzyme pos for [%d, %d] x [%d, %d]\n", qb, qe, sb, se);
    const int* q_reloci_array = (qdir == FWD)
                                ?
                                kv_data(qvep_list->fwd_vdf_endpoint_list)
                                :
                                kv_data(qvep_list->rev_vdf_endpoint_list);
    const int q_reloci_cnt = kv_size(qvep_list->rev_vdf_endpoint_list);
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[sid].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[sid].enzyme_loci_cnt;

    int r_q_intv_c = 0;
    int r_q_intv_i = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, qe, &r_q_intv_c);  
    int r_s_intv_c = 0;
    int r_s_intv_i = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, se, &r_s_intv_c);
    vector<int> qr_list;
    vector<int> sr_list;
    s_add_enzyme_pos_right(q_reloci_array, r_q_intv_c, r_q_intv_i, qb, qe, kEnzymeNbhd, qr_list);
    s_add_enzyme_pos_right(s_reloci_array, r_s_intv_c, r_s_intv_i, sb, se, kEnzymeNbhd, sr_list);
   // HBN_LOG("query enzyme sites:"); for (auto qi : qr_list) fprintf(stderr, "%d\t", qi); fprintf(stderr, "\n");
   // HBN_LOG("subject enzyme sites:"); for (auto si : sr_list) fprintf(stderr, "%d\t", si); fprintf(stderr, "\n");
    int added_pca = 0;
    for (const auto& qend : qr_list) {
        if (qend <= qb || qend > qsize) continue;
//fprintf(stderr, "process enzyme qend %d\n", qend);
        int best_send = -1;
        double best_ddf = 2.0;
        for (const auto& send : sr_list) {
            if (send <= sb || send > ssize) continue;
//fprintf(stderr, "process enzyme send %d\n", send);
            hbn_assert(qb < qend);
            hbn_assert(sb < send);
            int q_d = qend - qb;
            int s_d = send - sb;
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            if (ddf < best_ddf) {
                best_ddf = ddf;
                best_send = send;
            }
        }
        if (best_send == -1) continue;
#if 1
        int x_qb = qb, x_qe = qend, x_sb = sb, x_se = best_send;
        double pi = perc_identity;
            ++added_pca;
            PoreCAlign pca;
            pca.qid = qid;
            pca.qdir = qdir;
            pca.qoff = x_qb;
            pca.qend = x_qe;
            pca.qsize = qsize;
            pca.sid = sid;
            pca.soff = x_sb;
            pca.send = x_se;
            pca.ssize = ssize;
            pca.pi = pi;
            pca_list.push_back(pca);
#else
        int x_qb = qb, x_qe = qend, x_sb = sb, x_se = best_send;
        double pi = .0;
        if (s_realign_subseq(tbck_data, qdir, fwd_query, rev_query, qsize,
            fwd_subject, sid, &x_qb, &x_qe, &x_sb, &x_se, &pi)) {
            //HBN_LOG("find perfect pca");
            //fprintf(stderr, "[%d, %d] x [%d, %d], %g\n", x_qb, x_qe, x_sb, x_se, pi);
            ++added_pca;
            PoreCAlign pca;
            pca.qid = qid;
            pca.qdir = qdir;
            pca.qoff = x_qb;
            pca.qend = x_qe;
            pca.qsize = qsize;
            pca.sid = sid;
            pca.soff = x_sb;
            pca.send = x_se;
            pca.ssize = ssize;
            pca.pi = pi;
            pca_list.push_back(pca);
        }
#endif
    }    
    return added_pca;
}

struct EnzymeStartOffsetPair
{
    int qoff, soff;
    double pi;
};

static void
s_find_enzyme_pos(const int* reloci_array,
    const int reloci_cnt,
    const int start,
    const int end,
    int* enzyme_start,
    int* enzyme_end)
{
    int s_intv_c = 0;
    int s_intv_i = offset_to_enzyme_intv_idx(reloci_array, reloci_cnt, start, &s_intv_c);     
    int e_intv_c = 0;
    int e_intv_i = offset_to_enzyme_intv_idx(reloci_array, reloci_cnt, end, &e_intv_c);
    hbn_assert(s_intv_c == e_intv_c);
    hbn_assert(s_intv_i <= e_intv_i);

    if (s_intv_i == e_intv_i) {
        *enzyme_start = reloci_array[s_intv_i];
        *enzyme_end = reloci_array[e_intv_i + 1];
        return;
    } 
    
    if (s_intv_i + 1 == e_intv_i) {
        int s_l_d = start - reloci_array[s_intv_i];
        int s_r_d = reloci_array[s_intv_i+1] - start;
        hbn_assert(s_l_d >= 0);
        hbn_assert(s_r_d >= 0);
        int e_s = (s_l_d < s_r_d) ? reloci_array[s_intv_i] : reloci_array[s_intv_i+1];

        int e_l_d = end - reloci_array[e_intv_i];
        int e_r_d = reloci_array[e_intv_i+1] - end;
        hbn_assert(e_l_d >= 0);
        hbn_assert(e_r_d >= 0);
        int e_e = (e_l_d < e_r_d) ? reloci_array[e_intv_i] : reloci_array[e_intv_i+1];

        if (e_s == e_e) {
            if (s_r_d < e_l_d) {
                e_e = reloci_array[e_intv_i+1];
            } else {
                e_s = reloci_array[s_intv_i];
            }
        }
        *enzyme_start = e_s;
        *enzyme_end = e_e;
        return;
    }

    {
        int s_l_d = start - reloci_array[s_intv_i];
        int s_r_d = reloci_array[s_intv_i+1] - start;
        hbn_assert(s_l_d >= 0);
        hbn_assert(s_r_d >= 0);
        int e_s = (s_l_d < s_r_d) ? reloci_array[s_intv_i] : reloci_array[s_intv_i+1];

        int e_l_d = end - reloci_array[e_intv_i];
        int e_r_d = reloci_array[e_intv_i+1] - end;
        hbn_assert(e_l_d >= 0);
        hbn_assert(e_r_d >= 0);
        int e_e = (e_l_d < e_r_d) ? reloci_array[e_intv_i] : reloci_array[e_intv_i+1];

        hbn_assert(e_s < e_e);
        *enzyme_start = e_s;
        *enzyme_end = e_e;        
    }
}

static void
s_set_enzyme_pos_for_one_pca(RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    PoreCAlign* pca,
    int* cov_stats)
{
    fprintf(stderr, "set enzyme pos for [%d, %d] x [%d, %d], %g\n", pca->qoff, pca->qend, pca->soff, pca->send, pca->pi);
    const int* q_reloci_array = (pca->qdir == FWD)
                                ?
                                kv_data(qvep_list->fwd_vdf_endpoint_list)
                                :
                                kv_data(qvep_list->rev_vdf_endpoint_list);
    const int q_reloci_cnt = kv_size(qvep_list->rev_vdf_endpoint_list);
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_cnt;
    int enzyme_qb = 0;
    int enzyme_qe = 0;
    int enzyme_sb = 0;
    int enzyme_se = 0;

    s_find_enzyme_pos(q_reloci_array, q_reloci_cnt, pca->qoff, pca->qend, &enzyme_qb, &enzyme_qe);
    s_find_enzyme_pos(s_reloci_array, s_reloci_cnt, pca->soff, pca->send, &enzyme_sb, &enzyme_se);

    fprintf(stderr, "eqb = %d, eqe = %d, esb = %d, ese = %d\n", enzyme_qb, enzyme_qe, enzyme_sb, enzyme_se);

    const int kMaxD = 10;
    int l_enzyme_match = (enzyme_qb == 0)
                         ||
                         (enzyme_sb == 0)
                         ||
                         (pca->qoff == enzyme_qb)
                         ||
                         abs(pca->soff - enzyme_sb) <= kMaxD;
    int r_enzyme_match = (enzyme_se == pca->qsize)
                         ||
                         (enzyme_se == pca->ssize)
                         ||
                         (pca->qend == enzyme_qe)
                         ||
                         abs(pca->send - enzyme_se) <= kMaxD;

    if (l_enzyme_match == 0 && r_enzyme_match == 0) {
        pca->qid = -1;
    } else {
        pca->enzyme_qoff = enzyme_qb;
        pca->enzyme_qend = enzyme_qe;
        pca->enzyme_soff = enzyme_sb;
        pca->enzyme_send = enzyme_se;
        set_pca_l_enzyme_match(*pca);
        set_pca_r_enzyme_match(*pca);
        dump_pca(fprintf, stderr, *pca, -1);
        HBN_LOG("l_match = %d, r_match = %d", pca->l_enzyme_match, pca->r_enzyme_match);
        for (int p = pca->qoff; p < pca->qend; ++p) ++cov_stats[p];
    }
}

int 
align_enzyme_ends(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const int qid,
    const u8* fwd_query,
    const u8* rev_query,
    int qdir,
    int qsize,
    int sid,
    const u8* fwd_subject,
    int ssize,
    int qb,
    int qe,
    int sb,
    int se,
    double perc_identity,
    int* cov_stats,
    vector<PoreCAlign>& pca_list)
{
    //fprintf(stderr, "===========> find enzyme ends for [%d, %d] x [%d, %d], %g\n", qb, qe, sb, se, perc_identity);
    const int* q_reloci_array = (qdir == FWD)
                                ?
                                kv_data(qvep_list->fwd_vdf_endpoint_list)
                                :
                                kv_data(qvep_list->rev_vdf_endpoint_list);
    const int q_reloci_cnt = kv_size(qvep_list->rev_vdf_endpoint_list);
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[sid].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[sid].enzyme_loci_cnt;

    int l_q_intv_c = 0;
    int l_q_intv_i = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, qb, &l_q_intv_c);  
    int l_s_intv_c = 0;
    int l_s_intv_i = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, sb, &l_s_intv_c);
    vector<int> ql_list;
    vector<int> sl_list;
    s_add_enzyme_pos_left(q_reloci_array, l_q_intv_c, l_q_intv_i, qb, qe, kEnzymeNbhd, ql_list);
    s_add_enzyme_pos_left(s_reloci_array, l_s_intv_c, l_s_intv_i, sb, se, kEnzymeNbhd, sl_list);

    vector<EnzymeStartOffsetPair> start_offset_pair_list;
    for (const auto& qoff : ql_list) {
        if (qoff < 0 || qoff >= qsize) continue;
       //fprintf(stderr, "process enzyme read pos %d\n", qoff);
        int best_soff = -1;
        double best_ddf = 2.0;
        for (const auto& soff : sl_list) {
            if (soff < 0 || soff >= ssize) continue;
//fprintf(stderr, "  process enzyme subject pos %d\n", soff);
            hbn_assert(qe > qoff);
            hbn_assert(se > soff);
            int q_d = qe - qoff;
            int s_d = se - soff;
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            if (ddf < best_ddf) {
                best_ddf = ddf;
                best_soff = soff;
            }
        }
        if (best_soff == -1) continue;
#if 0
        int x_qb = qoff, x_qe = qe, x_sb = best_soff, x_se = se;
        double pi = .0;
        if (s_realign_subseq(tbck_data, qdir, fwd_query, rev_query, qsize,
            fwd_subject, sid, &x_qb, &x_qe, &x_sb, &x_se, &pi)) {
            EnzymeStartOffsetPair p;
            p.qoff = x_qb;
            p.soff = x_sb;
            p.pi = pi;
            start_offset_pair_list.push_back(p);
            //fprintf(stderr, "find enzyme offset [%d, %d] x [%d, %d], %g\n", x_qb, x_qe, x_sb, x_se, pi);
        }
#else
        int x_qb = qoff, x_qe = qe, x_sb = best_soff, x_se = se;
        double pi = perc_identity;
        EnzymeStartOffsetPair p;
        p.qoff = x_qb;
        p.soff = x_sb;
        p.pi = pi;
        start_offset_pair_list.push_back(p);
        //fprintf(stderr, "find enzyme offset [%d, %d] x [%d, %d], %g\n", x_qb, x_qe, x_sb, x_se, pi);
#endif
    }

    if (start_offset_pair_list.empty()) {
            EnzymeStartOffsetPair p;
            p.qoff = qb;
            p.soff = sb;
            p.pi = perc_identity;
            start_offset_pair_list.push_back(p);
    }	

    int added_pca = 0;
    for (auto& p : start_offset_pair_list) {
        int n = 
        s_align_right_enzyme_ends(tbck_data, reloci_list, qvep_list,
            qid, fwd_query, rev_query, qdir, qsize,
            sid, fwd_subject, ssize, p.qoff, qe, p.soff, se, p.pi,
            cov_stats, pca_list);
        added_pca += n;
        if (!n) {
            PoreCAlign pca;
            pca.qid = qid;
            pca.qdir = qdir;
            pca.qoff = p.qoff;
            pca.qend = qe;
            pca.qsize = qsize;
            pca.sid = sid;
            pca.soff = p.soff;
            pca.send = se;
            pca.ssize = ssize;
            pca.pi = p.pi;
            pca_list.push_back(pca);
            ++added_pca;
        }
    }

    if (1) {
        PoreCAlign pca;
        pca.qid = qid;
        pca.qdir = qdir;
        pca.qoff = qb;
        pca.qend = qe;
        pca.qsize = qsize;
        pca.sid = sid;
        pca.soff = sb;
        pca.send = se;
        pca.ssize = ssize;
        pca.pi = perc_identity;
        pca_list.push_back(pca);
        ++added_pca;        
    }

    int n = pca_list.size();
    for (int i = 0; i < n; ++i) s_set_enzyme_pos_for_one_pca(reloci_list, qvep_list, &pca_list[i], cov_stats);
    int m = 0;
    for (int i = 0; i < n; ++i) if (pca_list[i].qid >= 0) pca_list[m++] = pca_list[i];
    pca_list.resize(m);
    return m;
}
