#include "fix_pca_offset_with_enzyme_info.hpp"

#include <algorithm>
#include <vector>
#include <cmath>

#include "../../algo/hbn_traceback_aux.h"
#include "../../algo/edlib.h"

using namespace std;

static const int kMaxW = 50;

/////////////////////

static void
s_find_broken_perfect_points(QueryVdfEndPointList* qvep_list, PoreCAlign* pca_a, int pca_c)
{
    const int max_offset = qvep_list->query_size + 1;
    kv_resize(u8, qvep_list->perfect_match_offset_point_list, max_offset);
    kv_resize(u8, qvep_list->perfect_match_end_point_list, max_offset);
    u8* perfect_match_offset_point_list = kv_data(qvep_list->perfect_match_offset_point_list);
    u8* perfect_match_end_point_list = kv_data(qvep_list->perfect_match_end_point_list);
    const int query_size = qvep_list->query_size;
    const int enzyme_size = qvep_list->enzyme->enzyme_size;

    fill(perfect_match_offset_point_list, perfect_match_offset_point_list + max_offset, 0);
    fill(perfect_match_end_point_list, perfect_match_end_point_list + max_offset, 0);
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        //dump_pca(fprintf, stderr, *pca, i);
        if (pca->enzyme_qoff > 0 && pca->enzyme_qoff == pca->qoff) {
            if (pca->qdir == FWD) {
                perfect_match_offset_point_list[pca->enzyme_qoff] = 1;
            } else {
                int offset = query_size - pca->enzyme_qoff;
                if (offset >= enzyme_size) offset -= enzyme_size;
                perfect_match_end_point_list[offset] = 1;
            }
        }
        if (pca->enzyme_qend != query_size && pca->enzyme_qend == pca->qend) {
            if (pca->qdir == FWD) {
                perfect_match_end_point_list[pca->enzyme_qend] = 1;
            } else {
                int offset = query_size - pca->enzyme_qend;
                if (offset >= enzyme_size) offset -= enzyme_size;
                perfect_match_offset_point_list[offset] = 1;
            }
        }
    }
}

static bool 
s_is_ddfs_supported(int qb, int qe, int sb, int se)
{
    const int kMaxW = 50;
    const double kEpsilon = 0.2;
    int q_d = qe - qb;
    int s_d = se - sb;
    int d_d = abs(q_d - s_d);
    //HBN_LOG("[%d, %d] x [%d, %d] q_d = %d, s_d = %d, d_d = %d", qb, qe, sb, se, q_d, s_d, d_d);
    if (d_d > kMaxW) return false;
    bool r = fabs(1.0 - 1.0 * q_d / s_d) <= kEpsilon;
    return r;
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
    int org_ql,
    int org_sl,
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

    //HBN_LOG("long edlib, org_q_l = %d, org_s_l = %d, q_l = %d, s_l = %d", 
    //    org_ql, org_sl, qe - qb, se - sb);
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
s_is_duplicate_query_offsets(const PoreCAlign* org_pca_a,
    const int org_pca_c,
    int qdir, int qb, int qe)
{
    for (int i = 0; i < org_pca_c; ++i) {
        const PoreCAlign* pca = org_pca_a + i;
        if (pca->qdir != qdir) continue;
        if (pca->qoff == qb && pca->qend == qe) return TRUE;
    }
    return FALSE;
}

void
s_fix_left_imperfect_pca(HbnTracebackData* tbck_data,
    const u8* fwd_query,
    const u8* rev_query,
    const u8* subject,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const double min_perc_identity,
    const PoreCAlign* org_pca_a,
    const int org_pca_c,
    PoreCAlign* pca,
    std::vector<PoreCAlign>& new_pca_list)
{
    //HBN_LOG("fix left imperfect pca");
    //dump_pca(fprintf, stderr, *pca, -1);
    const u8* perfect_match_offset_point_list = kv_data(qvep_list->perfect_match_offset_point_list);
    const u8* perfect_match_end_point_list = kv_data(qvep_list->perfect_match_end_point_list);
    const int* f_q_vdf_a = kv_data(qvep_list->fwd_vdf_endpoint_list);
    const int* r_q_vdf_a = kv_data(qvep_list->rev_vdf_endpoint_list);
    const int* q_vdf_a = (pca->qdir == FWD) ? f_q_vdf_a : r_q_vdf_a;
    const int q_vdf_c = kv_size(qvep_list->fwd_vdf_endpoint_list);
    const int* s_vdfi_a = reloci_list->reloci_array
                          +
                          reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_offset;
    const int s_vdfi_c = reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_cnt;
    const int query_size = pca->qsize;
    const int enzyme_size = reloci_list->enzyme.enzyme_size;

    int q_v_i = -1, q_v_c = 0;
    q_v_i = offset_to_enzyme_intv_idx(q_vdf_a, q_vdf_c, pca->qoff, &q_v_c);
    int s_v_i = -1, s_v_c = 0;
    s_v_i = offset_to_enzyme_intv_idx(s_vdfi_a, s_vdfi_c, pca->soff, &s_v_c);

    int i_from = 0;
    int i_to = q_v_i;
    while (i_to < q_v_c && q_vdf_a[i_to] < pca->qend) ++i_to;
    for (int i = i_from; i < i_to; ++i) {
        int x_qb = q_vdf_a[i];
if (abs(x_qb - pca->qoff) > 5000) continue;
        if (s_is_duplicate_query_offsets(org_pca_a, org_pca_c, pca->qdir, x_qb, pca->qend)) continue;
        if (s_is_duplicate_query_offsets(new_pca_list.data(), new_pca_list.size(),
            pca->qdir, x_qb, pca->qend)) continue;
        if (pca->qdir == FWD) {
            if (!perfect_match_end_point_list[x_qb]) continue;
        } else {
            int f_x_qb = query_size - x_qb;
            if (f_x_qb > 0 && f_x_qb != query_size) f_x_qb -= enzyme_size;
            hbn_assert(f_x_qb >= 0 && f_x_qb <= query_size);
            if (!perfect_match_offset_point_list[f_x_qb]) continue;
        }
        int j_from = s_v_i;
        while (j_from) {
            int x_sb = s_vdfi_a[j_from - 1];
            int q_d = pca->qend - x_qb;
            int s_d = pca->send - x_sb;
            hbn_assert(q_d > 0);
            hbn_assert(s_d > 0);
            if (s_d > q_d + kMaxW) break;
            --j_from;
        }
        int j_to = s_v_i;
        while (j_to < s_v_c && s_vdfi_a[j_to] < pca->send) ++j_to;

        for (int j = j_from; j < j_to; ++j) {
            int x_sb = s_vdfi_a[j];
if (abs(x_sb - pca->soff) > 5000) continue;
            if (!s_is_ddfs_supported(x_qb, pca->qend, x_sb, pca->send)) continue;
            //HBN_LOG("find feasible left end point %d, %d", x_qb, x_sb);
            int a_qb = x_qb, a_qe = pca->qend, a_sb = x_sb, a_se = pca->send;
            double a_pi = 0.0;
            if (!s_realign_subseq(tbck_data, pca->qdir, fwd_query, rev_query, query_size,
                subject, pca->sid, pca->qend - pca->qoff, pca->send - pca->soff, &a_qb, &a_qe, &a_sb, &a_se, &a_pi)) {
                continue;
            }
            if (a_pi < min_perc_identity) continue;

            PoreCAlign new_pca = *pca;
            new_pca.qoff = a_qb;
            new_pca.enzyme_qoff = x_qb;
            new_pca.soff = a_sb;
            new_pca.soff = a_sb;
            new_pca.enzyme_soff = x_sb;
            new_pca.pi = a_pi;
            if(set_pca_chain_offset(&new_pca, qvep_list->enzyme->enzyme_size)) {
                hbn_assert(new_pca.chain_qoff < new_pca.chain_qend);
                hbn_assert(new_pca.soff < new_pca.send);
                new_pca_list.push_back(new_pca);
                //HBN_LOG("find perfect pca");
                //dump_pca(fprintf, stderr, new_pca, -1);
                break;
            }
        }
    }
}

void
s_fix_right_imperfect_pca(HbnTracebackData* tbck_data,
    const u8* fwd_query,
    const u8* rev_query,
    const u8* subject,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const double min_perc_identity,
    const PoreCAlign* org_pca_a,
    const int org_pca_c,
    PoreCAlign* pca,
    std::vector<PoreCAlign>& new_pca_list)
{
    //HBN_LOG("fix right imperfect pca");
    //dump_pca(fprintf, stderr, *pca, -1);
    const u8* perfect_match_offset_point_list = kv_data(qvep_list->perfect_match_offset_point_list);
    const u8* perfect_match_end_point_list = kv_data(qvep_list->perfect_match_end_point_list);
    const int* f_q_vdf_a = kv_data(qvep_list->fwd_vdf_endpoint_list);
    const int* r_q_vdf_a = kv_data(qvep_list->rev_vdf_endpoint_list);
    const int* q_vdf_a = (pca->qdir == FWD) ? f_q_vdf_a : r_q_vdf_a;
    const int q_vdf_c = kv_size(qvep_list->fwd_vdf_endpoint_list);
    const int* s_vdfi_a = reloci_list->reloci_array
                          +
                          reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_offset;
    const int s_vdfi_c = reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_cnt;
    const int query_size = pca->qsize;
    const int enzyme_size = reloci_list->enzyme.enzyme_size;

    int q_v_i = -1, q_v_c = 0;
    q_v_i = offset_to_enzyme_intv_idx(q_vdf_a, q_vdf_c, pca->qend, &q_v_c);
    int s_v_i = -1, s_v_c = 0;
    s_v_i = offset_to_enzyme_intv_idx(s_vdfi_a, s_vdfi_c, pca->send, &s_v_c);

    int i_from = q_v_i;
    while (i_from && q_vdf_a[i_from-1] > pca->qoff) --i_from;
    int i_to = q_v_c;
    for (int i = i_from; i < i_to; ++i) {
        int x_qe = q_vdf_a[i];
if (abs(x_qe - pca->qend) > 5000) continue;
        if (x_qe <= pca->qoff) continue;
        if (s_is_duplicate_query_offsets(org_pca_a, org_pca_c, pca->qdir, pca->qoff, x_qe)) continue;
        if (s_is_duplicate_query_offsets(new_pca_list.data(), new_pca_list.size(),
            pca->qdir, pca->qoff, x_qe)) continue;
        if (pca->qdir == FWD) {
            if (!perfect_match_offset_point_list[x_qe]) continue;
        } else {
            int f_x_qe = query_size - x_qe;
            if (f_x_qe != 0 && f_x_qe != query_size) f_x_qe -= enzyme_size;
            if (!(f_x_qe >= 0 && f_x_qe <= query_size)) dump_pca(fprintf, stderr, *pca, -1);
            hbn_assert(f_x_qe >= 0 && f_x_qe <= query_size);
            if (!perfect_match_end_point_list[f_x_qe]) continue;
        }

        int j_from = s_v_i;
        while (j_from && s_vdfi_a[j_from-1] > pca->soff) --j_from;
        int j_to = s_v_i;
        while (j_to < s_v_c) {
            int q_d = x_qe - pca->qoff;
            hbn_assert(q_d > 0);
            int s_d = s_vdfi_a[j_to] - pca->soff;
            if (s_d > q_d + kMaxW) break;
            ++j_to;
        }

        for (int j = j_from; j < j_to; ++j) {
            int x_se = s_vdfi_a[j];
if (abs(x_se - pca->send) > 5000) continue;
            if (!s_is_ddfs_supported(pca->qoff, x_qe, pca->soff, x_se)) continue;
            //HBN_LOG("find feasible right end point %d, %d", x_qe, x_se);
            int a_qb = pca->qoff, a_qe = x_qe, a_sb = pca->soff, a_se = x_se;
            double a_pi = 0.0;
            if (!s_realign_subseq(tbck_data, pca->qdir, fwd_query, rev_query, query_size,
                subject, pca->sid, pca->qend - pca->qoff, pca->send - pca->soff, &a_qb, &a_qe, &a_sb, &a_se, &a_pi)) {
                continue;
            }
            if (a_pi < min_perc_identity) continue;

            PoreCAlign new_pca = *pca;
            new_pca.qend = a_qe;
            new_pca.enzyme_qend = x_qe;
            new_pca.send = a_se;
            new_pca.enzyme_send = x_se;
            new_pca.pi = a_pi;
            if(set_pca_chain_offset(&new_pca, qvep_list->enzyme->enzyme_size)) {
                hbn_assert(new_pca.chain_qoff < new_pca.chain_qend);
                hbn_assert(new_pca.soff < new_pca.send);
                new_pca_list.push_back(new_pca);
                //HBN_LOG("find perfect pca");
                //dump_pca(fprintf, stderr, new_pca, -1);
                break;
            }
        }
    }
}

static void
s_extend_pefect_right_end(HbnTracebackData* tbck_data,
    const u8* fwd_query,
    const u8* rev_query,
    const u8* subject,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const double min_perc_identity,
    const PoreCAlign* org_pca_a,
    const int org_pca_c,
    PoreCAlign* pca,
    vector<PoreCAlign>& added_pca_list)
{
    u8* perfect_match_offset_point_list = kv_data(qvep_list->perfect_match_offset_point_list);
    u8* perfect_match_end_point_list = kv_data(qvep_list->perfect_match_end_point_list);
    const int* s_vdfi_a = reloci_list->reloci_array
                          +
                          reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_offset;
    const int s_vdfi_c = reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_cnt;
    const int query_size = qvep_list->query_size;
    const int enzyme_size = qvep_list->enzyme->enzyme_size;

    if (pca->qdir == REV) {
        int offset = query_size - pca->enzyme_qend;
        if (offset >= enzyme_size) offset -= enzyme_size;
        hbn_assert(perfect_match_offset_point_list[offset]);
        if (perfect_match_end_point_list[offset]) return;
        int fwd_qoff = -1;
        for (int x = offset - 1; x > offset - 100 && x > 0 ; --x) {
            if (perfect_match_end_point_list[x]) {
                if (!perfect_match_offset_point_list[x]) {
                    fwd_qoff = x;
                    break;
                }
            }
        }
        if (fwd_qoff == -1) return;
        int rev_qoff = query_size - fwd_qoff;
        if (rev_qoff >= enzyme_size) rev_qoff -= enzyme_size;

        int s_v_i = -1, s_v_c = 0;
        s_v_i = offset_to_enzyme_intv_idx(s_vdfi_a, s_vdfi_c, pca->send, &s_v_c);
        int fwd_soff = -1;
        for (int x = s_v_i; x < s_v_c; ++x) {
            int soff = s_vdfi_a[x];
            if (pca->send >= soff) continue;
            if (soff - pca->send > 100) break;
            if (s_is_ddfs_supported(pca->qoff, rev_qoff, pca->soff, soff)) {
                fwd_soff = soff;
                break;
            }
        }
        if (fwd_soff == -1) return;
        //HBN_LOG("find target offset pair %d, %d", rev_qoff, fwd_soff);
        int a_qb = pca->qoff, a_qe = rev_qoff, a_sb = pca->soff, a_se = fwd_soff;
        double a_pi = 0.0;
        if (!s_realign_subseq(tbck_data, pca->qdir, fwd_query, rev_query, query_size,
            subject, pca->sid, pca->qend - pca->qoff, pca->send - pca->soff, &a_qb, &a_qe, &a_sb, &a_se, &a_pi)) {
            return;
        }
        if (a_pi < min_perc_identity) return;

        PoreCAlign new_pca = *pca;
        new_pca.qend = rev_qoff;
        new_pca.enzyme_qend = rev_qoff;
        new_pca.send = fwd_soff;
        new_pca.enzyme_send = fwd_soff;
        new_pca.pi = a_pi;
        if(set_pca_chain_offset(&new_pca, qvep_list->enzyme->enzyme_size)) {
            hbn_assert(new_pca.chain_qoff < new_pca.chain_qend);
            hbn_assert(new_pca.soff < new_pca.send);
            added_pca_list.push_back(new_pca);
            perfect_match_offset_point_list[fwd_qoff] = 1;
            //HBN_LOG("find extended pca");
            //dump_pca(fprintf, stderr, new_pca, -1);
        }
    } else {
        int offset = pca->enzyme_qend;
        hbn_assert(perfect_match_end_point_list[offset]);
        if (perfect_match_offset_point_list[offset]) return;
        int fwd_qoff = -1;
        for (int x = offset + 1; x < query_size - 1 && x < offset + 100; ++x) {
            if (perfect_match_offset_point_list[x]) {
                if (!perfect_match_end_point_list[x]) {
                    fwd_qoff = x;
                    break;
                }
            }
        }
        if (fwd_qoff == -1) return;

        int s_v_i = -1, s_v_c = 0;
        s_v_i = offset_to_enzyme_intv_idx(s_vdfi_a, s_vdfi_c, pca->send, &s_v_c);
        int fwd_soff = -1;
        for (int x = s_v_i; x < s_v_c; ++x) {
            int soff = s_vdfi_a[x];
            if (pca->send >= soff) continue;
            if (soff - pca->send > 100) break;
            if (s_is_ddfs_supported(pca->qoff, fwd_qoff, pca->soff, soff)) {
                fwd_soff = soff;
                break;
            }
        }
        if (fwd_soff == -1) return;
        //HBN_LOG("find target offset pair %d, %d", fwd_qoff, fwd_soff);
        int a_qb = pca->qoff, a_qe = fwd_qoff, a_sb = pca->soff, a_se = fwd_soff;
        double a_pi = 0.0;
        if (!s_realign_subseq(tbck_data, pca->qdir, fwd_query, rev_query, query_size,
            subject, pca->sid, pca->qend - pca->qoff, pca->send - pca->soff, &a_qb, &a_qe, &a_sb, &a_se, &a_pi)) {
            return;
        }
        if (a_pi < min_perc_identity) return;

        PoreCAlign new_pca = *pca;
        new_pca.qend = fwd_qoff;
        new_pca.enzyme_qend = fwd_qoff;
        new_pca.send = fwd_soff;
        new_pca.enzyme_send = fwd_soff;
        new_pca.pi = a_pi;
        if(set_pca_chain_offset(&new_pca, qvep_list->enzyme->enzyme_size)) {
            hbn_assert(new_pca.chain_qoff < new_pca.chain_qend);
            hbn_assert(new_pca.soff < new_pca.send);
            added_pca_list.push_back(new_pca);
            perfect_match_offset_point_list[fwd_qoff] = 1;
            //HBN_LOG("find extended pca");
            //dump_pca(fprintf, stderr, new_pca, -1);
        }      
    }
}

static void
s_extend_pefect_left_end(HbnTracebackData* tbck_data,
    const u8* fwd_query,
    const u8* rev_query,
    const u8* subject,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const double min_perc_identity,
    const PoreCAlign* org_pca_a,
    const int org_pca_c,
    PoreCAlign* pca,
    vector<PoreCAlign>& added_pca_list)
{
    u8* perfect_match_offset_point_list = kv_data(qvep_list->perfect_match_offset_point_list);
    u8* perfect_match_end_point_list = kv_data(qvep_list->perfect_match_end_point_list);
    const int* s_vdfi_a = reloci_list->reloci_array
                          +
                          reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_offset;
    const int s_vdfi_c = reloci_list->seq_reloci_info_array[pca->sid].enzyme_loci_cnt;
    const int query_size = qvep_list->query_size;
    const int enzyme_size = qvep_list->enzyme->enzyme_size;

    if (pca->qdir == REV) {
        int offset = query_size - pca->enzyme_qoff;
        if (offset >= enzyme_size) offset -= enzyme_size;
        hbn_assert(perfect_match_end_point_list[offset]);
        if (perfect_match_offset_point_list[offset]) return;
        int fwd_qoff = -1;
        for (int x = offset + 1; x < query_size - 1 && x < offset + 100; ++x) {
            if (perfect_match_offset_point_list[x]) {
                if (!perfect_match_end_point_list[x]) {
                    fwd_qoff = x;
                    break;
                }
            }
        }
        if (fwd_qoff == -1) return;
        int rev_qoff = query_size - fwd_qoff;
        if (rev_qoff >= enzyme_size) rev_qoff -= enzyme_size;
        int s_v_i = -1, s_v_c = 0;
        s_v_i = offset_to_enzyme_intv_idx(s_vdfi_a, s_vdfi_c, pca->soff, &s_v_c);
        int fwd_soff = -1;
        for (int x = s_v_i; x >= 0; --x) {
            int soff = s_vdfi_a[x];
            if (pca->soff <= soff) continue;
            if (pca->soff - soff > 100) break;
            if (s_is_ddfs_supported(rev_qoff, pca->qend, soff, pca->send)) {
                fwd_soff = soff;
                break;
            }
        }
        if (fwd_soff == -1) return;
        //HBN_LOG("find target offset pair %d, %d", rev_qoff, fwd_soff);
        int a_qb = rev_qoff, a_qe = pca->qend, a_sb = fwd_soff, a_se = pca->send;
        double a_pi = 0.0;
        if (!s_realign_subseq(tbck_data, pca->qdir, fwd_query, rev_query, query_size,
            subject, pca->sid, pca->qend - pca->qoff, pca->send - pca->soff, &a_qb, &a_qe, &a_sb, &a_se, &a_pi)) {
            return;
        }
        if (a_pi < min_perc_identity) return;

        PoreCAlign new_pca = *pca;
        new_pca.qoff = rev_qoff;
        new_pca.enzyme_qoff = rev_qoff;
        new_pca.soff = fwd_soff;
        new_pca.enzyme_soff = fwd_soff;
        new_pca.pi = a_pi;
        if(set_pca_chain_offset(&new_pca, qvep_list->enzyme->enzyme_size)) {
            hbn_assert(new_pca.chain_qoff < new_pca.chain_qend);
            hbn_assert(new_pca.soff < new_pca.send);
            added_pca_list.push_back(new_pca);
            perfect_match_offset_point_list[fwd_qoff] = 1;
            //HBN_LOG("find extended pca");
            //dump_pca(fprintf, stderr, new_pca, -1);
        }
    } else {
        int offset = pca->enzyme_qoff;
        hbn_assert(perfect_match_offset_point_list[offset]);
        if (perfect_match_end_point_list[offset]) return;
        int fwd_qoff = -1;
        for (int x = fwd_qoff - 1; x > 0 && x > fwd_qoff - 100; --x) {
            if (perfect_match_end_point_list[x]) {
                if (!perfect_match_offset_point_list[x]) {
                    fwd_qoff = x;
                    break;
                }
            }
        }
        if (fwd_qoff == -1) return;
        int s_v_i = -1, s_v_c = 0;
        s_v_i = offset_to_enzyme_intv_idx(s_vdfi_a, s_vdfi_c, pca->soff, &s_v_c);
        int fwd_soff = -1;
        for (int x = s_v_i; x >= 0; --x) {
            int soff = s_vdfi_a[x];
            if (pca->soff <= soff) continue;
            if (pca->soff - soff > 100) break;
            if (s_is_ddfs_supported(fwd_qoff, pca->qend, soff, pca->send)) {
                fwd_soff = soff;
                break;
            }
        }
        if (fwd_soff == -1) return;
        //HBN_LOG("find target offset pair %d, %d", fwd_qoff, fwd_soff);
        int a_qb = fwd_qoff, a_qe = pca->qend, a_sb = fwd_soff, a_se = pca->send;
        double a_pi = 0.0;
        if (!s_realign_subseq(tbck_data, pca->qdir, fwd_query, rev_query, query_size,
            subject, pca->sid, pca->qend - pca->qoff, pca->send - pca->soff, &a_qb, &a_qe, &a_sb, &a_se, &a_pi)) {
            return;
        }
        if (a_pi < min_perc_identity) return;   
        PoreCAlign new_pca = *pca;
        new_pca.qoff = fwd_qoff;
        new_pca.enzyme_qoff = fwd_qoff;
        new_pca.soff = fwd_soff;
        new_pca.enzyme_soff = fwd_soff;
        new_pca.pi = a_pi; 
        if(set_pca_chain_offset(&new_pca, qvep_list->enzyme->enzyme_size)) {
            hbn_assert(new_pca.chain_qoff < new_pca.chain_qend);
            hbn_assert(new_pca.soff < new_pca.send);
            added_pca_list.push_back(new_pca);
            perfect_match_offset_point_list[fwd_qoff] = 1;
            //HBN_LOG("find extended pca");
            //dump_pca(fprintf, stderr, new_pca, -1);
        }   
    }
}

static void
s_set_perfect_points(QueryVdfEndPointList* qvep_list, PoreCAlign* pca_a, int pca_c)
{
    const int max_offset = qvep_list->query_size + 1;
    kv_resize(u8, qvep_list->perfect_match_offset_point_list, max_offset);
    kv_resize(u8, qvep_list->perfect_match_end_point_list, max_offset);
    u8* perfect_match_offset_point_list = kv_data(qvep_list->perfect_match_offset_point_list);
    u8* perfect_match_end_point_list = kv_data(qvep_list->perfect_match_end_point_list);
    const int query_size = qvep_list->query_size;
    const int enzyme_size = qvep_list->enzyme->enzyme_size;

    fill(perfect_match_offset_point_list, perfect_match_offset_point_list + max_offset, 0);
    fill(perfect_match_end_point_list, perfect_match_end_point_list + max_offset, 0);
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        if (pca->enzyme_qoff > 0 && pca->enzyme_qoff == pca->qoff) {
            if (pca->qdir == FWD) {
                perfect_match_offset_point_list[pca->enzyme_qoff] = 1;
            } else {
                int offset = query_size - pca->enzyme_qoff;
                if (offset >= enzyme_size) offset -= enzyme_size;
                perfect_match_end_point_list[offset] = 1;
            }
        }
        if (pca->enzyme_qend != query_size && pca->enzyme_qend == pca->qend) {
            if (pca->qdir == FWD) {
                perfect_match_end_point_list[pca->enzyme_qend] = 1;
            } else {
                int offset = query_size - pca->enzyme_qend;
                if (offset >= enzyme_size) offset -= enzyme_size;
                perfect_match_offset_point_list[offset] = 1;
            }
        }
    }
}

void
fix_imperfect_ends(HbnTracebackData* tbck_data,
    const u8* fwd_query,
    const u8* rev_query,
    const u8* subject,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    PoreCAlign* pca_a,
    int pca_c,
    const double min_perc_identity,
    std::vector<PoreCAlign>& new_pca_list)
{
    const int verbose = 0;
    s_find_broken_perfect_points(qvep_list, pca_a, pca_c);

    int r = 0;
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        r = (pca->enzyme_qoff == 0 && pca->enzyme_qend == pca->qend)
            ||
            (pca->enzyme_qend == pca->qsize && pca->qoff == pca->enzyme_qoff);
        if (r) {
            if (verbose) {
                HBN_LOG("head or tail pca");
                dump_pca(fprintf, stderr, *pca, i);
            }
            continue;
        }
        
        r = pca->is_perfect;
        if (r) {
            if (verbose) {
                HBN_LOG("perfect pca");
                dump_pca(fprintf, stderr, *pca, i);
            }
            continue;
        }

        r = pca->lp;
        if (!r) {
            if (verbose) {
                HBN_LOG("left imperfect pca");
                dump_pca(fprintf, stderr, *pca, i);
            }
            s_fix_left_imperfect_pca(tbck_data, fwd_query, rev_query, subject,
                reloci_list, qvep_list, min_perc_identity, pca_a, pca_c, pca, new_pca_list);
            continue;            
        }

        r = pca->rp;
        if (!r) {
            if (verbose) {
                HBN_LOG("right imperfect pca");
                dump_pca(fprintf, stderr, *pca, i);
            }
            s_fix_right_imperfect_pca(tbck_data, fwd_query, rev_query, subject,
                reloci_list, qvep_list, min_perc_identity, pca_a, pca_c, pca, new_pca_list);
            continue;            
        }
    }

    new_pca_list.assign(pca_a, pca_a + pca_c);
    pca_a = new_pca_list.data();
    pca_c = new_pca_list.size();
    //remove_duplicate_pca(pca_a, &pca_c);
    int n = 0;
    for (int i = 0; i < pca_c; ++i) if (pca_a[i].qoff < pca_a[i].qend) pca_a[n++] = pca_a[i];
    pca_c = n;
    new_pca_list.resize(pca_c);

    vector<PoreCAlign> added_pca_list;
    pca_a = new_pca_list.data();
    pca_c = new_pca_list.size();
    s_set_perfect_points(qvep_list, pca_a, pca_c);
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        if (pca->enzyme_qend != pca->qsize && pca->qend == pca->enzyme_qend) {
            //HBN_LOG("extend perfect right end for");
            //dump_pca(fprintf, stderr, *pca, i);
            s_extend_pefect_right_end(tbck_data,
                fwd_query,
                rev_query,
                subject,
                reloci_list,
                qvep_list,
                min_perc_identity,
                pca_a,
                pca_c,
                pca,
                added_pca_list);
        }
        if (pca->enzyme_qoff > 0 && pca->enzyme_qoff == pca->qoff) {
            //HBN_LOG("extend pefect left end for");
            //dump_pca(fprintf, stderr, *pca, i);
            s_extend_pefect_left_end(tbck_data,
                fwd_query,
                rev_query,
                subject,
                reloci_list,
                qvep_list,
                min_perc_identity,
                pca_a,
                pca_c,
                pca,
                added_pca_list);            
        }
    }
    new_pca_list.insert(new_pca_list.end(), added_pca_list.begin(), added_pca_list.end());
}
