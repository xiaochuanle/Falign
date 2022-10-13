#include "extend_one_chain.hpp"

#include "gap_filling.hpp"

#include <algorithm>
#include <cmath>

using namespace std;

static const int kMatLen = 8;
static const int kExtendBlock = 50;
static const int kMaxExtendBlock = 90;
static const int kMaxW = 50;
static const double kEpsilon = 0.2;
static const int kEnzymeNbhd = 50;

BOOL 
refine_pi_for_one_pca(HbnTracebackData* tbck_data,
    const u8* fwd_read,
    const u8* rev_read,
    const int read_size,
    const u8* subject,
    PoreCAlign* pca)
{
    const u8* read = (pca->qdir == FWD) ? fwd_read : rev_read;
    const u8* q = read + pca->qoff;
    int ql = pca->qend - pca->qoff;
    const u8* s = subject + pca->soff;
    int sl = pca->send - pca->soff;
    int dist = max(ql, sl) * 0.5;
    int x = edlib_nw_dist(tbck_data->edlib, q, ql, s, sl, dist);
    if (x == -1) return FALSE;

    double pi = 100.0 * (1.0 - 2.0 * x / (ql + sl));
    pca->pi = pi;
    return TRUE;
}

double
compute_pi_for_one_align(const char* qas1,
    const char* sas1,
    const int as_size1,
    const char* qas2,
    const char* sas2,
    const int as_size2)
{
    int mat = 0, len = 0;
    for (int i = 0; i < as_size1; ++i, ++len) if (qas1[i] == sas1[i]) ++mat;
    for (int i = 0; i < as_size2; ++i, ++len) if (qas2[i] == sas2[i]) ++mat;
    if (!len) return 0.0;
    return 100.0 * mat / len;
}

///// extension
static int
overhang_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    int* _qend,
    int* _send,
    kstring_t* qaln,
    kstring_t* saln)
{
    *_qend = 0;
    *_send = 0;
    ks_clear(*qaln);
    ks_clear(*saln);
    int qoff, qend, soff, send;
    double ident_perc;
    int r = dalign_align(dalign,
                query,
                0,
                query_length,
                subject,
                0,
                subject_length,
                1,
                10.0,
                &qoff,
                &qend,
                &soff,
                &send,
                &ident_perc,
                NULL,
                NULL);
    if (!r) return r;
    hbn_assert(qoff == 0 && soff == 0);

    r = edlib_nw(edlib,
            query + qoff,
            qend - qoff,
            subject + soff,
            send - soff,
            -1,
            qaln,
            saln);
    *_qend = qend;
    *_send = send;
    return r;
}

static void
right_extend_align(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    int* qend,
    const int query_length,
    const u8* subject,
    int* send,
    const int subject_length,
    kstring_t* qabuf,
    kstring_t* sabuf,
    char** qae,
    char** sae)
{
    *qend -= kMatLen; *send -= kMatLen; *qae -= kMatLen; *sae -= kMatLen;
    int qt = *qend;
    int st = *send;
    bool last_block = false;
    bool done = false;
    while (1) {
        int qblk, sblk;
        if (query_length - qt <= kMaxExtendBlock || subject_length - st <= kMaxExtendBlock) {
            qblk = (subject_length - st) * 1.2;
            qblk = min<int>(qblk, query_length - qt);
            sblk = (query_length - qt) * 1.2;
            sblk = min<int>(sblk, subject_length - st);
            last_block = true;
        } else {
            qblk = kExtendBlock;
            sblk = kExtendBlock;
            last_block = false;
        }
        if (qblk == 0 || sblk == 0) break;
    //HBN_LOG("qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, query_length, subject_length);
        int qfae = 0, sfae = 0;
        //edlib_shw(edlib, query + qt, qblk, subject + st, sblk, &qfae, &sfae, qabuf, sabuf);
        overhang_extend(dalign, edlib, query + qt, qblk, subject + st, sblk, &qfae, &sfae, qabuf, sabuf);
        if (qblk - qfae > 30 || sblk - sfae > 30) done = true;
        int acnt = 0, qcnt = 0, scnt = 0;
        int as_size = ks_size(*qabuf), k = as_size - 1, m = 0;
        while (k >= 0) {
            char qc = ks_A(*qabuf, k);
            char sc = ks_A(*sabuf, k);
            if (qc != GAP_CHAR) ++qcnt;
            if (sc != GAP_CHAR) ++scnt;
            m = (qc == sc) ? (m+1) : 0;
            ++acnt;
            if (m == kMatLen) break;
            --k;
        }
        //HBN_LOG("qblk = %d, qfae = %d, sblk = %d, sfae = %d, m = %d, k = %d, qcnt = %d, scnt = %d, acnt = %d", qblk, qfae, sblk, sfae, m, k, qcnt, scnt, acnt);
        if (m != kMatLen || k < 1) {
            as_size = 0;
            for (int i = 0; i < qblk && i < sblk; ++i, ++qt, ++st) {
                int qc = query[qt];
                int sc = subject[st];
                qc = DECODE_RESIDUE(qc);
                sc = DECODE_RESIDUE(sc);
                if (qc != sc) break;
                **qae = qc; ++(*qae);
                **sae = sc; ++(*sae);
                ++as_size;
            }
            done = true;
        } else {
            as_size -= acnt;
            qt += (qfae - qcnt);
            st += (sfae - scnt);
            if (done) {
                as_size += kMatLen;
                qt += kMatLen;
                st += kMatLen;
            }
            memcpy(*qae, ks_s(*qabuf), as_size);
            memcpy(*sae, ks_s(*sabuf), as_size);
            *(qae) += as_size;
            *(sae) += as_size;
        }
        //HBN_LOG("qt = %d, st = %d", qt, st);
        if (done) break;
    }
    *qend = qt;
    *send = st;
}

static void
left_extend_align(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    int* qoff,
    const u8* subject,
    int* soff,
    kstring_t* qabuf,
    kstring_t* sabuf,
    char** qas,
    char** sas)
{
    *qoff += kMatLen; *soff += kMatLen; *qas += kMatLen; *sas += kMatLen;

    int qf = *qoff;
    int sf = *soff;
    bool last_block = false;
    bool done = false;
    vector<u8> qfrag, sfrag;
    while (1) {
        int qblk, sblk;
        if (qf <= kMaxExtendBlock || sf <= kMaxExtendBlock) {
            qblk = sf * 1.2;
            qblk = min<int>(qblk, qf);
            sblk = qf * 1.2;
            sblk = min<int>(sblk, sf);
            last_block = true;
        } else {
            qblk = kExtendBlock;
            sblk = kExtendBlock;
            last_block = false;
        }
        if (qblk == 0 || sblk == 0) break;
        qfrag.clear(); for (int i = 1; i <= qblk; ++i) qfrag.push_back(query[qf-i]);
        sfrag.clear(); for (int i = 1; i <= sblk; ++i) sfrag.push_back(subject[sf-i]);

        int qfae = 0, sfae = 0;
        edlib_shw(edlib, qfrag.data(), qblk, sfrag.data(), sblk, &qfae, &sfae, qabuf, sabuf);
        if (qblk - qfae > 30 || sblk - sfae > 30) done = true;
        int acnt = 0, qcnt = 0, scnt = 0;
        int as_size = ks_size(*qabuf), k = as_size - 1, m = 0;
        while (k >= 0) {
            char qc = ks_A(*qabuf, k);
            char sc = ks_A(*sabuf, k);
            if (qc != GAP_CHAR) ++qcnt;
            if (sc != GAP_CHAR) ++scnt;
            m = (qc == sc) ? (m+1) : 0;
            ++acnt;
            if (m == kMatLen) break;
            --k;
        }
        //HBN_LOG("m = %d, k = %d, qcnt = %d, scnt = %d, acnt = %d", m, k, qcnt, scnt, acnt);
        if (m != kMatLen || k < 1) {
            for (int i = 0; i < qblk && i < sblk; ++i) {
                int qc = query[qf-1];
                int sc = subject[sf-1];
                qc = DECODE_RESIDUE(qc);
                sc = DECODE_RESIDUE(sc);
                if (qc != sc) break;
                --(*qas); **qas = qc;
                --(*sas); **sas = sc;
                --qf; --sf;
            }
            done = true;
        } else {
            as_size -= acnt;
            qf -= (qfae - qcnt);
            sf -= (sfae - scnt);
            if (done) {
                as_size += kMatLen;
                qf -= kMatLen;
                sf -= kMatLen;
            }
            for (int i = 0; i < as_size; ++i) {
                int qc = ks_A(*qabuf, i);
                --(*qas); **qas = qc;
                int sc = ks_A(*sabuf, i);
                --(*sas); **sas = sc;
            }
        }
        //HBN_LOG("sf = %d, qf = %d", qf, sf);
        if (done) break;
    }
    *qoff = qf;
    *soff = sf;
}

///// fix enzyme ends

struct enzyme_offset {
    int qoff, soff;
    int enzyme_qoff, enzyme_soff;

    enzyme_offset() {
        qoff = -1;
        soff = -1;
        enzyme_qoff = -1;
        enzyme_soff = -1;
    }

    void dump(FILE* out, const char* sep = NULL) {
        fprintf(out, "[%d (%d), %d (%d)]", qoff, enzyme_qoff, soff, enzyme_soff);
        if (sep) fprintf(out, "%s", sep);
    }
};

static void
s_add_nearby_enzyme_pos(const int* reloci_array,
    const int intv_cnt,
    const int intv_i,
    const int offset,
    const int min_enzyme_offset,
    const int max_enzyme_offset,
    const int max_dist,
    vector<int>& offset_list)
{
    offset_list.clear();
    hbn_assert(offset >= reloci_array[intv_i] && offset <= reloci_array[intv_i + 1]);
    hbn_assert(offset >= min_enzyme_offset && offset <= max_enzyme_offset);

    int i = intv_i;
    while (i >= 0) {
        if (offset - reloci_array[i] > max_dist) break;
        if (reloci_array[i] < min_enzyme_offset) break;
        offset_list.push_back(reloci_array[i]);
        --i;
    }

    i = intv_i;
    while (i < intv_cnt) {
        if (reloci_array[i+1] - offset > max_dist) break;
        if (reloci_array[i+1] > max_enzyme_offset) break;
        offset_list.push_back(reloci_array[i+1]);
        ++i;
    }
    sort(offset_list.begin(), offset_list.end());
}

static bool 
s_resolve_unfixed_start_offset(HbnTracebackData* tbck_data,
    const u8* query,
    const u8* subject,
    enzyme_offset* eop)
{
    hbn_assert(eop->qoff == -1 || eop->soff == -1);
    if (eop->qoff == -1) {
        hbn_assert(eop->soff < tbck_data->send);
        if (eop->soff >= tbck_data->soff) {
            int qi = tbck_data->qoff;
            int si = tbck_data->soff;
            int ai = 0;
            while (si < eop->soff) {
                if (tbck_data->qas[ai] != GAP_CHAR) ++qi;
                if (tbck_data->sas[ai] != GAP_CHAR) ++si;
                ++ai;
            }
            hbn_assert(si == eop->soff);
            hbn_assert(qi < tbck_data->qend);
            eop->qoff = qi;
            return true;
        }
        hbn_assert(eop->soff < tbck_data->soff);
        int sf = eop->soff;
        int st = tbck_data->soff + kMatLen;
        hbn_assert(sf < st);
        int qt = tbck_data->qoff + kMatLen;
        int qf = qt - (st - sf) - 50; qf = max(0, qf);
        hbn_assert(sf >= 0 && sf < st && st <= tbck_data->send);
        hbn_assert(qf >= 0 && qf < qt && qt <= tbck_data->qend);
        vector<u8> qsubseq, ssubseq;
        for (int i = qt - 1; i >= qf; --i) qsubseq.push_back(query[i]);
        for (int i = st - 1; i >= sf; --i) ssubseq.push_back(subject[i]);
        int qe = 0, se = 0;
        //HBN_LOG("1 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        overhang_extend(tbck_data->dalign, tbck_data->edlib, qsubseq.data(), qsubseq.size(),
            ssubseq.data(), ssubseq.size(), &qe, &se, &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
        if (ssubseq.size() != se) return false;
        eop->qoff = tbck_data->qoff + kMatLen - qe;
        hbn_assert(eop->qoff >= 0);
        return true;
    }
    if (eop->soff == -1) {
        hbn_assert(eop->qoff < tbck_data->qend);
        if (eop->qoff >= tbck_data->qoff) {
            int qi = tbck_data->qoff;
            int si = tbck_data->soff;
            int ai = 0;
            while (qi < eop->qoff) {
                if (tbck_data->qas[ai] != GAP_CHAR) ++qi;
                if (tbck_data->sas[ai] != GAP_CHAR) ++si;
                ++ai;
            }
            hbn_assert(qi == eop->qoff);
            hbn_assert(si < tbck_data->send);
            eop->soff = si;
            return true;
        }

        int qf = eop->qoff;
        int qt = tbck_data->qoff + kMatLen;
        hbn_assert(qf < qt);
        int st = tbck_data->soff + kMatLen;
        int sf = st - (qt - qf) - 50; sf = max(0, sf);
        hbn_assert(sf >= 0 && sf < st && st <= tbck_data->send);
        hbn_assert(qf >= 0 && qf < qt && qt <= tbck_data->qend);
        vector<u8> qsubseq, ssubseq;
        int qe = 0, se = 0;
        //HBN_LOG("2 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        overhang_extend(tbck_data->dalign, tbck_data->edlib, qsubseq.data(), qsubseq.size(),
            ssubseq.data(), ssubseq.size(), &qe, &se, &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
        if (qsubseq.size() != qe) return false;
        eop->soff = tbck_data->soff + kMatLen - se;
        return true;      
    }
    return false;
}

static void
s_add_start_offset_pairs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const u8* query,
    const u8* subject,
    int qdir,
    int sid,
    const int qoff,
    const int soff,
    const int b_qend,
    const int b_send,
    vector<enzyme_offset>& leop_list)
{
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
    int l_q_intv_i = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, qoff, &l_q_intv_c);  
    int l_s_intv_c = 0;
    int l_s_intv_i = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, soff, &l_s_intv_c);

    vector<int> ql_list;
    vector<int> sl_list;
    s_add_nearby_enzyme_pos(q_reloci_array, l_q_intv_c, l_q_intv_i, qoff, 0, b_qend - 1, kEnzymeNbhd, ql_list);
    s_add_nearby_enzyme_pos(s_reloci_array, l_s_intv_c, l_s_intv_i, soff, 0, b_send - 1, kEnzymeNbhd, sl_list);
    hbn_assert(b_qend <= tbck_data->qend);
    hbn_assert(b_send <= tbck_data->send);
    for (auto qi : ql_list) hbn_assert(qi < b_qend);
    for (auto si : sl_list) hbn_assert(si < b_send);

    for (auto qi : ql_list) {
        int q_d = qi - qoff;
        for (auto si : sl_list) {
            int s_d = si - soff;
            if (q_d != s_d) continue;
            enzyme_offset eop;
            eop.qoff = qi;
            eop.soff = si;
            eop.enzyme_qoff = qi;
            eop.enzyme_soff = si;
            leop_list.push_back(eop);
        }
    } 
    if (!leop_list.empty()) {
        int min_d = numeric_limits<int>::max();
        enzyme_offset eop;
        for (auto& x : leop_list) {
            int d = abs(x.qoff - qoff);
            if (d < min_d) {
                min_d = d;
                eop = x;
            }
        }
        leop_list.clear();
        leop_list.push_back(eop);
    }
    if (!leop_list.empty()) return;

    for (auto qi : ql_list) {
        int q_d = b_qend - qi;
        hbn_assert(q_d > 0);
        for (auto si : sl_list) {
            int s_d = b_send - si;
            hbn_assert(s_d > 0);
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            enzyme_offset eop;
            eop.qoff = qi;
            eop.soff = si;
            eop.enzyme_qoff = qi;
            eop.enzyme_soff = si;
            leop_list.push_back(eop);
        }
    } 
    if (!leop_list.empty()) return;

#if 0
    int min_q_d = numeric_limits<int>::max();
    int best_qoff = -1;
    for (auto qi : ql_list) {
        int q_d = abs(qi - qoff);
        if (q_d < min_q_d) {
            min_q_d = q_d;
            best_qoff = qi; 
        }
    }
    if (best_qoff >= 0) {
        enzyme_offset eop;
        eop.qoff = best_qoff;
        eop.enzyme_qoff = best_qoff;
        eop.soff = -1;
        eop.enzyme_soff = ((soff - s_reloci_array[l_s_intv_i]) <= (s_reloci_array[l_s_intv_i+1] - soff)) ? s_reloci_array[l_s_intv_i] : s_reloci_array[l_s_intv_i+1];
        //fprintf(stderr, "fix soff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_start_offset(tbck_data, query, subject, &eop)) {
            //fprintf(stderr, "fixed soff for\t"); eop.dump(stderr, "\n");
            leop_list.push_back(eop);
        }
    }

    int min_s_d = numeric_limits<int>::max();
    int best_soff = -1;
    for (auto si : sl_list) {
        int s_d = abs(si - soff);
        if (s_d < min_s_d) {
            min_s_d = s_d;
            best_soff = si;
        }
    }
    if (best_soff >= 0) {
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qoff - q_reloci_array[l_q_intv_i]) <= (q_reloci_array[l_q_intv_i+1] - qoff)) ? q_reloci_array[l_q_intv_i] : q_reloci_array[l_q_intv_i+1];
        eop.soff = best_soff;
        eop.enzyme_soff = best_soff;
        //fprintf(stderr, "fix qoff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_start_offset(tbck_data, query, subject, &eop)) {
            //fprintf(stderr, "fixed qoff for\t"); eop.dump(stderr, "\n");
            leop_list.push_back(eop);
        }
    }
#else
    for (auto qi : ql_list) {
        enzyme_offset eop;
        eop.qoff = qi;
        eop.enzyme_qoff = qi;
        eop.soff = -1;
        eop.enzyme_soff = ((soff - s_reloci_array[l_s_intv_i]) <= (s_reloci_array[l_s_intv_i+1] - soff)) ? s_reloci_array[l_s_intv_i] : s_reloci_array[l_s_intv_i+1];
        //fprintf(stderr, "fix soff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_start_offset(tbck_data, query, subject, &eop)) {
            //fprintf(stderr, "fixed soff for\t"); eop.dump(stderr, "\n");
            leop_list.push_back(eop);
        }
    }

    for (auto si :sl_list) {
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qoff - q_reloci_array[l_q_intv_i]) <= (q_reloci_array[l_q_intv_i+1] - qoff)) ? q_reloci_array[l_q_intv_i] : q_reloci_array[l_q_intv_i+1];
        eop.soff = si;
        eop.enzyme_soff = si;
        //fprintf(stderr, "fix qoff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_start_offset(tbck_data, query, subject, &eop)) {
            //fprintf(stderr, "fixed qoff for\t"); eop.dump(stderr, "\n");
            leop_list.push_back(eop);
        }     
    }
#endif 

    if (leop_list.empty()) {
        enzyme_offset eop;
        eop.qoff = qoff;
        eop.enzyme_qoff = ((qoff - q_reloci_array[l_q_intv_i]) <= (q_reloci_array[l_q_intv_i+1] - qoff)) ? q_reloci_array[l_q_intv_i] : q_reloci_array[l_q_intv_i+1];
        eop.soff = soff;
        eop.enzyme_soff = ((soff - s_reloci_array[l_s_intv_i]) <= (s_reloci_array[l_s_intv_i+1] - soff)) ? s_reloci_array[l_s_intv_i] : s_reloci_array[l_s_intv_i+1];
        leop_list.push_back(eop);
    }
}

static bool 
s_resolve_unfixed_end_offset(HbnTracebackData* tbck_data,
    const u8* query,
    const int query_size,
    const u8* subject,
    const int subject_size,
    enzyme_offset* eop)
{   
    hbn_assert(eop->qoff == -1 || eop->soff == -1);
    if (eop->qoff == -1) {
        hbn_assert(eop->soff > tbck_data->soff);
        if (eop->soff <= tbck_data->send) {
            int qi = tbck_data->qoff;
            int si = tbck_data->soff;
            int ai = 0;
            while (si < eop->soff) {
                if (tbck_data->qas[ai] != GAP_CHAR) ++qi;
                if (tbck_data->sas[ai] != GAP_CHAR) ++si;
                ++ai;
            }
            hbn_assert(si == eop->soff);
            hbn_assert(qi <= tbck_data->qend);
            eop->qoff = qi;
            return true;
        }
        int sf = tbck_data->send - kMatLen;
        int st = eop->soff;
        hbn_assert(sf < st);
        int qf = tbck_data->qend - kMatLen;
        int qt = qf + (st - sf) + 50; qt = min(qt, query_size);
        hbn_assert(sf >= tbck_data->soff && sf < st && st <= subject_size);
        hbn_assert(qf >= tbck_data->qoff && qf < qt && qt <= query_size);
        int qe = 0, se = 0;
        //HBN_LOG("3 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        overhang_extend(tbck_data->dalign, tbck_data->edlib, query + qf, qt - qf, subject + sf, st - sf, &qe, &se, &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
        if (se != (st - sf)) return 0;
        eop->qoff = tbck_data->qend - kMatLen + qe;
        return true;
    }
    if (eop->soff == -1) {
        hbn_assert(eop->qoff > tbck_data->qoff);
        if (eop->qoff <= tbck_data->qend) {
            int qi = tbck_data->qoff;
            int si = tbck_data->soff;
            int ai = 0;
            while (qi < eop->qoff) {
                if (tbck_data->qas[ai] != GAP_CHAR) ++qi;
                if (tbck_data->sas[ai] != GAP_CHAR) ++si;
                ++ai;
            }
            hbn_assert(qi == eop->qoff);
            hbn_assert(si <= tbck_data->send);
            eop->soff = si;
            return true;            
        }
        int qf = tbck_data->qend - kMatLen;
        int qt = eop->qoff;
        hbn_assert(qf < qt);
        int sf = tbck_data->send - kMatLen;
        int st = sf + (qt - qf) + 50; st = min(st, subject_size);
        hbn_assert(sf >= tbck_data->soff && sf < st && st <= subject_size);
        hbn_assert(qf >= tbck_data->qoff && qf < qt && qt <= query_size);
        int qe = 0, se = 0;
        //HBN_LOG("4 qf = %d, qt = %d, sf = %d, st = %d", qf, qt, sf, st);
        overhang_extend(tbck_data->dalign, tbck_data->edlib, query + qf, qt - qf, subject + sf, st - sf, &qe, &se, &tbck_data->ext_qabuf, &tbck_data->ext_sabuf);
        if (qe != (qt - qf)) return 0;
        eop->soff = tbck_data->send - kMatLen + se;
        return true;
    }
    return false;
}

static void
s_add_end_offset_pairs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const u8* query,
    const int query_size,
    const u8* subject,
    const int subject_size,
    int qdir,
    int sid,
    const int qend,
    const int send,
    const int b_qoff,
    const int b_soff,
    vector<enzyme_offset>& reop_list)
{
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
    int r_q_intv_i = offset_to_enzyme_intv_idx(q_reloci_array, q_reloci_cnt, qend, &r_q_intv_c);  
    int r_s_intv_c = 0;
    int r_s_intv_i = offset_to_enzyme_intv_idx(s_reloci_array, s_reloci_cnt, send, &r_s_intv_c);
    vector<int> qr_list;
    vector<int> sr_list;
    s_add_nearby_enzyme_pos(q_reloci_array, r_q_intv_c, r_q_intv_i, qend, b_qoff + 1, query_size, kEnzymeNbhd, qr_list);
    s_add_nearby_enzyme_pos(s_reloci_array, r_s_intv_c, r_s_intv_i, send, b_soff + 1, subject_size, kEnzymeNbhd, sr_list);

    for (auto qi : qr_list) hbn_assert(qi > b_qoff);
    for (auto si : sr_list) hbn_assert(si > b_soff);

    //fprintf(stderr, "qoff: "); for (auto qi : qr_list) fprintf(stderr, "%d\t", qi); fprintf(stderr, "\n");
    //fprintf(stderr, "soff: "); for (auto si : sr_list) fprintf(stderr, "%d\t", si); fprintf(stderr, "\n");

    for (auto qi : qr_list) {
        int q_d = qi - qend;
        for (auto si : sr_list) {
            int s_d = si - send;
            if (q_d != s_d) continue;
            enzyme_offset eop;
            eop.qoff = qi;
            eop.enzyme_qoff = qi;
            eop.soff = si;
            eop.enzyme_soff = si;
            reop_list.push_back(eop);
        }
    } 
    if (!reop_list.empty()) {
        int min_d = numeric_limits<int>::max();
        enzyme_offset eop;
        for (auto& x : reop_list) {
            int d = abs(x.qoff - qend);
            if (d < min_d) {
                min_d = d;
                eop = x;
            }
        }
        reop_list.clear();
        reop_list.push_back(eop);
    }
    if (!reop_list.empty()) return;

    for (auto qi : qr_list) {
        int q_d = qi - b_qoff;
        //HBN_LOG("qi = %d, q_d = %d, b_qoff = %d", qi, q_d, b_qoff);
        hbn_assert(q_d > 0);
        for (auto si : sr_list) {
            int s_d = si - b_soff;
            //HBN_LOG("\tsi = %d, s_d = %d, b_soff = %d", si, s_d, b_soff);
            hbn_assert(s_d > 0);
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            enzyme_offset eop;
            eop.qoff = qi;
            eop.enzyme_qoff = qi;
            eop.soff = si;
            eop.enzyme_soff = si;
            reop_list.push_back(eop);
        }
    } 
    if (!reop_list.empty()) return;

#if 0
    int min_q_d = numeric_limits<int>::max();
    int best_qend = -1;
    for (auto qi : qr_list) {
        int q_d = abs(qi - qend);
        if (q_d < min_q_d) {
            min_q_d = q_d;
            best_qend = qi; 
        }
    }
    if (best_qend >= 0) {
        enzyme_offset eop;
        eop.qoff = best_qend;
        eop.enzyme_qoff = best_qend;
        eop.soff = -1;
        eop.enzyme_soff = ((send - s_reloci_array[r_s_intv_i]) < (s_reloci_array[r_s_intv_i+1] - send)) ? s_reloci_array[r_s_intv_i] : s_reloci_array[r_s_intv_i+1];
        //fprintf(stderr, "resolve soff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, &eop)) {
            //fprintf(stderr, "resolved soff for\t"); eop.dump(stderr, "\n");
            reop_list.push_back(eop);
        }
    }

    int min_s_d = numeric_limits<int>::max();
    int best_send = -1;
    for (auto si : sr_list) {
        int s_d = abs(si - send);
        if (s_d < min_s_d) {
            min_s_d = s_d;
            best_send = si;
        }
    }
    if (best_send >= 0) {
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qend - q_reloci_array[r_q_intv_i]) < (q_reloci_array[r_q_intv_i+1] - qend)) ? q_reloci_array[r_q_intv_i] : q_reloci_array[r_q_intv_i+1];
        eop.soff = best_send;
        eop.enzyme_soff = best_send;
        //fprintf(stderr, "resolve qoff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, &eop)) {
            //fprintf(stderr, "resolved qoff for\t"); eop.dump(stderr, "\n");
            reop_list.push_back(eop);
        }
    }
#else
    for (auto qi :qr_list) {
        enzyme_offset eop;
        eop.qoff = qi;
        eop.enzyme_qoff = qi;
        eop.soff = -1;
        eop.enzyme_soff = ((send - s_reloci_array[r_s_intv_i]) < (s_reloci_array[r_s_intv_i+1] - send)) ? s_reloci_array[r_s_intv_i] : s_reloci_array[r_s_intv_i+1];
        //fprintf(stderr, "resolve soff for\t"); eop.dump(stderr, "\n");
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, &eop)) {
            //fprintf(stderr, "resolved soff for\t"); eop.dump(stderr, "\n");
            reop_list.push_back(eop);
        }
    }

    for (auto si : sr_list) {
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qend - q_reloci_array[r_q_intv_i]) < (q_reloci_array[r_q_intv_i+1] - qend)) ? q_reloci_array[r_q_intv_i] : q_reloci_array[r_q_intv_i+1];
        eop.soff = si;
        eop.enzyme_soff = si;
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, &eop)) {
            reop_list.push_back(eop);
        }
    }
#endif 

    if (reop_list.empty()) {
        enzyme_offset eop;
        eop.qoff = qend;
        eop.enzyme_qoff = ((qend - q_reloci_array[r_q_intv_i]) < (q_reloci_array[r_q_intv_i+1] - qend)) ? q_reloci_array[r_q_intv_i] : q_reloci_array[r_q_intv_i+1];
        eop.soff = send;
        eop.enzyme_soff = ((send - s_reloci_array[r_s_intv_i]) < (s_reloci_array[r_s_intv_i+1] - send)) ? s_reloci_array[r_s_intv_i] : s_reloci_array[r_s_intv_i+1];
        reop_list.push_back(eop);
    }
}

static void
s_fix_enzyme_align_ends(const int query_id,
    const int query_dir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    const double perc_identity,
    const int b_qoff,
    const int b_qend,
    const int b_soff,
    const int b_send,
    HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    int* cov_stats,
    std::vector<PoreCAlign>& align_list)
{
    const int* q_reloci_array = (query_dir == FWD)
                                ?
                                kv_data(qvep_list->fwd_vdf_endpoint_list)
                                :
                                kv_data(qvep_list->rev_vdf_endpoint_list);
    const int q_reloci_cnt = kv_size(qvep_list->rev_vdf_endpoint_list);
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[subject_id].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[subject_id].enzyme_loci_cnt;
    const u8* query = (query_dir == FWD) ? fwd_query : rev_query;

    vector<enzyme_offset> leop_list;
    s_add_start_offset_pairs(tbck_data, reloci_list, qvep_list, query, subject, query_dir, subject_id, tbck_data->qoff, tbck_data->soff, b_qend, b_send, leop_list);
    //HBN_LOG("left enzyme offset pairs:");
    //for (auto& eop : leop_list) eop.dump(stderr, "\t"); fprintf(stderr, "\n");

    vector<enzyme_offset> reop_list;
    s_add_end_offset_pairs(tbck_data, reloci_list, qvep_list, query, query_size, subject, subject_size, query_dir, subject_id, tbck_data->qend, tbck_data->send, b_qoff, b_soff, reop_list);
    //HBN_LOG("right enzyme offset pairs:");
    //for (auto& eop : reop_list) eop.dump(stderr, "\t"); fprintf(stderr, "\n");

    PoreCAlign pca;
    pca.qid = query_id;
    pca.qdir = query_dir,
    pca.qsize = query_size;
    pca.sid = subject_id;
    pca.ssize = subject_size;
    for (auto& sx : leop_list) {
        pca.qoff = sx.qoff;
        pca.enzyme_qoff = sx.enzyme_qoff;
        pca.soff = sx.soff;
        pca.enzyme_soff = sx.enzyme_soff;
        for (auto& ex : reop_list) {
            pca.qend = ex.qoff;
            pca.enzyme_qend = ex.enzyme_qoff;
            pca.send = ex.soff;
            pca.enzyme_send = ex.enzyme_soff;
            if (pca.qoff >= pca.qend || pca.soff >= pca.send) continue;
            //dump_pca(fprintf, stderr, pca, -1);
            if (!refine_pi_for_one_pca(tbck_data, fwd_query, rev_query, query_size, subject, &pca)) continue;
            if (!set_pca_chain_offset(&pca, qvep_list->enzyme->enzyme_size)) continue;
            hbn_assert(pca.chain_qoff < pca.chain_qend);
            hbn_assert(pca.soff < pca.send);
            //dump_pca(fprintf, stderr, pca, -1);
            align_list.push_back(pca);
            int x = pca.qoff;
            int y = pca.qend;
            if (query_dir == REV) {
                x = query_size - pca.qend;
                y = query_size - pca.qoff;
                for (int i = x; i < y; ++i) ++cov_stats[i];
            }
        }
    }
}

//// main function

void
extend_one_chain(const int query_id,
    const int query_dir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    DDFS_Seed* sa,
    int sc,
    const double perc_identity,
    HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    int* cov_stats,
    std::vector<PoreCAlign>& align_list)
{
    const u8* query = (query_dir == FWD) ? fwd_query : rev_query;
    int r = porec_compute_traceback(tbck_data, query, query_size, subject, subject_size, sa, sc, 20, perc_identity);
    if (!r) return;
    HbnTracebackDataDump(fprintf, stderr, tbck_data);
    int b_qoff = tbck_data->qoff;
    int b_qend = tbck_data->qend;
    int b_soff = tbck_data->soff;
    int b_send = tbck_data->send;
    hbn_assert(b_qend <= query_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    hbn_assert(b_send <= subject_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    HbnTracebackData* data = tbck_data;
    int as_size = 0;

    right_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qend, query_size,
        subject, &tbck_data->send, subject_size, &tbck_data->ext_qabuf, &tbck_data->ext_sabuf, &tbck_data->qae, &tbck_data->sae);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qend >= b_qend, "qid = %d, qend = %d, b_qend = %d, send = %d, b_send = %d", query_id, 
        tbck_data->qend, b_qend, tbck_data->send, b_send);
    hbn_assert(tbck_data->send >= b_send);

    left_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qoff, subject, &tbck_data->soff,
        &tbck_data->ext_qabuf, &tbck_data->ext_sabuf, &tbck_data->qas, &tbck_data->sas);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qoff <= b_qoff);
    hbn_assert(tbck_data->soff <= b_soff);
        
    as_size = data->qae - data->qas;
    data->ident_perc = calc_ident_perc(data->qas, data->sas, as_size, &data->dist, &data->score);
    if (data->ident_perc < perc_identity) return;
   HbnTracebackDataDump(fprintf, stderr, tbck_data);

    s_fix_enzyme_align_ends(query_id, query_dir, fwd_query, rev_query, query_size,
        subject_id, subject, subject_size, perc_identity,
        b_qoff, b_qend, b_soff, b_send, tbck_data, reloci_list, qvep_list, cov_stats, align_list);
}

///////////////////////////////////////

