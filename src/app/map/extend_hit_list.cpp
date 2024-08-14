#include "extend_hit_list.hpp"
#include "../../corelib/pdqsort.h"

#include <algorithm>
#include <cmath>

using namespace std;

static const int kMatLen = 8;
static const int kExtendBlock = 50;
static const int kMaxExtendBlock = 90;
static const int kMaxW = 50;
static const double kEpsilon = 0.2;
static const int kEnzymeNbhd = 30;

static BOOL 
chain_seed_list_is_contained_in_align(const frag_align_struct* align_list,
    const int align_cnt,
    HbnInitHit* hit)
{
    int qbeg = hit->qbeg;
    int qend = hit->qend;
    int sbeg = hit->sbeg;
    int send = hit->send;
    for (int i = 0; i < align_cnt; ++i) {
        const frag_align_struct* ai = align_list + i;
        if (ai->qdir != hit->qdir) continue;
        if (ai->sid != hit->sid) continue;
        const int E = 10;
        int r = (qbeg + E >= ai->qoff)
                &&
                (qend <= ai->qend + E)
                &&
                (sbeg + E >= ai->soff)
                &&
                (send <= ai->send + E);
        if (r) {
            return TRUE;
        }
    }

    return FALSE;
}

static BOOL
chain_seed_list_is_contained_in_pca(PoreCAlign* pca_a,
    const int pca_c,
    HbnInitHit* hit)
{
    int qbeg = hit->qbeg;
    int qend = hit->qend;
    int sbeg = hit->sbeg;
    int send = hit->send;
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        if (pca->qdir != hit->qdir) continue;
        if (pca->sid != hit->sid) continue;
        const int E = 10;
        int r = (qbeg + E >= pca->qoff)
                &&
                (qend <= pca->qend + E)
                &&
                (sbeg + E >= pca->soff)
                &&
                (send <= pca->send + E);
        if (r) {
            return TRUE;
        }
    }

    return FALSE;
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
    string& qaln,
    string& saln)
{
    *_qend = 0;
    *_send = 0;
    qaln.clear();
    saln.clear();
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
                &ident_perc);
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
    string& qabuf,
    string& sabuf,
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
            qblk = (subject_length - st) + 50;;
            qblk = min<int>(qblk, query_length - qt);
            sblk = (query_length - qt) + 50;
            sblk = min<int>(sblk, subject_length - st);
            last_block = true;
        } else {
            qblk = kExtendBlock;
            sblk = kExtendBlock;
            last_block = false;
        }
        if (qblk == 0 || sblk == 0) break;
    //HBN_LOG("1 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, query_length, subject_length);
        int qfae = 0, sfae = 0;
        //edlib_shw(edlib, query + qt, qblk, subject + st, sblk, &qfae, &sfae, qabuf, sabuf);
        qabuf.clear();
        sabuf.clear();
        overhang_extend(dalign, edlib, query + qt, qblk, subject + st, sblk, &qfae, &sfae, qabuf, sabuf);
        if (qblk - qfae > 30 || sblk - sfae > 30) done = true;
        int acnt = 0, qcnt = 0, scnt = 0;
        int as_size = qabuf.size(), k = as_size - 1, m = 0;
        while (k >= 0) {
            char qc = qabuf[k];
            char sc = sabuf[k];
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
            memcpy(*qae, qabuf.c_str(), as_size);
            memcpy(*sae, sabuf.c_str(), as_size);
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
    string& qabuf,
    string& sabuf,
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
            qblk = min<int>(sf + 50, qf);
            sblk = min<int>(qf + 50, sf);
            last_block = true;
        } else {
            qblk = kExtendBlock;
            sblk = kExtendBlock;
            last_block = false;
        }
        if (qblk == 0 || sblk == 0) break;
        qfrag.clear(); for (int i = 1; i <= qblk; ++i) qfrag.push_back(query[qf-i]);
        sfrag.clear(); for (int i = 1; i <= sblk; ++i) sfrag.push_back(subject[sf-i]);
        //fprintf(stderr, "1 qf = %d, sf = %d, qblk = %d, sblk = %d\n", qf, sf, qblk, sblk);

        int qfae = 0, sfae = 0;
        qabuf.clear();
        sabuf.clear();
        edlib_shw(edlib, qfrag.data(), qblk, sfrag.data(), sblk, &qfae, &sfae, qabuf, sabuf);
        if (qblk - qfae > 30 || sblk - sfae > 30) done = true;
        int acnt = 0, qcnt = 0, scnt = 0;
        int as_size = qabuf.size(), k = as_size - 1, m = 0;
        while (k >= 0) {
            char qc = qabuf[k];
            char sc = sabuf[k];
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
                int qc = qabuf[i];
                --(*qas); **qas = qc;
                int sc = sabuf[i];
                --(*sas); **sas = sc;
            }
        }
        //HBN_LOG("sf = %d, qf = %d", qf, sf);
        if (done) break;
    }

    *qoff = qf;
    *soff = sf;
}

int
extend_one_chain(const int query_id,
    const int query_dir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    HbnInitHit* hit,
    const double perc_identity,
    HbnTracebackData* tbck_data,
    frag_align_list_struct& align_list)
{
    const u8* query = (query_dir == FWD) ? fwd_query : rev_query;
    int r = porec_compute_traceback(tbck_data, hit->qbeg, hit->qend, hit->sbeg, hit->send, query, query_size, subject, subject_size, 20, perc_identity);
    if (!r) return 0;
    //HbnTracebackDataDump(fprintf, stderr, tbck_data);
    int b_qoff = tbck_data->qoff;
    int b_qend = tbck_data->qend;
    int b_soff = tbck_data->soff;
    int b_send = tbck_data->send;
    hbn_assert(b_qend <= query_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    hbn_assert(b_send <= subject_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    HbnTracebackData* data = tbck_data;
    int as_size = 0;

    right_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qend, query_size,
        subject, &tbck_data->send, subject_size, tbck_data->ext_qabuf, tbck_data->ext_sabuf, &tbck_data->qae, &tbck_data->sae);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qend >= b_qend, "qid = %d, qend = %d, b_qend = %d, send = %d, b_send = %d", query_id, 
        tbck_data->qend, b_qend, tbck_data->send, b_send);
    hbn_assert(tbck_data->send >= b_send);

    left_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qoff, subject, &tbck_data->soff,
        tbck_data->ext_qabuf, tbck_data->ext_sabuf, &tbck_data->qas, &tbck_data->sas);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qoff <= b_qoff);
    hbn_assert(tbck_data->soff <= b_soff);
        
    as_size = data->qae - data->qas;
    data->ident_perc = calc_ident_perc(data->qas, data->sas, as_size, &data->dist, &data->score);
    if (data->ident_perc < perc_identity) return 0;
   //HbnTracebackDataDump(fprintf, stderr, tbck_data);

   frag_align_struct frag_align;
   frag_align.is_valid = 1;

   frag_align.qdir = query_dir;
   frag_align.qoff = tbck_data->qoff;
   frag_align.qend = tbck_data->qend;
   frag_align.qsize = query_size;
   frag_align.sid = subject_id;
   frag_align.soff = tbck_data->soff;
   frag_align.send = tbck_data->send;
   frag_align.ssize = subject_size;

   frag_align.score = tbck_data->score;
   frag_align.chain_score = hit->score;
   frag_align.pi = tbck_data->ident_perc;
   frag_align.qas_offset = align_list.align_strings.size();
   align_list.align_strings.append(tbck_data->qas, as_size);
   frag_align.sas_offset = align_list.align_strings.size();
   align_list.align_strings.append(tbck_data->sas, as_size);
   frag_align.as_size = as_size;

   frag_align.bqoff = b_qoff;
   frag_align.bqend = b_qend;
   frag_align.bsoff = b_soff;
   frag_align.bsend = b_send;

   frag_align.sqoff = hit->qbeg;
   frag_align.sqend = hit->qend;
   frag_align.ssoff = hit->sbeg;
   frag_align.ssend = hit->send;

   frag_align.left_align_id = -1;
   frag_align.right_align_id = -1;

   align_list.frag_align_list.push_back(frag_align);
   return 1;
}

///// fix enzyme ends

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

#if 1
    int x = edlib_nw_dist(tbck_data->edlib, q, ql, s, sl, dist, &pca->score);
    if (x == -1) return FALSE;

    double pi = 100.0 * (1.0 - 2.0 * x / (ql + sl));
    pca->pi = pi;
    return TRUE;
#else
    int qb, qe, sb, se;
    double pi;
int r = 
    nw_ksw2_extd2(tbck_data->ksw,
        0,
        q,
        0,
        ql,
        ql,
        0,
        s, 
        0,
        sl,
        sl,
        10,
        0.0,
        dist,
        &qb,
        &qe,
        &sb,
        &se,
        &pi,
        &tbck_data->qabuf,
        &tbck_data->sabuf);
    const char* qas = ks_s(tbck_data->qabuf);
    const char* sas = ks_s(tbck_data->sabuf);
    int as_size = ks_size(tbck_data->qabuf);
    pca->pi = calc_ident_perc(qas, sas, as_size, NULL, &pca->score);
    return TRUE;
#endif 
}

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
    pdqsort(offset_list.begin(), offset_list.end());
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
            ssubseq.data(), ssubseq.size(), &qe, &se, tbck_data->ext_qabuf, tbck_data->ext_sabuf);
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
            ssubseq.data(), ssubseq.size(), &qe, &se, tbck_data->ext_qabuf, tbck_data->ext_sabuf);
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
                                qvep_list->fwd_vdf_endpoint_list.data()
                                :
                                qvep_list->rev_vdf_endpoint_list.data();
    const int q_reloci_cnt = qvep_list->rev_vdf_endpoint_list.size();
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

    for (auto& qi : ql_list) {
        int q_d = b_qend - abs(qi);
        hbn_assert(q_d > 0);
        for (auto& si : sl_list) {
            int s_d = b_send - abs(si);
            hbn_assert(s_d > 0);
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            enzyme_offset eop;
            eop.qoff = abs(qi);
            eop.soff = abs(si);
            eop.enzyme_qoff = abs(qi);
            eop.enzyme_soff = abs(si);
            leop_list.push_back(eop);   
            qi = -abs(qi);
            si = -abs(si);        
        }
    }
    for (auto& qi : ql_list) {
        if (qi < 0) continue;
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
    for (auto& si : sl_list) {
        if (si < 0) continue;
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
        overhang_extend(tbck_data->dalign, tbck_data->edlib, query + qf, qt - qf, subject + sf, st - sf, &qe, &se, tbck_data->ext_qabuf, tbck_data->ext_sabuf);
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
        overhang_extend(tbck_data->dalign, tbck_data->edlib, query + qf, qt - qf, subject + sf, st - sf, &qe, &se, tbck_data->ext_qabuf, tbck_data->ext_sabuf);
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
                                qvep_list->fwd_vdf_endpoint_list.data()
                                :
                                qvep_list->rev_vdf_endpoint_list.data();
    const int q_reloci_cnt = qvep_list->rev_vdf_endpoint_list.size();
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

    for (auto& qi : qr_list) {
        int q_d = abs(qi) - b_qoff;
        //HBN_LOG("qi = %d, q_d = %d, b_qoff = %d", qi, q_d, b_qoff);
        hbn_assert(q_d > 0);
        for (auto& si : sr_list) {
            int s_d = abs(si) - b_soff;
            //HBN_LOG("\tsi = %d, s_d = %d, b_soff = %d", si, s_d, b_soff);
            hbn_assert(s_d > 0);
            int d_d = abs(q_d - s_d);
            if (d_d > kMaxW) continue;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > kEpsilon) continue;
            enzyme_offset eop;
            eop.qoff = abs(qi);
            eop.enzyme_qoff = abs(qi);
            eop.soff = abs(si);
            eop.enzyme_soff = abs(si);
            reop_list.push_back(eop);
            qi = -abs(qi);
            si = -abs(si);
        }
    } 
    for (auto& qi : qr_list) {
        if (qi < 0) continue;
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
    for (auto& si : sr_list) {
        if (si < 0) continue;
        enzyme_offset eop;
        eop.qoff = -1;
        eop.enzyme_qoff = ((qend - q_reloci_array[r_q_intv_i]) < (q_reloci_array[r_q_intv_i+1] - qend)) ? q_reloci_array[r_q_intv_i] : q_reloci_array[r_q_intv_i+1];
        eop.soff = si;
        eop.enzyme_soff = si;
        if (s_resolve_unfixed_end_offset(tbck_data, query, query_size, subject, subject_size, &eop)) {
            reop_list.push_back(eop);
        }        
    }

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
    const int chain_score,
    HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    std::vector<PoreCAlign>& align_list)
{
    const int* q_reloci_array = (query_dir == FWD)
                                ?
                                qvep_list->fwd_vdf_endpoint_list.data()
                                :
                                qvep_list->rev_vdf_endpoint_list.data();
    const int q_reloci_cnt = qvep_list->rev_vdf_endpoint_list.size();
    const int* s_reloci_array = reloci_list->reloci_array
                                +
                                reloci_list->seq_reloci_info_array[subject_id].enzyme_loci_offset;
    const int s_reloci_cnt = reloci_list->seq_reloci_info_array[subject_id].enzyme_loci_cnt;
    const u8* query = (query_dir == FWD) ? fwd_query : rev_query;

    vector<enzyme_offset> leop_list;
    s_add_start_offset_pairs(tbck_data, reloci_list, qvep_list, query, subject, query_dir, subject_id, tbck_data->qoff, tbck_data->soff, b_qend, b_send, leop_list);

    //fprintf(stderr, "==============\n");
    //HbnTracebackDataDump(fprintf, stderr, tbck_data);
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
            pca.chain_score = chain_score;
            align_list.push_back(pca);
        }
    }
}

static void
s_fix_enzyme_ends(HbnTracebackData* tbck_data,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    const int query_id,
    const int query_dir,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    int b_qoff,
    int b_qend,
    int b_soff,
    int b_send,
    int qoff,
    int qend,
    int soff,
    int send,
    int chain_score,
    const char* qas,
    const char* sas,
    const int as_size,
    vector<PoreCAlign>& pca_list)
{
    tbck_data->init(qoff, query_size, soff, subject_size);
    tbck_data->qoff = qoff;
    tbck_data->qend = qend;
    tbck_data->soff = soff;
    tbck_data->send = send;
    tbck_data->qas = tbck_data->qabuf.data() + min<int>(qoff, soff) * 2 + 200;
    tbck_data->sas = tbck_data->sabuf.data() + min<int>(qoff, soff) * 2 + 200;
    memcpy(tbck_data->qas, qas, as_size);
    tbck_data->qae = tbck_data->qas + as_size;
    memcpy(tbck_data->sas, sas, as_size);
    tbck_data->sae = tbck_data->sas + as_size;

    s_fix_enzyme_align_ends(query_id,
        query_dir,
        fwd_query,
        rev_query,
        query_size,
        subject_id,
        subject,
        subject_size,
        opts->perc_identity,
        b_qoff,
        b_qend,
        b_soff,
        b_send,
        chain_score,
        tbck_data,
        reloci_list,
        qvep_list,
        pca_list);
}

void
extend_hit_list(HbnTracebackData* tbck_data,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    HbnUnpackedDatabase* subjects,
    const int query_id,
    const char* query_name,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    HbnInitHit* hita,
    int hitc,
    const int max_hsps,
    int* cov_stats,
    vector<PoreCAlign>& pca_list)
{
#if 0
    HBN_LOG("mapping %d hits", hitc);
    vector<HbnInitHit> hit_list;
    for (int i = 0; i < hitc; ++i) {
        HbnInitHit hit = hita[i];
        if (hit.qdir == REV) {
            int x = query_size - hit.qend;
            int y = query_size - hit.qbeg;
            hit.qbeg = x;
            hit.qend = y;
        }
        hit_list.push_back(hit);
    }
    pdqsort(hit_list.begin(), hit_list.end(), [](const HbnInitHit& a, const HbnInitHit& b) { return (a.sid < b.sid) || (a.sid == b.sid && a.qbeg < b.qbeg); });
    for (int i = 0; i < hitc; ++i) {
        HbnInitHit& hit = hit_list[i];
        const char* sname = subjects->SeqName(hit.sid);
        fprintf(stderr, "%d\t%d:%s\t[%d, %d, %d] x [%d, %d], %d, qsize = %d\n", 
            i, hit.sid, sname, hit.qdir, hit.qbeg, hit.qend, hit.sbeg, hit.send, hit.score, hit.qsize);
    }       
#endif 

    const int verbose = 0;
    frag_align_list_struct align_list;
    pdqsort(hita, hita + hitc, [](const HbnInitHit& a, const HbnInitHit& b) { return a.score > b.score; });
    fill(cov_stats, cov_stats + query_size, 0);
    int added_align = 0;
    for (int i = 0; i < hitc && added_align < 50; ++i) {
        HbnInitHit* hit = hita + i;
        hita[0].qid = -1;
        const u8* subject = subjects->GetSequence(hit->sid); 
        const int subject_size = subjects->SeqSize(hit->sid); 
        const char* subject_name = subjects->SeqName(hit->sid);
        if (chain_seed_list_is_contained_in_align(
                align_list.frag_align_list.data(), 
                align_list.frag_align_list.size(),
                hit)) continue;
        if (verbose) {
            fprintf(stderr, "%d\t%d:%s\t[%d, %d] x [%d, %d], %d, qsize = %d\n", 
                i, hita[i].sid, subject_name, hita[i].qbeg, hita[i].qend, 
                hita[i].sbeg, hita[i].send, hita[i].score, hita[i].qsize);
        }
        int r = 
        extend_one_chain(query_id,
            hit->qdir,
            fwd_query,
            rev_query,
            query_size,
            hit->sid,
            subject,
            subject_size,
            hit,
            opts->perc_identity,
            tbck_data,
            align_list);
        if (r && verbose) dump_frag_align(fprintf, stderr, align_list.frag_align_list.back(), -1);
        if (!r) continue;
        int qoff = align_list.frag_align_list.back().qoff;
        int qend = align_list.frag_align_list.back().qend;
        if (hita[i].qdir == REV) {
            int x = query_size - qend;
            int y = query_size - qoff;
            qoff = x;
            qend = y;
        }
        for (int k = qoff; k < qend; ++k) ++cov_stats[k];
        ++added_align;
    }

    for (int i = 0; i < hitc; ++i) {
        HbnInitHit* hit = hita + i;
        if (hit->qid == -1) continue;
        const u8* subject = subjects->GetSequence(hit->sid);
        const int subject_size = subjects->SeqSize(hit->sid); 
        const char* subject_name = subjects->SeqName(hit->sid);
        int qoff = hit->qbeg;
        int qend = hit->qend;
        if (hita[i].qdir == REV) {
            int x = query_size - qend;
            int y = query_size - qoff;
            qoff = x;
            qend = y;            
        }
        int cov = 0;
        for (int k = qoff; k < qend; ++k) if (cov_stats[k]) ++cov;
        if (hita[i].qend - hita[i].qbeg - cov < 2) continue;
        if (verbose) {
            fprintf(stderr, "%d\t%d:%s\t[%d, %d] x [%d, %d], %d, qsize = %d\n", 
                i, hita[i].sid, subject_name, hita[i].qbeg, hita[i].qend, 
                hita[i].sbeg, hita[i].send, hita[i].score, hita[i].qsize);
        }
        int r = 
        extend_one_chain(query_id,
            hit->qdir,
            fwd_query,
            rev_query,
            query_size,
            hit->sid,
            subject,
            subject_size,
            hit,
            opts->perc_identity,
            tbck_data,
            align_list);
        if (r && verbose) dump_frag_align(fprintf, stderr, align_list.frag_align_list.back(), -1);
        if (!r) continue;
        qoff = align_list.frag_align_list.back().qoff;
        qend = align_list.frag_align_list.back().qend;
        if (hita[i].qdir == REV) {
            int x = query_size - qend;
            int y = query_size - qoff;
            qoff = x;
            qend = y;
        }
        for (int k = qoff; k < qend; ++k) ++cov_stats[k];        
    }

    frag_align_struct* fasa = align_list.frag_align_list.data();
    int fasc = align_list.frag_align_list.size();

#if 1
    pdqsort(fasa, fasa + fasc, [](const frag_align_struct& a, const frag_align_struct& b) { return a.ssoff < b.ssoff; });
    int score_list[fasc]; for (int i = 0; i < fasc; ++i) score_list[i] = fasa[i].score;
    int pred_list[fasc]; fill(pred_list, pred_list + fasc, -1);
    int succ_list[fasc]; fill(succ_list, succ_list + fasc, 0);
    int avail_list[fasc]; fill(avail_list, avail_list + fasc, 1);

    for (int i = 1; i < fasc; ++i) {
        frag_align_struct* fi = fasa + i;
        int max_score = fasa[i].score;
        int max_j = -1;
        for (int j = i - 1; j >= 0; --j) {
            frag_align_struct* fj = fasa + j;
            if (fi->qdir != fj->qdir || fi->sid != fj->sid) continue;
            if (fj->sqend > fi->sqoff || fj->ssend > fi->ssoff) continue;
	        if (fi->sqoff - fj->sqend > 300 || fi->ssoff - fj->ssend > 300) continue;
            int q_d = fi->sqend - fj->sqoff;
            int s_d = fi->ssend - fj->ssoff;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > 0.2) continue;
            int score = score_list[j] + fasa[i].score;
            if (score > max_score) {
                max_score = score;
                max_j = j;
            }
        }
        score_list[i] = max_score;
        pred_list[i] = max_j;
    }

    for (int i = 0; i < fasc; ++i) {
        if (pred_list[i] >= 0) succ_list[pred_list[i]] = 1;
    }
    vector<pair<int, int>> idx_and_score_list;
    for (int i = 0; i < fasc; ++i) {
        if (succ_list[i] == 0 && score_list[i] > 0) {
            idx_and_score_list.emplace_back(i, score_list[i]);
        }
    }
    pdqsort(idx_and_score_list.begin(), 
         idx_and_score_list.end(),
         [](const pair<int, int>& a, const pair<int, int>& b)->bool {
             return (a.second > b.second) || (a.second == b.second && a.first < b.first);
         });

    vector<frag_align_struct> colin_align_list;
    for (auto& ias : idx_and_score_list) {
        if (!avail_list[ias.first]) continue;
        int p = ias.first;
        colin_align_list.clear();
        while (p >= 0) {
            avail_list[p] = 0;
            colin_align_list.push_back(align_list.frag_align_list[p]);
            p = pred_list[p];
        }
        reverse(colin_align_list.begin(), colin_align_list.end());
        if (colin_align_list.size() < 2) continue;
        int sid = colin_align_list.front().sid;
        int qdir = colin_align_list.front().qdir;
        const u8* subject = subjects->GetSequence(sid);
        const int subject_size = subjects->SeqSize(sid);
        HbnInitHit hit;
        hit.qid = query_id;
        hit.qdir = qdir;
        hit.qbeg = colin_align_list.front().sqoff;
        hit.qend = colin_align_list.back().sqend;
        hit.qsize = query_size;
        hit.sbeg = colin_align_list.front().ssoff;
        hit.send = colin_align_list.back().ssend;
        hit.ssize = subject_size;
        hit.score = 0;
        for (auto& f : colin_align_list) hit.score += f.chain_score;
        extend_one_chain(query_id,
            qdir,
            fwd_query,
            rev_query,
            query_size,
            sid,
            subject,
            subject_size,
            &hit,
            opts->perc_identity,
            tbck_data,
            align_list);
    }

    fasa = align_list.frag_align_list.data();
    fasc = align_list.frag_align_list.size();
    pdqsort(fasa, fasa + fasc, [](const frag_align_struct& a, const frag_align_struct& b) { return a.score > b.score; });
    //for (int i = 0; i < fasc; ++i) dump_frag_align(fprintf, stderr, fasa[i], i);
    for (int i = 0; i < fasc; ++i) {
        frag_align_struct* fi = fasa + i;
        if (!fi->is_valid) continue;
        for (int j = i + 1; j < fasc; ++j) {
            frag_align_struct* fj = fasa + j;
            if (!fj->is_valid) continue;
            if (fi->qdir != fj->qdir) continue;
            if (fi->sid != fj->sid) continue;
            const int E = 10;
            int r = (fj->qoff + E >= fi->qoff) && (fj->qend <= fi->qend + E) && (fj->soff + E >= fi->soff) && (fj->send <= fi->send + E);
            if (r) fj->is_valid = 0;
        }
    }
    int n = 0;
    for (int i = 0; i < fasc; ++i) if (fasa[i].is_valid) fasa[n++] = fasa[i];
    fasc = n;
    //sort(fasa, fasa + fasc, [](const frag_align_struct& a, const frag_align_struct& b) { return a.qoff < b.qoff; });
    //for (int i = 0; i < fasc; ++i) dump_frag_align(fprintf, stderr, fasa[i], i);
#endif

    for (int i = 0; i < fasc; ++i) {
        const char* qas = align_list.align_strings.c_str() + fasa[i].qas_offset;
        const char* sas = align_list.align_strings.c_str() + fasa[i].sas_offset;
        int as_size = fasa[i].as_size;
        s_fix_enzyme_ends(tbck_data,
            opts,
            reloci_list,
            qvep_list,
            query_id,
            fasa[i].qdir,
            query_name,
            fwd_query,
            rev_query,
            query_size,
            fasa[i].sid,
            subjects->GetSequence(fasa[i].sid),
            subjects->SeqSize(fasa[i].sid),
            fasa[i].bqoff,
            fasa[i].bqend,
            fasa[i].bsoff,
            fasa[i].bsend,
            fasa[i].qoff,
            fasa[i].qend,
            fasa[i].soff,
            fasa[i].send,
            fasa[i].chain_score,
            qas,
            sas,
            as_size,
            pca_list);
    }
    //for (auto& pca : pca_list) dump_chain_pca(fprintf, stderr, pca, -1);
}

void
align_subseq(const int query_id,
    const int query_dir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    const int qb, 
    const int qe,
    const int sb,
    const int se,
    const int chain_score,
    const double perc_identity,
    HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    std::vector<PoreCAlign>& align_list)
{
    const u8* query = (query_dir == FWD) ? fwd_query : rev_query;
    int r = porec_compute_traceback(tbck_data, qb, qe, sb, se, query, query_size, subject, subject_size, 20, perc_identity);
    if (!r) return;
    //HbnTracebackDataDump(fprintf, stderr, tbck_data);
    int b_qoff = tbck_data->qoff;
    int b_qend = tbck_data->qend;
    int b_soff = tbck_data->soff;
    int b_send = tbck_data->send;
    hbn_assert(b_qend <= query_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    hbn_assert(b_send <= subject_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    HbnTracebackData* data = tbck_data;
    int as_size = 0;

    right_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qend, query_size,
        subject, &tbck_data->send, subject_size, tbck_data->ext_qabuf, tbck_data->ext_sabuf, &tbck_data->qae, &tbck_data->sae);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qend >= b_qend, "qid = %d, qend = %d, b_qend = %d, send = %d, b_send = %d", query_id, 
        tbck_data->qend, b_qend, tbck_data->send, b_send);
    hbn_assert(tbck_data->send >= b_send);

    left_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qoff, subject, &tbck_data->soff,
        tbck_data->ext_qabuf, tbck_data->ext_sabuf, &tbck_data->qas, &tbck_data->sas);
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
   //HbnTracebackDataDump(fprintf, stderr, tbck_data);

    s_fix_enzyme_align_ends(query_id, query_dir, fwd_query, rev_query, query_size,
        subject_id, subject, subject_size, perc_identity,
        b_qoff, b_qend, b_soff, b_send, chain_score, tbck_data, reloci_list, qvep_list, align_list);
}

bool
align_subseq_enzyme_inference(const int query_id,
    const int query_dir,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    const int subject_id,
    const u8* subject,
    const int subject_size,
    const int qb, 
    const int qe,
    const int sb,
    const int se,
    const double perc_identity,
    HbnTracebackData* tbck_data)
{
    const u8* query = (query_dir == FWD) ? fwd_query : rev_query;
    int r = porec_compute_traceback(tbck_data, qb, qe, sb, se, query, query_size, subject, subject_size, 20, perc_identity);
    if (!r) return false;
    //HbnTracebackDataDump(fprintf, stderr, tbck_data);
    int b_qoff = tbck_data->qoff;
    int b_qend = tbck_data->qend;
    int b_soff = tbck_data->soff;
    int b_send = tbck_data->send;
    hbn_assert(b_qend <= query_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    hbn_assert(b_send <= subject_size, "b_qend = %d, qsize = %d, b_send = %d, ssize = %d", b_qend, query_size, b_send, subject_size);
    HbnTracebackData* data = tbck_data;
    int as_size = 0;

    right_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qend, query_size,
        subject, &tbck_data->send, subject_size, tbck_data->ext_qabuf, tbck_data->ext_sabuf, &tbck_data->qae, &tbck_data->sae);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qend >= b_qend, "qid = %d, qend = %d, b_qend = %d, send = %d, b_send = %d", query_id, 
        tbck_data->qend, b_qend, tbck_data->send, b_send);
    hbn_assert(tbck_data->send >= b_send);

    left_extend_align(tbck_data->dalign, tbck_data->edlib, query, &tbck_data->qoff, subject, &tbck_data->soff,
        tbck_data->ext_qabuf, tbck_data->ext_sabuf, &tbck_data->qas, &tbck_data->sas);
    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
    hbn_assert(tbck_data->qoff <= b_qoff);
    hbn_assert(tbck_data->soff <= b_soff);
        
    as_size = data->qae - data->qas;
    data->ident_perc = calc_ident_perc(data->qas, data->sas, as_size, &data->dist, &data->score);
    if (data->ident_perc < perc_identity) return false;
   //HbnTracebackDataDump(fprintf, stderr, tbck_data);

    return true;
}