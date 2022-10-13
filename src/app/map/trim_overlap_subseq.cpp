#include "trim_overlap_subseq.hpp"

#include <algorithm>
#include <string>
#include <vector>

#include "../../algo/align_format/format_aux.h"
#include "../../algo/align_format/paf.h"
#include "../../algo/hbn_traceback_aux.h"
#include "../../corelib/ksort.h"

using namespace std;

static const int kMinFragSize = 20;

typedef struct {
    int is_valid;
    PoreCAlign pca;
    int qas_offset;
    int sas_offset;
    int as_size;
    int qdir;
    int qoff, qend;
    int e_qoff, e_qend;
    int sdir;
    int soff, send;
    int e_soff, e_send;
} TrimPca;

#define dump_tpca(output_func, stream, tpca) output_func(stream, "[%d, %d, %d, %d] x [%d, %d, %d, %d]\n", \
	tpca.pca.qid, tpca.qdir, tpca.qoff, tpca.qend,\
	tpca.pca.sid, tpca.sdir, tpca.soff, tpca.send)

bool operator < (const TrimPca& lhs, const TrimPca& rhs)
{
	if (lhs.qoff >= rhs.qoff && lhs.qend <= rhs.qend) return true;
	if (lhs.pca.sid == rhs.pca.sid && lhs.soff >= rhs.soff && lhs.send <= rhs.send) return true;
	return false;
}

struct TrimPcaList
{
    vector<TrimPca> tpca_list;
    string align_strings;

    void add(PoreCAlign* pca, 
            int qoff, int qend, int soff, int send,
            const char* qas, const char* sas, const int as_size) {
        TrimPca tpca;
        tpca.is_valid = 1;
        tpca.pca = *pca;
        tpca.qas_offset = align_strings.size();
        align_strings.append(qas, as_size);
        tpca.sas_offset = align_strings.size();
        align_strings.append(sas, as_size);
        tpca.as_size = as_size;

        tpca.qdir = pca->qdir;
        tpca.qoff = qoff;
        tpca.qend = qend;
        tpca.e_qoff = pca->enzyme_qoff;
        tpca.e_qend = pca->enzyme_qend;
        tpca.sdir = FWD;
        tpca.soff = soff;
        tpca.send = send;
        tpca.e_soff = pca->enzyme_soff;
        tpca.e_send = pca->enzyme_send;
        tpca_list.push_back(tpca);
    }
};

void
set_trim_pca_list_for_fwd_query(TrimPcaList* src, TrimPcaList* dst, const int enzyme_size)
{
    dst->align_strings.clear();
    dst->tpca_list.clear();
    TrimPca dpca;
    for (auto& tpca : src->tpca_list) {
        if (!tpca.is_valid) continue;
        if (tpca.qdir == FWD) {
            dpca = tpca;
            const char* qas = src->align_strings.c_str() + tpca.qas_offset;
            const char* sas = src->align_strings.c_str() + tpca.sas_offset;
            const int as_size = tpca.as_size;
            dpca.qas_offset = dst->align_strings.size();
            dst->align_strings.append(qas, as_size);
            dpca.sas_offset = dst->align_strings.size();
            dst->align_strings.append(sas, as_size);
            dpca.as_size = as_size;
            dst->tpca_list.push_back(dpca);
            continue;
        }

        dpca.pca = tpca.pca;
        dpca.is_valid = 1;

        dpca.qdir = FWD;
        dpca.qoff = tpca.pca.qsize - tpca.qend;
        dpca.qend = tpca.pca.qsize - tpca.qoff;
        dpca.e_qoff = 
        (tpca.e_qend == tpca.pca.qsize) ? 0 : (tpca.pca.qsize - tpca.e_qend - enzyme_size);
        dpca.e_qend = 
        (tpca.e_qoff == 0) ? tpca.pca.qsize : (tpca.pca.qsize - tpca.e_qoff - enzyme_size);
        
        dpca.sdir = 1 - tpca.sdir;
        dpca.soff = tpca.pca.ssize - tpca.send;
        dpca.send = tpca.pca.ssize - tpca.soff;
        dpca.e_soff = 
        (tpca.e_send == tpca.pca.ssize) ? 0 : (tpca.pca.ssize - tpca.e_send - enzyme_size);
        dpca.e_send = 
        (tpca.e_soff == 0) ? tpca.pca.ssize : (tpca.pca.ssize - tpca.e_soff - enzyme_size);

        const char* qas = src->align_strings.c_str() + tpca.qas_offset;
        const char* sas = src->align_strings.c_str() + tpca.sas_offset;
        const int as_size = tpca.as_size;
        int as_i;

        dpca.qas_offset = dst->align_strings.size();
        as_i = as_size;
        while (as_i) {
            --as_i;
            int c = qas[as_i];
            if (c != GAP_CHAR) {
                c = 3 - nst_nt4_table[c];
                c = DECODE_RESIDUE(c);
            }
            dst->align_strings += c;
        }

        dpca.sas_offset = dst->align_strings.size();
        as_i = as_size;
        while (as_i) {
            --as_i;
            int c = sas[as_i];
            if (c != GAP_CHAR) {
                c = 3 - nst_nt4_table[c];
                c = DECODE_RESIDUE(c);
            }
            dst->align_strings += c;
        }
        dpca.as_size = as_size;
        dst->tpca_list.push_back(dpca);
    }
}

void
set_trim_pca_list_for_fwd_subject(TrimPcaList* src, TrimPcaList* dst, const int enzyme_size)
{
    dst->align_strings.clear();
    dst->tpca_list.clear();
    TrimPca dpca;
    for (auto& tpca : src->tpca_list) {
        if (!tpca.is_valid) {
            //HBN_LOG("invalid pca");
            //dump_pca(fprintf, stderr, tpca.pca, -1);
            continue;
        }
        if (tpca.sdir == FWD) {
            dpca = tpca;
            const char* qas = src->align_strings.c_str() + tpca.qas_offset;
            const char* sas = src->align_strings.c_str() + tpca.sas_offset;
            const int as_size = tpca.as_size;
            dpca.qas_offset = dst->align_strings.size();
            dst->align_strings.append(qas, as_size);
            dpca.sas_offset = dst->align_strings.size();
            dst->align_strings.append(sas, as_size);
            dpca.as_size = as_size;
            dst->tpca_list.push_back(dpca);
            continue;
        }

        dpca.pca = tpca.pca;
        dpca.is_valid = 1;

        dpca.sdir = FWD;
        dpca.soff = tpca.pca.ssize - tpca.send;
        dpca.send = tpca.pca.ssize - tpca.soff;
        dpca.e_soff = 
        (tpca.e_send == tpca.pca.ssize) ? 0 : (tpca.pca.ssize - tpca.e_send - enzyme_size);
        dpca.e_send = 
        (tpca.e_soff == 0) ? tpca.pca.ssize : (tpca.pca.ssize - tpca.e_soff - enzyme_size);

        dpca.qdir = 1 - tpca.qdir;
        dpca.qoff = tpca.pca.qsize - tpca.qend;
        dpca.qend = tpca.pca.qsize - tpca.qoff;
        dpca.e_qoff = 
        (tpca.e_qend == tpca.pca.qsize) ? 0 : (tpca.pca.qsize - tpca.e_qend - enzyme_size);
        dpca.e_qend = 
        (tpca.e_qoff == 0) ? tpca.pca.qsize : (tpca.pca.qsize - tpca.e_qoff - enzyme_size);

        const char* qas = src->align_strings.c_str() + tpca.qas_offset;
        const char* sas = src->align_strings.c_str() + tpca.sas_offset;
        int as_size = tpca.as_size;
        int as_i;

        dpca.qas_offset = dst->align_strings.size();
        as_i = as_size;
        while (as_i) {
            --as_i;
            int c = qas[as_i];
            if (c != GAP_CHAR) {
                c = 3 - nst_nt4_table[c];
                c = DECODE_RESIDUE(c);
            }
            dst->align_strings += c;
        }

        dpca.sas_offset = dst->align_strings.size();
        as_i = as_size;
        while (as_i) {
            --as_i;
            int c = sas[as_i];
            if (c != GAP_CHAR) {
                c = 3 - nst_nt4_table[c];
                c = DECODE_RESIDUE(c);
            }
            dst->align_strings += c;
        }
        dpca.as_size = as_size;
        dst->tpca_list.push_back(dpca);
    }
}

static void
s_divide_enzyme_offset(RestrictEnzyme* enzyme,
    PoreCAlign* pca,
    int* qb_,
    int* qe_,
    int* sb_,
    int* se_)
{
    if (pca->qdir == FWD) {
        int qb = pca->qoff;
        int qe = pca->qend;
        int sb = pca->soff;
        int se = pca->send;
        if (pca->enzyme_qoff > 0 && pca->qoff == pca->enzyme_qoff && pca->soff == pca->enzyme_soff) {
            qb += enzyme->break_loci;
            sb += enzyme->break_loci;
        }
        if (pca->enzyme_qend != pca->qsize && pca->qend == pca->enzyme_qend && pca->send == pca->enzyme_send) {
            qe += enzyme->break_loci;
            se += enzyme->break_loci;
        }
        *qb_ = qb;
        *qe_ = qe;
        *sb_ = sb;
        *se_ = se;
    } else {
        int qb = pca->qsize - pca->qend;
        int qe = pca->qsize - pca->qoff;
        int sb = pca->soff;
        int se = pca->send;
        if (pca->enzyme_qoff > 0 && pca->qoff == pca->enzyme_qoff && pca->soff == pca->enzyme_soff) {
            qe -= enzyme->enzyme_size;
            sb += enzyme->enzyme_size;
            qe += enzyme->break_loci;
            sb -= enzyme->break_loci;
        }        
        if (pca->enzyme_qend != pca->qsize && pca->qend == pca->enzyme_qend && pca->send == pca->enzyme_send) {
            qb -= enzyme->enzyme_size;
            se += enzyme->enzyme_size;
            qb += enzyme->break_loci;
            se -= enzyme->break_loci;
        } 
        *qb_ = pca->qsize - qe;
        *qe_ = pca->qsize - qb;
        *sb_ = sb;
        *se_ = se;       
    }
}

void
get_pca_align_string_edlib(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    SeqReader* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    PoreCAlign* pca,
    TrimPcaList* tpca_list)
{
    bool rr = (pca->qoff >=0 && pca->qoff < pca->qend && pca->qend <= pca->qsize)
             &&
             (pca->soff >= 0 && pca->soff < pca->send && pca->send <= pca->ssize);
    if (!rr) {
        fprintf(stderr, "%s\n", query_name);
        dump_pca(fprintf, stderr, *pca, -1);
        abort();
    }
    //HBN_LOG("get align string for");
    //dump_pca(fprintf, stderr, *pca, -1);
    const u8* query = (pca->qdir == FWD) ? fwd_query : rev_query;
    const u8* subject = SeqReader_Seq(subjects, pca->sid, FWD);
    const int subject_size = SeqReader_SeqSize(subjects, pca->sid);
    int qb, qe, sb, se;
    s_divide_enzyme_offset(&reloci_list->enzyme, pca, &qb, &qe, &sb, &se);
    int dist = max(qe - qb, se - sb);
    dist = dist * 0.5;// (100.0 - pca->pi) * 1.2 / 100.0;
    
    int a_qb = 0, a_qe = 0, a_sb = 0, a_se = 0;
    double a_pi = 0.0, a_epi = 0.0;
    const char* qas = NULL;
    const char* qae = NULL;
    const char* sas = NULL;
    const char* sae = NULL;
    int as_size = 0;

    int r = edlib_nw(tbck_data->edlib,
        query + qb,
        qe - qb,
        subject + sb,
        se - sb,
        dist,
        &tbck_data->qabuf,
        &tbck_data->sabuf);
    if (!r) return;
    a_qb = qb;
    a_qe = qe;
    a_sb = sb;
    a_se = se;

    qas = ks_s(tbck_data->qabuf);
    sas = ks_s(tbck_data->sabuf);
    as_size = ks_size(tbck_data->qabuf);
    qae = qas + as_size;
    sae = sas + as_size;  
    validate_aligned_string(__FILE__, __func__, __LINE__,
        0, query, a_qb, a_qe, qas,
        0, subject, a_sb, a_se, sas, as_size, TRUE); 

    tpca_list->add(pca, a_qb, a_qe, a_sb, a_se, qas, sas, as_size);
    //fprintf(stderr, "[%d, %d] x [%d, %d]\n", a_qb, a_qe, a_sb, a_se);
}

void
get_pca_align_string(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    SeqReader* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    PoreCAlign* pca,
    TrimPcaList* tpca_list)
{
    bool rr = (pca->qoff >=0 && pca->qoff < pca->qend && pca->qend <= pca->qsize)
             &&
             (pca->soff >= 0 && pca->soff < pca->send && pca->send <= pca->ssize);
    if (!rr) {
        fprintf(stderr, "%s\n", query_name);
        dump_pca(fprintf, stderr, *pca, -1);
        abort();
    }
    //HBN_LOG("get align string for");
    //dump_pca(fprintf, stderr, *pca, -1);
    const u8* query = (pca->qdir == FWD) ? fwd_query : rev_query;
    const u8* subject = SeqReader_Seq(subjects, pca->sid, FWD);
    const int subject_size = SeqReader_SeqSize(subjects, pca->sid);
    int qb, qe, sb, se;
    s_divide_enzyme_offset(&reloci_list->enzyme, pca, &qb, &qe, &sb, &se);
    int dist = max(qe - qb, se - sb);
    dist = dist * (100.0 - pca->pi) * 1.2 / 100.0;
    if (dist == 0) dist = max(qe - qb, se - sb) * 0.2;
    
    int a_qb = 0, a_qe = 0, a_sb = 0, a_se = 0;
    double a_pi = 0.0, a_epi = 0.0;
    const char* qas = NULL;
    const char* qae = NULL;
    const char* sas = NULL;
    const char* sae = NULL;
    int as_size = 0;

    int r = 0;
    r = nw_ksw2_extd2(tbck_data->ksw,
            0,
            query,
            qb,
            qe,
            query_size,
            0,
            subject,
            sb,
            se,
            subject_size,
            1,
            0.0,
            dist,
            &a_qb,
            &a_qe,
            &a_sb,
            &a_se,
            &a_pi,
            &tbck_data->qabuf,
            &tbck_data->sabuf);
    if (r) {
        qas = ks_s(tbck_data->qabuf);
        sas = ks_s(tbck_data->sabuf);
        as_size = ks_size(tbck_data->qabuf);
        a_epi = calc_effective_ident_perc(qas, sas, as_size);
        //if (a_epi < pca->pi - 3.0) r = 0;
    }
    if (!r) {
        //HBN_LOG("redo ksw, a_epi = %g, pi = %g", a_epi, pca->pi);
        //dump_pca(fprintf, stderr, *pca, -1);
        dist = max(qe - qb, se - sb);
        r = nw_ksw2_extd2(tbck_data->ksw,
                0,
                query,
                qb,
                qe,
                query_size,
                0,
                subject,
                sb,
                se,
                subject_size,
                1,
                0.0,
                dist,
                &a_qb,
                &a_qe,
                &a_sb,
                &a_se,
                &a_pi,
                &tbck_data->qabuf,
                &tbck_data->sabuf);        
    }
    if (!r) return;

    qas = ks_s(tbck_data->qabuf);
    sas = ks_s(tbck_data->sabuf);
    as_size = ks_size(tbck_data->qabuf);
    qae = qas + as_size;
    sae = sas + as_size;  
    validate_aligned_string(__FILE__, __func__, __LINE__,
        0, query, a_qb, a_qe, qas,
        0, subject, a_sb, a_se, sas, as_size, TRUE); 

    tpca_list->add(pca, a_qb, a_qe, a_sb, a_se, qas, sas, as_size);
    //fprintf(stderr, "[%d, %d] x [%d, %d]\n", a_qb, a_qe, a_sb, a_se);
}

static void
s_trim_qend(TrimPcaList* tpca_list, TrimPca* tpca, int trim_qend_to)
{
    hbn_assert(tpca->qend > trim_qend_to);
    if (tpca->qoff >= trim_qend_to) {
        tpca->is_valid = 0;
        return;
    }
    int qend = tpca->qend;
    int send = tpca->send;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_i = tpca->as_size;
    while (as_i) {
        --as_i;
        if (qas[as_i] != GAP_CHAR) --qend;
        if (sas[as_i] != GAP_CHAR) --send;
        if (qend == trim_qend_to) break;
    }
    hbn_assert(qend == trim_qend_to);
    tpca->qend = qend;
    tpca->send = send;
    tpca->as_size = as_i;
}

static void
s_trim_qoff(TrimPcaList* tpca_list, TrimPca* tpca, int trim_qoff_to)
{
    hbn_assert(tpca->qoff < trim_qoff_to);
    if (tpca->qend <= trim_qoff_to) {
        tpca->is_valid = 0;
        return;
    }
    int qoff = tpca->qoff;
    int soff = tpca->soff;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    //HBN_LOG("before:");
    //dump_align_string(qas, sas, tpca->as_size, stderr);
    int as_i = 0;
    while (as_i < tpca->as_size) {
        if (qas[as_i] != GAP_CHAR) ++qoff;
        if (sas[as_i] != GAP_CHAR) ++soff;
        ++as_i;
        if (qoff == trim_qoff_to) break;
    }    
    hbn_assert(qoff == trim_qoff_to);
    hbn_assert(as_i <= tpca->as_size);
    tpca->qoff = qoff;
    tpca->soff = soff;
    tpca->qas_offset += as_i;
    tpca->sas_offset += as_i;
    tpca->as_size -= as_i;

    //HBN_LOG("after:");
    //qas += as_i;
    //sas += as_i;
    //dump_align_string(qas, sas, tpca->as_size, stderr);
}

static void
s_trim_send(TrimPcaList* tpca_list, TrimPca* tpca, int trim_send_to)
{
    hbn_assert(tpca->send > trim_send_to);
    if (tpca->soff >= trim_send_to) {
        tpca->is_valid = 0;
        return;
    }
    int qend = tpca->qend;
    int send = tpca->send;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_i = tpca->as_size;
    while (as_i) {
        --as_i;
        if (qas[as_i] != GAP_CHAR) --qend;
        if (sas[as_i] != GAP_CHAR) --send;
        if (send == trim_send_to) break;
    }
    hbn_assert(send == trim_send_to);
    tpca->qend = qend;
    tpca->send = send;
    tpca->as_size = as_i;
}

static void
s_trim_soff(TrimPcaList* tpca_list, TrimPca* tpca, int trim_soff_to)
{
    hbn_assert(tpca->soff < trim_soff_to);
    if (tpca->send <= trim_soff_to) {
        tpca->is_valid = 0;
        return;
    }
    int qoff = tpca->qoff;
    int soff = tpca->soff;
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_i = 0;
    while (as_i < tpca->as_size) {
        if (qas[as_i] != GAP_CHAR) ++qoff;
        if (sas[as_i] != GAP_CHAR) ++soff;
        ++as_i;
        if (soff == trim_soff_to) break;
    }    
    hbn_assert(as_i <= tpca->as_size);
    hbn_assert(soff == trim_soff_to);
    tpca->qoff = qoff;
    tpca->soff = soff;
    tpca->qas_offset += as_i;
    tpca->sas_offset += as_i;
    tpca->as_size -= as_i;
}

static void
s_trim_query_overlap_subseq(TrimPcaList* tpca_list)
{
    int n_tpca = tpca_list->tpca_list.size();
    for (int i = 0; i < n_tpca - 1; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        TrimPca* pj = &tpca_list->tpca_list[i+1];
        if (pi->is_valid == 0 || pj->is_valid == 0) continue;
        if (pi->qend <= pj->qoff) continue;
	if ((*pi) < (*pj)) {
		pi->is_valid = 0;
		break;
	}
	if ((*pj) < (*pi)) {
		pj->is_valid = 0;
		continue;
	}	

        int i_q_e = (pi->qend == pi->e_qend);
        int i_s_e = (pi->send == pi->e_send);
        int j_q_s = (pj->qoff == pj->e_qoff);
        int j_s_s = (pj->soff == pj->e_soff);

        if (i_q_e && i_s_e && j_q_s && j_s_s) {
            if (pi->qend - pi->qoff > pj->qend - pj->qoff) {
                s_trim_qoff(tpca_list, pj, pi->qend);
            } else {
                s_trim_qend(tpca_list, pi, pj->qoff);
            }
            continue;
        }

        if ((i_q_e || i_s_e) && (j_q_s == 0 && j_s_s == 0)) {
            s_trim_qoff(tpca_list, pj, pi->qend);
            continue;
        }

        if ((i_q_e == 0 && i_s_e == 0) && (j_q_s || j_s_s)) {
            s_trim_qend(tpca_list, pi, pj->qoff);
            continue;
        }

        if (pi->qend - pi->qoff > pj->qend - pj->qoff) {
            s_trim_qoff(tpca_list, pj, pi->qend);
        } else {
            s_trim_qend(tpca_list, pi, pj->qoff);
        }
    }

    for (int i = 0; i < n_tpca; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        if (!pi->is_valid) continue;
        hbn_assert(pi->qoff < pi->qend);
        if (pi->qend - pi->qoff < kMinFragSize) pi->is_valid = 0;
    }
}

static void
s_trim_subject_overlap_subseq(TrimPcaList* tpca_list)
{
    const int n_tpca = tpca_list->tpca_list.size();
    for (int i = 0; i < n_tpca; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        if (!pi->is_valid) continue;
        for (int j = i + 1; j < n_tpca; ++j) {
            TrimPca* pj = &tpca_list->tpca_list[j];
            if (!pj->is_valid) continue;
            if (pi->pca.sid != pj->pca.sid) continue;

            TrimPca* pl = (pi->soff < pj->soff) ? pi : pj;
            TrimPca* pr = (pi->soff < pj->soff) ? pj : pi;
            if (pl->send <= pr->soff) continue;

            int l_q_e = (pl->qend == pl->e_qend);
            int l_s_e = (pl->send == pl->e_send);
            int r_q_s = (pr->qoff == pr->e_qoff);
            int r_s_s = (pr->soff == pr->e_soff);

            if (l_q_e && l_s_e && r_q_s && r_s_s) {
                if (pl->send - pl->soff < pr->send - pr->soff) {
                    s_trim_send(tpca_list, pl, pr->soff);
                } else {
                    s_trim_soff(tpca_list, pr, pl->send);
                }
                continue;
            }

            if ((l_q_e || l_s_e) && (r_q_s == 0 && r_s_s)) {
                s_trim_soff(tpca_list, pr, pl->send);
                continue;
            }

            if ((l_q_e == 0 && l_s_e == 0) && (r_q_s || r_s_s)) {
                s_trim_send(tpca_list, pl, pr->soff);
                continue;
            }

            if (pl->send - pl->soff > pr->send - pr->soff) {
                s_trim_soff(tpca_list, pr, pl->send);
            } else {
                s_trim_send(tpca_list, pl, pr->soff);
            }
        }
    }    

    for (int i = 0; i < n_tpca; ++i) {
        TrimPca* pi = &tpca_list->tpca_list[i];
        if (!pi->is_valid) continue;
        if (pi->qend - pi->qoff < kMinFragSize) pi->is_valid = 0;
    }
}

static void
s_dump_repeat_pca(PoreCAlign* pca, const char* query_name, SeqReader* subjects, kstring_t* out)
{
    ksprintf(out, "hm:Z:");
    const char* subject_name = SeqReader_SeqName(subjects, pca->sid);
    const char strand = (pca->qdir == FWD) ? '+' : '-';
    ksprintf(out, "%s:%c:%d:%d:%s:%d:%d:%g", query_name, strand, pca->chain_qoff, pca->chain_qend,
        subject_name, pca->soff, pca->send, pca->pi);
}

static void
dump_sam_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    if (qoff) ksprintf(out, "%dS", qoff);
    dump_cigar_string(qaln, saln, aln_size, out);
    if (qend < qsize) ksprintf(out, "%dS", qsize - qend);
}

static BOOL 
validate_cigar_and_seq_size(const char* cigar, const int query_size)
{
    const int c_n = strlen(cigar);
    int c_i = 0;
    int q_base = 0;
    while (c_i < c_n) {
        size_t j = c_i + 1;
        while (j < c_n && isdigit(cigar[j])) ++j;
        char op = cigar[j];
        int num = atoi(cigar + c_i);
        switch (op)
        {
        case 'S':
            q_base += num;
            break;
        case 'M':
            q_base += num;
            break;
        case 'D':
            break;
        case 'I':
            q_base += num;
            break;
        default:
            break;
        }
        c_i = j + 1;
    }    
    if (q_base != query_size) {
        HBN_LOG("cigar and seq_size is inconsistent: %d v.s. %d", q_base, query_size);
        return false;
    }
    return true;
}

void
s_dump_sam_pca(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    SeqReader* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const char* query_qv,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    TrimPcaList* tpca_list,
    TrimPca* tpca,
    const EChainType chain_type,
    const int extended_subject,
    kstring_t* out)
{
    const char tab = '\t';
    int flag = 0;
    if (tpca->pca.qdir ==  REV) flag |= 0x10; // reverse query strand

    ksprintf(out, "%s", query_name); /// 1) query name --- string  [!-?A-~]{1,254}
    kputc(tab, out);
    ksprintf(out, "%d", flag); /// 2) flag --- int  [0, 2^16-1]
    kputc(tab, out);
    const char* subject_name = SeqReader_SeqName(subjects, tpca->pca.sid);
    ksprintf(out, "%s", subject_name); /// 3) subject name --- string  \*|[:rname:^*=][:rname:]*
    kputc(tab, out);
    ksprintf(out, "%d", tpca->soff + 1); /// 4) left most subject position (1-based) --- int  [0, 2^31-1]
    kputc(tab, out);
    ksprintf(out, "%d", tpca->pca.map_q); /// 5) mapq --- int [0, 2^8-1]
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    const int as_size = tpca->as_size;
    kputc(tab, out);
    dump_sam_cigar(tpca->qoff, tpca->qend, query_size, qas, sas, as_size, out); /// 6) cigar --- string  \*|([0-9]+[MIDNSHPX=])+
    ks_dinit(cigar);
    dump_sam_cigar(tpca->qoff, tpca->qend, query_size, qas, sas, as_size, &cigar);
    validate_cigar_and_seq_size(ks_s(cigar), query_size);
    ks_destroy(cigar);
    kputc(tab, out);
    ksprintf(out, "*"); /// 7) reference name of the mate/next read --- string  \*|=|[:rname:^*=][:rname]*
    kputc(tab, out);
    ksprintf(out, "0"); /// 8) position of the mate/next read --- int [0, 2^31-1]
    kputc(tab, out);
    ksprintf(out, "%d", 0); /// 9) observed template LENgth --- int [-2^31 + 1, 2^31 - 1]
    kputc(tab, out);
    /// 10) segment SEQuence --- string \*|[A-Za-z=.]+
    if (!fwd_query) {
        ksprintf(out, "*");
    } else {
        for (int i = 0; i < query_size; ++i) {
            int c = (tpca->pca.qdir == FWD) ? fwd_query[i] : rev_query[i];
            c = DECODE_RESIDUE(c);
            kputc(c, out);
        }
    }
    /// 11) ASCII of Phred-scaled base QUALity+33 --- string [!-~]+
    kputc(tab, out);
    if (!query_qv) {
        ksprintf(out, "*"); 
    } else {
        for (int i = 0; i < query_size; ++i) {
            char c = (tpca->pca.qdir == FWD) ? query_qv[i] : query_qv[query_size - 1 - i];
            kputc(c, out);
        }
    }

    if (tpca->pca.qc >= 0) {
        kputc(tab, out);
        s_dump_repeat_pca(all_pca_a + tpca->pca.qc, query_name, subjects, out);
    }

    int num_ident = 0;
    for (int i = 0; i < as_size; ++i) if (qas[i] == sas[i]) ++num_ident;
    int map_score = 0;
    double pi = calc_ident_perc(qas, sas, as_size, NULL, &map_score);
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", map_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", as_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", map_score); ///  dp score
    kputc(tab, out);
    dump_md_string(qas, sas, as_size, out);
    /// qid
    kputc(tab, out);
    ksprintf(out, "qi:i:%d", query_id);
    /// qdir
    kputc(tab, out);
    ksprintf(out, "qd:i:%d", tpca->pca.qdir);
    /// qs
    kputc(tab, out);
    ksprintf(out, "qs:i:%d", tpca->qoff);
    /// qe
    kputc(tab, out);
    ksprintf(out, "qe:i:%d", tpca->qend);
    /// vqs
    kputc(tab, out);
    ksprintf(out, "qS:i:%d", tpca->pca.enzyme_qoff);
    /// vqe
    kputc(tab, out);
    ksprintf(out, "qE:i:%d", tpca->pca.enzyme_qend);
    // ql
    kputc(tab, out);
    ksprintf(out, "ql:i:%d", query_size);
    /// sid
    kputc(tab, out);
    ksprintf(out, "si:i:%d", tpca->pca.sid);
    /// ss
    kputc(tab, out);
    ksprintf(out, "ss:i:%d", tpca->soff);
    /// se
    kputc(tab, out);
    ksprintf(out, "se:i:%d", tpca->send);
    /// vss
    kputc('\t', out);
    ksprintf(out, "vS:i:%d", tpca->pca.enzyme_soff);
    /// vse
    kputc('\t', out);
    ksprintf(out, "vE:i:%d", tpca->pca.enzyme_send);
    ///sl
    kputc(tab, out);
    ksprintf(out, "sl:i:%d", tpca->pca.ssize);
    /// identity
    kputc(tab, out);
    ksprintf(out, "pi:f:%g", pi);
    /// global chain score
    kputc(tab, out);
    ksprintf(out, "gs:i:%d", tpca->pca.chain_score);
    /// complete map
    //kputc(tab, out);
    //const char* map_type = get_chain_type_name(chain_type);
    //ksprintf(out, "mt:Z:%s", map_type);
    kputc('\n', out);
}

static int
s_dump_one_pca(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    SeqReader* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const char* query_qv,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    TrimPcaList* tpca_list,
    TrimPca* tpca,
    const EChainType chain_type,
    const int extended_subject,
    kstring_t* out)
{
    PoreCAlign* pca = &tpca->pca;
    const char* subject_name = SeqReader_SeqName(subjects, pca->sid);
    const int subject_size = SeqReader_SeqSize(subjects, pca->sid);
    ksprintf(out, "%s", query_name);
    kputc('\t', out);
    ksprintf(out, "%d", query_size);
    kputc('\t', out);
    if (pca->qdir == FWD) {
        ksprintf(out, "%d", tpca->qoff);
        kputc('\t', out);
        ksprintf(out, "%d", tpca->qend);
        kputc('\t', out);
    } else {
        ksprintf(out, "%d", query_size - tpca->qend);
        kputc('\t', out);
        ksprintf(out, "%d", query_size - tpca->qoff);
        kputc('\t', out);
    }
    char strand = (pca->qdir == FWD) ? '+' : '-';
    kputc(strand, out);
    kputc('\t', out);

    ksprintf(out, "%s", subject_name);
    kputc('\t', out);
    ksprintf(out, "%d", subject_size);
    kputc('\t', out);
    ksprintf(out, "%d", tpca->soff);
    kputc('\t', out);
    ksprintf(out, "%d", tpca->send);
    kputc('\t', out);

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    const int as_size = tpca->as_size;
    int mat = 0;
    for (int i = 0; i < as_size; ++i) mat += (qas[i] == sas[i]);
    ksprintf(out, "%d", mat);
    kputc('\t', out);
    ksprintf(out, "%d", as_size);
    kputc('\t', out);
    ksprintf(out, "%d", tpca->pca.map_q);

    kputc('\t', out);
    if (pca->qdir == FWD) {
        ksprintf(out, "qS:i:%d", pca->enzyme_qoff);
    } else {
        int enzyme_size = reloci_list->enzyme.enzyme_size;
        int x = (pca->enzyme_qend == pca->qsize) ? 0 : (pca->qsize - pca->enzyme_qend - enzyme_size);
        ksprintf(out, "qS:i:%d", x);
    }

    kputc('\t', out);
    if (pca->qdir == FWD) {
        ksprintf(out, "qE:i:%d", pca->enzyme_qend);
    } else {
        int enzyme_size = reloci_list->enzyme.enzyme_size;
        int x = (pca->enzyme_qoff == 0) ? pca->qsize : (pca->qsize - pca->enzyme_qoff - enzyme_size);
        ksprintf(out, "qE:i:%d", x);
    }

    kputc('\t', out);
    ksprintf(out, "vS:i:%d", pca->enzyme_soff);
    kputc('\t', out);
    ksprintf(out, "vE:i:%d", pca->enzyme_send);

    double pi = calc_ident_perc(qas, sas, as_size, NULL, NULL);

    kputc('\t', out);
    ksprintf(out, "pi:f:%g", pi);

    /// global chain score
    kputc('\t', out);
    ksprintf(out, "gs:i:%d", tpca->pca.chain_score);

    kputc('\t', out);
    ksprintf(out, "Es:i:%d", extended_subject);

    //const char* map_type = get_chain_type_name(chain_type);
    //kputc('\t', out);
    //ksprintf(out, "mt:Z:%s", map_type);

    kputc('\t', out);
    print_paf_cigar(tpca->qoff, tpca->qend, query_size, qas, sas, as_size, out);
    kputc('\t', out);
    dump_md_string(qas, sas, as_size, out);

    if (pca->qc >= 0) {
        kputc('\t', out);
        s_dump_repeat_pca(all_pca_a + pca->qc, query_name, subjects, out);
    }

    kputc('\n', out);

    return 1;
}

bool s_is_complete_chain(const int* vdfa, const int vdfc, const PoreCAlign* pca_a, const int pca_c);

void
trim_overlap_subseqs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    SeqReader* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const char* query_qv,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    PoreCAlign* pca_a,
    int pca_c,
    const EChainType chain_type,
    int extended_subject,
    EOutputFmt outfmt,
    kstring_t* out)
{
    if (1)
    {
        const int* vdfa = kv_data(qvep_list->fwd_vdf_endpoint_list);
        const int vdfc = kv_size(qvep_list->fwd_vdf_endpoint_list);
        bool is_complete_chain = s_is_complete_chain(vdfa, vdfc, pca_a, pca_c);
        if (!is_complete_chain) {
            for (int i = 0; i < pca_c; ++i) {
                if (pca_a[i].chain_qend - pca_a[i].chain_qoff >= 500) continue;
                int ld = abs(pca_a[i].soff - pca_a[i].enzyme_soff) <= 20;
                int rd = abs(pca_a[i].send - pca_a[i].enzyme_send) <= 20;
                if (ld || rd) continue;
                //fprintf(stderr, "garbage pca\n");
                //dump_chain_pca(fprintf, stderr, pca_a[i], i);
                pca_a[i].qid = -1;
            }
        } else if (1) {
		for (int i = 0; i < pca_c; ++i) {
            if (pca_a[i].qend - pca_a[i].qoff > 30) continue;
			int lsd = abs(pca_a[i].soff - pca_a[i].enzyme_soff) <= 20;	
			int rsd = abs(pca_a[i].send - pca_a[i].enzyme_send) <= 20;
			int lqd = 0;//abs(pca_a[i].qoff - pca_a[i].enzyme_qoff) <= 20;	
			int rqd = 0;//abs(pca_a[i].qend - pca_a[i].enzyme_qend) <= 20;
			if (lsd == 0 && rsd == 0 && lqd == 0 && rqd == 0) {
				//fprintf(stderr, "garbage complte pca\n");
				//dump_chain_pca(fprintf, stderr, pca_a[i], i);
				pca_a[i].qid = -1;
			}
		}
	}
            int n = 0;
            for (int i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
            pca_c = n;
    }

    TrimPcaList all_tpca_list;
    for (int i = 0; i < pca_c; ++i) {
        PoreCAlign* pca = pca_a + i;
        //dump_chain_pca(fprintf, stderr, pca_a[i], i);
        get_pca_align_string(tbck_data, reloci_list, subjects, 
            query_name, query_id, fwd_query, rev_query, query_size,
            pca, &all_tpca_list);
    }

    TrimPcaList qry_tpca_list;
    set_trim_pca_list_for_fwd_query(&all_tpca_list, &qry_tpca_list, reloci_list->enzyme.enzyme_size);
    sort(qry_tpca_list.tpca_list.begin(),
        qry_tpca_list.tpca_list.end(),
        [](const TrimPca& a, const TrimPca& b)->bool { return a.qoff < b.qoff; });
    s_trim_query_overlap_subseq(&qry_tpca_list);

    TrimPcaList sbj_tpca_list;
    set_trim_pca_list_for_fwd_subject(&qry_tpca_list, &sbj_tpca_list, reloci_list->enzyme.enzyme_size);
    int n_pca = sbj_tpca_list.tpca_list.size();
    for (int i = 0; i < n_pca; ++i) {
        TrimPca* tpca = &sbj_tpca_list.tpca_list[i];
        if (!tpca->is_valid) continue;    
        hbn_assert(tpca->sdir == FWD);  
        PoreCAlign* pca = &tpca->pca;
        const u8* subject = SeqReader_Seq(subjects, pca->sid, FWD);  
        const char* qas = sbj_tpca_list.align_strings.c_str() + tpca->qas_offset;
        const char* sas = sbj_tpca_list.align_strings.c_str() + tpca->sas_offset;
        int as_size = tpca->as_size;
        validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            pca->qid, (pca->qdir ==  FWD) ? fwd_query : rev_query, tpca->qoff, tpca->qend, qas,
            pca->sid, subject, tpca->soff, tpca->send, sas, as_size, TRUE);
    }
    s_trim_subject_overlap_subseq(&sbj_tpca_list);
    n_pca = sbj_tpca_list.tpca_list.size();
    for (int i = 0; i < n_pca; ++i) {
        TrimPca* tpca = &sbj_tpca_list.tpca_list[i];
        if (!tpca->is_valid) continue;      
        PoreCAlign* pca = &tpca->pca;
        const u8* subject = SeqReader_Seq(subjects, pca->sid, FWD);  
        const char* qas = sbj_tpca_list.align_strings.c_str() + tpca->qas_offset;
        const char* sas = sbj_tpca_list.align_strings.c_str() + tpca->sas_offset;
        int as_size = tpca->as_size;
        validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
            pca->qid, (pca->qdir ==  FWD) ? fwd_query : rev_query, tpca->qoff, tpca->qend, qas,
            pca->sid, subject, tpca->soff, tpca->send, sas, as_size, TRUE);
    }

    int seq_has_dummped = 0;
    for (int i = 0; i < n_pca; ++i) {
        TrimPca* tpca = &sbj_tpca_list.tpca_list[i];
        if (!tpca->is_valid) continue;
        if (chain_type == eMaxCovChain && tpca->pca.map_q < 10) {
            //HBN_LOG("1 repeat pca");
            //dump_chain_pca(fprintf, stderr, tpca->pca, -1);
		    continue;
	    }
        if (outfmt == eOutputFmt_PAF) {
            s_dump_one_pca(tbck_data, reloci_list, subjects, query_name, query_id, 
                fwd_query, rev_query, query_qv, query_size, all_pca_a, all_pca_c, &sbj_tpca_list, tpca, 
                chain_type, extended_subject, out);
        } else if (outfmt == eOutputFmt_SAM) {
            const u8* fwd_q = NULL;
            const u8* rev_q = NULL;
            const char* qv = NULL;
            if (!seq_has_dummped) {
                fwd_q = fwd_query;
                rev_q = rev_query;
                qv = query_qv;
                seq_has_dummped = 1;
            }
            s_dump_sam_pca(tbck_data, reloci_list, subjects, query_name, query_id, 
                fwd_q, rev_q, qv, query_size, all_pca_a, all_pca_c, &sbj_tpca_list, tpca, 
                chain_type, extended_subject, out);
        }
    }
}
