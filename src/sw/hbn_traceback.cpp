#include "hbn_traceback.hpp"

#include <math.h>
#include <cstring>

#include "hbn_traceback_aux.h"

using namespace std;

static const int kMatLen = 8;
static const int kMaxOverHang = 3000;

void
validate_mem(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const DDFS_Seed* cdpsa,
    const int cdpsc)
{
    //return;
    //HBN_LOG("validating mems from [%s, %s, %d]", HBN_LOG_ARGS_GENERIC);
    for (int i = 0; i < cdpsc; ++i) {
        DDFS_Seed s = cdpsa[i];
        //fprintf(stderr, "\tvalidating %d\t%d\t%d\t%d\t%d\n", i, s.qoff, s.soff, s.length, s.sdir);
        int qi = s.qoff;
        int si = s.soff;
        for (int k = 0; k < s.length; ++k, ++qi, ++si) {
            hbn_assert(read[qi] == subject[si], "[%s, %s, %d] at (%d, %d, %d): qi = %d, si = %d, q = %d, s = %d", 
            HBN_LOG_ARGS_GENERIC, s.qoff, s.soff, s.length, qi, si, read[qi], subject[si]);
        }
    }
}

static void
compute_trace_points(const DDFS_Seed* seed_array,
    const int seed_count,
    vector<DDFS_Seed>& trace_seeds)
{
    trace_seeds.clear();
    if (seed_count == 0) return;
    int i = 0;
    while (i < seed_count) {
        int f = i;
        while (f < seed_count - 1) {
            int g = f + 1;
            DDFS_Seed sf = seed_array[f];
            DDFS_Seed sg = seed_array[g];
            hbn_assert(sf.qoff < sg.qoff, 
                "f = %d, g = %d, sf = (%d, %d, %d), sg = (%d, %d, %d)", 
                f, g, sf.qoff, sf.soff, sf.length, sg.qoff, sg.soff, sg.length);
            hbn_assert(sf.soff < sg.soff,
                "seed_count = %d, f = %d, g = %d, sf = (%d, %d, %d), sg = (%d, %d, %d)", 
                seed_count, f, g, sf.qoff, sf.soff, sf.length, sg.qoff, sg.soff, sg.length);
            int dq = sg.qoff - sf.qoff;
            int ds = sg.soff - sf.soff;
            int r1 = abs(dq - ds) <= 100;
            int r2 = fabs(1.0 - 1.0 * dq / ds) <= 0.21;
            if ((!r1) || (!r2)) break;
            ++f;            
        }

        if (f > i) {
            for (int k = i; k <= f; ++k) {
                trace_seeds.push_back(seed_array[k]);
            }
        } else if (i == 0) {
            trace_seeds.push_back(seed_array[0]);
        }

        i = f + 1;
    }

    if (trace_seeds.back().qoff < seed_array[seed_count-1].qoff) {
        trace_seeds.push_back(seed_array[seed_count-1]);
    }
}

static void
apped_match_subseq(const u8* query,
    const int qfrom,
    const int qto,
    const u8* subject,
    const int sfrom,
    const int sto,
    char** qae,
    char** sae)
{
    //HBN_LOG("qf = %d, qt = %d, sf = %d, st = %d", qfrom, qto, sfrom, sto);
    hbn_assert(qto - qfrom == sto - sfrom);
    int n = qto - qfrom;
    for (int i = 0; i < n; ++i) {
        int c1 = query[qfrom + i];
        c1 = DECODE_RESIDUE(c1);
        **qae = c1;
        ++(*qae);
        int c2 = subject[sfrom + i];
        c2 = DECODE_RESIDUE(c2);
        **sae = c2;
        ++(*sae);
        hbn_assert(c1 == c2, "qfrom = %d, qto = %d, sfrom = %d, sto = %d, qi = %d, si = %d, c1 = %c, c2 = %c", 
            qfrom, qto, sfrom, sto,
            qfrom + i, sfrom + i, c1, c2);
    }
    **qae = '\0';
    **sae = '\0';
}

static int
run_nw(const u8* query, 
    const int query_length,
    const u8* subject, 
    const int subject_length,
    HbnTracebackData* data,
    string& qaln,
    string& saln,
    char** qae,
    char** sae)
{
    if (query_length < SMALL_EDLIB_MAX_SEQ_SIZE
        &&
        subject_length < SMALL_EDLIB_MAX_SEQ_SIZE) {
        small_edlib_nw(query, query_length,
            subject, subject_length,
            data->small_edlib_data,
            qaln,
            saln);
    } else {
        edlib_nw(data->edlib, query, query_length, subject, subject_length, -1, qaln, saln);
    }
    hbn_assert(qaln.size() == saln.size());
    hbn_assert((*qae) != NULL);
    hbn_assert((*sae) != NULL);
    memcpy(*qae, qaln.c_str(), qaln.size());
    *qae += qaln.size();
    **qae = '\0';
    memcpy(*sae, saln.c_str(), saln.size());
    *sae += saln.size();
    **sae = '\0';
//const char* qas = ks_s(*qaln);
//const char* sas = ks_s(*saln);
//int as_size = ks_size(*qaln);
//dump_align_string(qas, sas, as_size, stderr);
    return 1;
}

static int
overhang_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    string& qaln,
    string& saln)
{
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
                0.65,
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
    return r;
}

static int
left_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const u8* subject,
    vector<u8>& qsbuf,
    vector<u8>& ssbuf,
    string& qabuf,
    string& sabuf,
    int* qbeg,
    int* sbeg,
    char** qas,
    char** sas)
{
    int qls = *qbeg;
    int sls = *sbeg;
    int ls = hbn_min(qls, sls);
    if (ls > kMaxOverHang) return 0;
    if (ls == 0) return 0;
//fprintf(stderr, "--qls = %d, sls = %d\n", qls, sls);

    {
        int x = hbn_min(qls, sls * 1.3);
        int y = hbn_min(qls * 1.3, sls);
	if (x > y) {
		if (x - y < 50) x = y + 50 > qls ? qls : y + 50;
	} else {
		if (y - x < 50) y = x + 50 > sls ? sls : x + 50;
	}
        qls = x;
        sls = y;
    }
    
    qls += kMatLen;
    sls += kMatLen;
//fprintf(stderr, "qls = %d, sls = %d\n", qls, sls);

    qsbuf.clear();
    ssbuf.clear();
    for (ls = 0; ls < qls; ++ls) {
        int i = *qbeg + kMatLen - ls - 1;
        hbn_assert(i >= 0);
        qsbuf.push_back(query[i]);
    }
    for (ls = 0; ls < sls; ++ls) {
        int i = *sbeg + kMatLen - ls - 1;
        hbn_assert(i >= 0);
        ssbuf.push_back(subject[i]);
    }

    int r = overhang_extend(dalign, edlib, qsbuf.data(), qls, ssbuf.data(), sls, qabuf, sabuf);
    if (!r) return r;
    *qas += kMatLen;
    *sas += kMatLen;
    *qbeg += kMatLen;
    *sbeg += kMatLen;
    int qi = 0, si = 0;
    hbn_assert(qabuf.size() == sabuf.size());
    for (size_t i = 0; i < qabuf.size(); ++i) {
        char qc = qabuf[i];
        --(*qas);
        **qas = qc;
        if (qc != GAP_CHAR) ++qi;
        char sc = sabuf[i];
        --(*sas);
        **sas = sc;
        if (sc != GAP_CHAR) ++si;
//fprintf(stderr, "%c\t%c\n", qc, sc);
    }

    *qbeg -= qi;
    *sbeg -= si;
    return 1;
}

static int
right_extend(DalignData* dalign,
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
    int qrs = query_length - (*qend);
    int srs = subject_length - (*send);
    int rs = hbn_min(qrs, srs);
    if (rs > kMaxOverHang || rs == 0) return 0;

    {
        int x = hbn_min(qrs, srs * 1.3);
        int y = hbn_min(qrs * 1.3, srs);
	if (x > y) {
		if (x - y < 50) x = y + 50 > qrs ? qrs : y + 50;
	} else {
		if (y - x < 50) y = x + 50 > srs ? srs : x + 50;
	}
        qrs = x;
        srs = y;
    }

    qrs += kMatLen;
    srs += kMatLen;
    const u8* q = query + (*qend) - kMatLen;
    const u8* s = subject + (*send) - kMatLen;

    int r = overhang_extend(dalign, edlib, q, qrs, s, srs, qabuf, sabuf);
    if (!r) return r;

    *qend -= kMatLen;
    *send -= kMatLen;
    *qae -= kMatLen;
    *sae -= kMatLen;
    int qi = 0;
    int si = 0;
    hbn_assert(qabuf.size() == sabuf.size());
    for (size_t i = 0; i < qabuf.size(); ++i) {
        char qc = qabuf[i];
        **qae = qc; ++(*qae);
        if (qc != GAP_CHAR) ++qi;
        char sc = sabuf[i];
        **sae = sc; ++(*sae);
        if (sc != GAP_CHAR) ++si;
    }
    **qae = '\0';
    **sae = '\0';
    *qend += qi;
    *send += si;
    return 1;
}

int
hbn_traceback(HbnTracebackData* data,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    const DDFS_Seed* seed_array,
    const int seed_count,
    const int min_align_size,
    const double min_ident_perc,
    const int process_over_hang)
{
    if (seed_count == 0) { return 0; }
    validate_mem(HBN_LOG_ARGS_DEFAULT, query, subject, seed_array, seed_count);
    const int kTracebackExtendBlock = 1024;
    compute_trace_points(seed_array, seed_count, data->trace_seeds);
    const DDFS_Seed* tsa = seed_array;//kv_data(data->trace_seeds);
    int tsc = seed_count;//kv_size(data->trace_seeds);
    for (int i = 0; i < tsc - 1; ++i) {
        hbn_assert(tsa[i].qoff <= tsa[i+1].qoff && tsa[i].soff <= tsa[i+1].soff,
            "i = %d, iqb = %d, isb = %d, il = %d, jqb = %d, jsb = %d, jl = %d",
            i, tsa[i].qoff, tsa[i].soff, tsa[i].length, tsa[i+1].qoff, tsa[i+1].soff, tsa[i+1].length);
    }
    //for (int i = 0; i < seed_count; ++i) fprintf(stderr, "%d\t%d\t%d\t%d\n", i, seed_array[i].qoff, seed_array[i].soff, seed_array[i].length);
    //for (int i = 0; i < tsc; ++i) fprintf(stderr, "%d\t%d\t%d\t%d\n", i, tsa[i].qoff, tsa[i].soff, tsa[i].length);
    data->init_ul(tsa[0].qoff, query_length, tsa[0].soff, subject_length);
    char** qae = &data->qae;
    char** sae = &data->sae;
    const int E = 8;
    const int E2 = E * 2;
    int qfrom = 0, qto = 0;
    int sfrom = 0, sto = 0;
    int r = 1;
    if (tsa[0].length > E2) {
        qfrom = tsa[0].qoff;
        qto = tsa[0].qoff + tsa[0].length;
        sfrom = tsa[0].soff;
        sto = tsa[0].soff + tsa[0].length;
        apped_match_subseq(query, qfrom, qto - E, subject, sfrom, sto - E, qae, sae);
        qfrom = qto - E;
        sfrom = sto - E;
    } else {
        qfrom = tsa[0].qoff;
        sfrom = tsa[0].soff;
    }

    for (int i = 0; i < tsc - 1; ++i) {
        DDFS_Seed sj = tsa[i+1];
        //HBN_LOG("%d: (%d, %d, %d)", i, sj.qoff, sj.soff, sj.length);
        if (sj.length > E2) {
            qto = sj.qoff + E;
            sto = sj.soff + E;
            hbn_assert(qto >= qfrom && sto >= sfrom, 
                "i = %d, qfrom = %d, qto = %d, sfrom = %d, sto = %d, qf1 = %d, sf1 = %d, l1 = %d, qf2 = %d, sf2 = %d, l2 = %d",
                    i, qfrom, qto, sfrom, sto, tsa[i].qoff, tsa[i].soff, tsa[i].length, tsa[i+1].qoff, tsa[i+1].soff, tsa[i+1].length);
            run_nw(query + qfrom, qto - qfrom, subject + sfrom, sto - sfrom, data, data->ext_qabuf, data->ext_sabuf, qae, sae);
            hbn_assert(r);
            qfrom = qto;
            qto = sj.qoff + sj.length - E;
            sfrom = sto;
            sto = sj.soff + sj.length - E;
            apped_match_subseq(query, qfrom, qto, subject, sfrom, sto, qae, sae);
            qfrom = qto;
            sfrom = sto;
        } else {
            qto = sj.qoff + sj.length / 2;
            sto = sj.soff + sj.length / 2;
//fprintf(stderr, "qfrom = %d, qto = %d, sfrom = %d, sto = %d\n", qfrom, qto, sfrom, sto);
            run_nw(query + qfrom, qto - qfrom, subject + sfrom, sto - sfrom, data, data->ext_qabuf, data->ext_sabuf, qae, sae);
            hbn_assert(r);
            qfrom = qto;
            sfrom = sto;
        }
    }

    DDFS_Seed se = tsa[tsc-1];
    qto = se.qoff + se.length;
    sto = se.soff + se.length;
    apped_match_subseq(query, qfrom, qto, subject, sfrom, sto, qae, sae);
    hbn_assert(strlen(data->qas) == strlen(data->sas));

    data->qoff = tsa[0].qoff;
    data->qend = qto;
    data->qsize = query_length;
    data->soff = tsa[0].soff;
    data->send = sto;
    data->ssize = subject_length;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    
    int qcnt = 0, scnt = 0;
    char* qa = data->qas;
    char* sa = data->sas;
    int m = 0;
    while (qa < data->qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    if (m < kMatLen) { return 0; }
    data->qas = qa - kMatLen;
    data->sas = sa - kMatLen;
    data->qoff += qcnt - kMatLen;
    data->soff += scnt - kMatLen;   

    qcnt = 0;
    scnt = 0;
    qa = data->qae;
    sa = data->sae;
    m = 0;
    while (qa > data->qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    hbn_assert(m == kMatLen, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == data->qas) { return 0; }

    data->qae = qa + kMatLen;
    data->sae = sa + kMatLen;
    data->qend -= qcnt - kMatLen;
    data->send -= scnt - kMatLen;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);

    if (1) { //(data->qoff <= kMaxEdlibOverHang || data->soff <= kMaxEdlibOverHang) {
    edlib_extend(data->edlib,
        query + data->qoff - 1,
        data->qoff,
        subject + data->soff - 1,
        data->soff,
        kTracebackExtendBlock,
        FALSE,
        data->qfrag,
        data->sfrag,
        &qcnt,
        &scnt,
        data->ext_qabuf,
        data->ext_sabuf);

    hbn_assert(qcnt <= data->qoff);
    hbn_assert(scnt <= data->soff);
    hbn_assert(data->ext_qabuf.size() == data->ext_sabuf.size());
    for (size_t i = 0; i < data->ext_qabuf.size(); ++i) {
        int res = data->ext_qabuf[i];
        --data->qas;
        *data->qas = res;
        res = data->ext_sabuf[i];
        --data->sas;
        *data->sas = res;
    }
    data->qoff -= qcnt;
    data->soff -= scnt;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    }

    if (1) { //if (query_length - data->qend <= kMaxEdlibOverHang || subject_length - data->send <= kMaxEdlibOverHang) {
    edlib_extend(data->edlib,
        query + data->qend,
        query_length - data->qend,
        subject + data->send,
        subject_length - data->send,
        kTracebackExtendBlock,
        TRUE,
        data->qfrag,
        data->sfrag,
        &qcnt,
        &scnt,
        data->ext_qabuf,
        data->ext_sabuf);
    hbn_assert(data->qend + qcnt <= query_length);
    hbn_assert(data->send + scnt <= subject_length);
    hbn_assert(data->ext_qabuf.size() == data->ext_sabuf.size());
    memcpy(data->qae, data->ext_qabuf.c_str(), data->ext_qabuf.size());
    data->qae += data->ext_qabuf.size();
    *data->qae = '\0';
    memcpy(data->sae, data->ext_sabuf.c_str(), data->ext_sabuf.size());
    data->sae += data->ext_sabuf.size();
    *data->sae = '\0';
    data->qend += qcnt;
    data->send += scnt;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    }

    if (data->qas == data->qae) { return 0; }

    if (process_over_hang) {
        //HBN_LOG("before:");
        //HbnTracebackDataDump(fprintf, stderr, data);
        //left_extend(data->ksw, query, subject, &data->ksw->qfrag, &data->ksw->tfrag,
        //    &data->qoff, &data->soff, &data->qas, &data->sas);
        //right_extend(data->ksw, query, &data->qend, query_length,
        //    subject, &data->send, subject_length, &data->qae, &data->sae);

        left_extend(data->dalign, data->edlib, query, subject, data->ksw->qfrag, data->ksw->tfrag,
            data->ext_qabuf, data->ext_sabuf, &data->qoff, &data->soff, &data->qas, &data->sas);
        right_extend(data->dalign, data->edlib, query, &data->qend, query_length,
            subject, &data->send, subject_length, data->ext_qabuf, data->ext_sabuf,
            &data->qae, &data->sae);
        //HBN_LOG("after:");
        //HbnTracebackDataDump(fprintf, stderr, data);
    }
hbn_assert(data->qae - data->qas == data->sae - data->sas);

    data->ident_perc = calc_ident_perc(data->qas, data->sas, strlen(data->qas), &data->dist, &data->score);
    //hbn_assert(strlen(data->qas) == strlen(data->sas), "ql = %zu, sl = %zu", strlen(data->qas), strlen(data->sas));
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    
    r = (data->qae - data->qas >= min_align_size) && (data->ident_perc >= min_ident_perc);
    return r;
}

BOOL
truncate_align_bad_ends(const char* qaln,
    const char* saln,
    const int aln_size,
    int* qoff,
    int* qend,
    int* soff,
    int* send,
    const char** qas_,
    const char** qae_,
    const char** sas_,
    const char** sae_)
{
    const char* const qas = qaln;
    const char* const qae = qaln + aln_size;
    const char* const sas = saln;
    const char* const sae = saln + aln_size;
    int qcnt = 0, scnt = 0;
    const char* qa = qas;
    const char* sa = sas;
    int m = 0;

    while (qa < qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    if (m < kMatLen || qa == qae) return 0;
    *qas_ = qa - kMatLen;
    *sas_ = sa - kMatLen;
    *qoff += qcnt - kMatLen;
    *soff += scnt - kMatLen;  

    qcnt = 0;
    scnt = 0;
    qa = qae;
    sa = sae;
    m = 0;
    while (qa > qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    hbn_assert(m == kMatLen, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == qas) return 0;

    *qae_ = qa + kMatLen;
    *sae_ = sa + kMatLen;
    *qend -= qcnt - kMatLen;
    *send -= scnt - kMatLen; 
    return TRUE;
}

static BOOL
truncate_align_bad_ends_1(char* qaln,
    char* saln,
    const int aln_size,
    int* qoff,
    int* qend,
    int* soff,
    int* send,
    char** qas_,
    char** qae_,
    char** sas_,
    char** sae_)
{
    char* const qas = qaln;
   char* const qae = qaln + aln_size;
    char* const sas = saln;
    char* const sae = saln + aln_size;
    int qcnt = 0, scnt = 0;
    char* qa = qas;
    char* sa = sas;
    int m = 0;
    const int M = 8;

    while (qa < qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == M) break;
    }
    if (m < M || qa == qae) return 0;
    *qas_ = qa - M;
    *sas_ = sa - M;
    *qoff += qcnt - M;
    *soff += scnt - M;  

    qcnt = 0;
    scnt = 0;
    qa = qae;
    sa = sae;
    m = 0;
    while (qa > qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == M) break;
    }
    hbn_assert(m == M, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == qas) return 0;

    *qae_ = qa + M;
    *sae_ = sa + M;
    *qend -= qcnt - M;
    *send -= scnt - M; 
    return TRUE;
}

int 
porec_compute_traceback(HbnTracebackData* data,
    int qb,
    int qe,
    int sb,
    int se,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    const int min_align_size,
    const double min_ident_perc)
{
	int r = 1;
    data->init(qb, query_length, sb, subject_length);
#if 0
    int ql = qe - qb;
    int sl = se - sb;
    int max_l = hbn_max(ql, sl);
    int min_l = hbn_min(ql, sl);
    int tol = max_l * 0.25;
    if (tol < 100) tol = 100;
    if (tol > min_l) tol = min_l;
    if (!edlib_nw(data->edlib, query + qb, qe - qb, subject + sb, se - sb, tol, &data->dalign->aseq, &data->dalign->bseq)) return FALSE;
    const char* qas = ks_s(data->dalign->aseq);
    const char* sas = ks_s(data->dalign->bseq);
    int as_size = ks_size(data->dalign->aseq);
    memcpy(data->qae, qas, as_size);
    data->qae += as_size;
    *data->qae = '\0';
    memcpy(data->sae, sas, as_size);
    data->sae += as_size;
    *data->sae = '\0';
#else
    run_nw(query + qb, qe - qb, subject + sb, se - sb, data, 
        data->edlib->sqaln, data->edlib->staln, &data->qae, &data->sae);
#endif
    data->qoff = qb;
    data->qend = qe;
    data->qsize = query_length;
    data->soff = sb;
    data->send = se;
    data->ssize = subject_length;
    int as_size = data->qae - data->qas;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size,
        1);
    
    int qcnt = 0, scnt = 0;
    char* qa = data->qas;
    char* sa = data->sas;
    int m = 0;
    while (qa < data->qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    if (m < kMatLen) return 0;
    data->qas = qa - kMatLen;
    data->sas = sa - kMatLen;
    data->qoff += qcnt - kMatLen;
    data->soff += scnt - kMatLen;   

    qcnt = 0;
    scnt = 0;
    qa = data->qae;
    sa = data->sae;
    m = 0;
    while (qa > data->qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == kMatLen) break;
    }
    hbn_assert(m == kMatLen, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == data->qas) return 0;

    data->qae = qa + kMatLen;
    data->sae = sa + kMatLen;
    data->qend -= qcnt - kMatLen;
    data->send -= scnt - kMatLen;

    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend,
        data->qas,
        0, subject, data->soff, data->send, data->sas,
        strlen(data->qas),
        1);
    if (data->qas == data->qae) return 0;

    hbn_assert(data->qae - data->qas == data->sae - data->sas);

    as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);

    r = truncate_align_bad_ends_1(data->qas, data->sas, as_size,
            &data->qoff, &data->qend, &data->soff, &data->send,
            &data->qas, &data->qae, &data->sas, &data->sae);
    if (!r) return r;
    as_size = data->qae - data->qas;
    data->ident_perc = calc_ident_perc(data->qas, data->sas, as_size, &data->dist, &data->score);
    
    r = (data->qae - data->qas >= min_align_size) && (data->ident_perc >= min_ident_perc);
    return r;    
}
