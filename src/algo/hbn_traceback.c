#include "hbn_traceback.h"

#include <math.h>

#include "hbn_traceback_aux.h"

static const int kMatLen = 8;
static const int kMaxOverHang = 3000;

HbnTracebackData*
HbnTracebackDataNew()
{
    HbnTracebackData* data = (HbnTracebackData*)calloc(1, sizeof(HbnTracebackData));
    ks_init(data->qabuf);
    ks_init(data->sabuf);
    ks_init(data->ext_qabuf);
    ks_init(data->ext_sabuf);
    kv_init(data->qfrag);
    kv_init(data->sfrag);
    data->edlib = EdlibAlignDataNew();
    data->ksw = Ksw2DataNew();
    ksw2_extd2_set_params(data->ksw);
    data->ksw->band_width = kMaxOverHang * 0.2;
    data->ksw->zdrop = 40;
    data->dalign = DalignDataNew(0.35);
    data->small_edlib = small_edlib_align_struct_new();
    return data;
}

HbnTracebackData*
HbnTracebackDataFree(HbnTracebackData* data)
{
    ks_destroy(data->qabuf);
    ks_destroy(data->sabuf);
    ks_destroy(data->ext_qabuf);
    ks_destroy(data->ext_sabuf);
    kv_destroy(data->qfrag);
    kv_destroy(data->sfrag);
    EdlibAlignDataFree(data->edlib);
    DalignDataFree(data->dalign);
    Ksw2DataFree(data->ksw);
    small_edlib_align_struct_free(data->small_edlib);
    free(data);
    return NULL;
}

static void
HbnTracebackDataInit(HbnTracebackData* data,
    int qoff,
    int qsize,
    int soff,
    int ssize)
{
    int wrk_l = 4 * hbn_max(qsize, ssize);
    ks_reserve(&data->qabuf, wrk_l);
    ks_reserve(&data->sabuf, wrk_l);
    int wrk_ll = 2 * hbn_max(qsize, ssize);
    data->qas = data->qae = ks_s(data->qabuf) + wrk_ll;
    data->sas = data->sae = ks_s(data->sabuf) + wrk_ll;
    *data->qae = '\0';
    *data->sae = '\0';
    data->qoff = data->qend = 0;
    data->soff = data->send = 0;
    data->dist = 0;
    data->ident_perc = 0;
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
    kstring_t* qaln,
    kstring_t* saln,
    char** qae,
    char** sae)
{
    if (query_length < SMALL_EDLIB_MAX_SEQ_SIZE
        &&
        subject_length < SMALL_EDLIB_MAX_SEQ_SIZE) {
        small_edlib_nw(query, query_length,
            subject, subject_length,
            data->small_edlib,
            qaln,
            saln);
    } else {
        edlib_nw(data->edlib, query, query_length, subject, subject_length, -1, qaln, saln);
    }
    //edlib_nw(data->edlib, query, query_length, subject, subject_length, -1, qaln, saln);
    hbn_assert(ks_size(*qaln) == ks_size(*saln));
    hbn_assert((*qae) != NULL);
    hbn_assert((*sae) != NULL);
    memcpy(*qae, ks_s(*qaln), ks_size(*qaln));
    *qae += ks_size(*qaln);
    **qae = '\0';
    memcpy(*sae, ks_s(*saln), ks_size(*saln));
    *sae += ks_size(*saln);
    **sae = '\0';
    return 1;
}

static int
overhang_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    kstring_t* qaln,
    kstring_t* saln)
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
    return r;
}

static int
left_extend(DalignData* dalign,
    EdlibAlignData* edlib,
    const u8* query,
    const u8* subject,
    vec_u8* qsbuf,
    vec_u8* ssbuf,
    kstring_t* qabuf,
    kstring_t* sabuf,
    int* qbeg,
    int* sbeg,
    char** qas,
    char** sas)
{
    int qls = *qbeg;
    int sls = *sbeg;
    int ls = hbn_min(qls, sls);
    if (ls == 0) return 0;

    qls += kMatLen;
    sls += kMatLen;
    kv_clear(*qsbuf);
    kv_clear(*ssbuf);

    int q_l = hbn_min(qls, sls + 100);
    int s_l = hbn_min(sls, qls + 100);
    for (int i = 0; i < q_l; ++i) kv_push(u8, *qsbuf, query[qls-i-1]);
    for (int i = 0; i < s_l; ++i) kv_push(u8, *ssbuf, subject[sls-i-1]);
    int r = overhang_extend(dalign, edlib, kv_data(*qsbuf), q_l, kv_data(*ssbuf), s_l, qabuf, sabuf);
    if (!r) return r;

    *qas += kMatLen;
    *sas += kMatLen;
    *qbeg += kMatLen;
    *sbeg += kMatLen;
    int qi = 0, si = 0;
    hbn_assert(ks_size(*qabuf) == ks_size(*sabuf));
    for (size_t i = 0; i < ks_size(*qabuf); ++i) {
        char qc = ks_A(*qabuf, i);
        --(*qas);
        **qas = qc;
        if (qc != GAP_CHAR) ++qi;
        char sc = ks_A(*sabuf, i);
        --(*sas);
        **sas = sc;
        if (sc != GAP_CHAR) ++si;
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
    kstring_t* qabuf,
    kstring_t* sabuf,
    char** qae,
    char** sae)
{
    int qrs = query_length - (*qend);
    int srs = subject_length - (*send);
    int rs = hbn_min(qrs, srs);
    if (rs == 0) return 0;

    qrs += kMatLen;
    srs += kMatLen;

    int q_l = hbn_min(qrs, srs + 100);
    int s_l = hbn_min(srs, qrs + 100);
    const u8* q = query + (*qend) - kMatLen;
    const u8* s = subject + (*send) - kMatLen;
    int r = overhang_extend(dalign, edlib, q, q_l, s, s_l, qabuf, sabuf);
    if (!r) return r;

    *qend -= kMatLen;
    *send -= kMatLen;
    *qae -= kMatLen;
    *sae -= kMatLen;
    int qi = 0;
    int si = 0;
    hbn_assert(ks_size(*qabuf) == ks_size(*sabuf));
    for (size_t i = 0; i < ks_size(*qabuf); ++i) {
        char qc = ks_A(*qabuf, i);
        **qae = qc; ++(*qae);
        if (qc != GAP_CHAR) ++qi;
        char sc = ks_A(*sabuf, i);
        **sae = sc; ++(*sae);
        if (sc != GAP_CHAR) ++si;
    }
    **qae = '\0';
    **sae = '\0';
    *qend += qi;
    *send += si;
    return 1;
}

static BOOL
truncate_align_bad_ends_1(const char* qaln,
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
    HbnTracebackDataInit(data, qb, query_length, sb, subject_length);
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
    run_nw(query + qb, qe - qb, subject + sb, se - sb, data, &data->dalign->aseq, &data->dalign->bseq, &data->qae, &data->sae);
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
    size_t aln_l = data->qae - data->qas;
    if (aln_l >= (data->qabuf.m) || aln_l >= (data->sabuf.m)) {
	    HBN_LOG("aln_l = %zu, qbuf_l = %zu, sbuf_l = %zu", aln_l, data->qabuf.m, data->sabuf.m);
	    dump_align_string(data->qas, data->sas, aln_l, stderr);
    }

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

int
porec_align_overhang(HbnTracebackData* data,
    const u8* query,
    const int query_length,
    const u8* subject,
    const int subject_length,
    const int min_align_size,
    const double min_ident_perc)
{
    left_extend(data->dalign, data->edlib, query, subject, &data->ksw->qfrag, &data->ksw->tfrag,
        &data->ext_qabuf, &data->ext_sabuf, &data->qoff, &data->soff, &data->qas, &data->sas);

    right_extend(data->dalign, data->edlib, query, &data->qend, query_length,
        subject, &data->send, subject_length, &data->ext_qabuf, &data->ext_sabuf,
        &data->qae, &data->sae);

    hbn_assert(data->qae - data->qas == data->sae - data->sas);
    size_t aln_l = data->qae - data->qas;
    if (aln_l >= (data->qabuf.m) || aln_l >= (data->sabuf.m)) {
	    HBN_LOG("aln_l = %zu, qbuf_l = %zu, sbuf_l = %zu", aln_l, data->qabuf.m, data->sabuf.m);
	    dump_align_string(data->qas, data->sas, aln_l, stderr);
    }

    //data->ident_perc = calc_ident_perc(data->qas, data->sas, strlen(data->qas), &data->dist, &data->score);
    //hbn_assert(strlen(data->qas) == strlen(data->sas), "ql = %zu, sl = %zu", strlen(data->qas), strlen(data->sas));
    int as_size = data->qae - data->qas;
    validate_aligned_string(HBN_LOG_ARGS_DEFAULT,
        0, query, data->qoff, data->qend, data->qas,
        0, subject, data->soff, data->send, data->sas,
        as_size, 1);
//fprintf(stderr, "before truncate: [%d, %d] x [%d, %d]\n", data->qoff, data->qend, data->soff, data->send);
	    //dump_align_string(data->qas, data->sas, aln_l, stderr);

#if 1
    int r = truncate_align_bad_ends_1(data->qas, data->sas, as_size,
                &data->qoff, &data->qend, &data->soff, &data->send,
                &data->qas, &data->qae, &data->sas, &data->sae);
    if (!r) return r;
#endif
    as_size = data->qae - data->qas;
    data->ident_perc = calc_ident_perc(data->qas, data->sas, as_size, &data->dist, &data->score);
    
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
