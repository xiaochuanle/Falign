#include "smooth_pca_list.hpp"

#include <algorithm>
#include <cmath>

using namespace std;

static bool 
s_gf_p2_left_extend(int qb, int qe, int eqb, int eqe, int left_qf,
    PoreCAlign* pca,
    SeqReader* ref,
    const u8* fwd_read,
    const u8* rev_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    int qf = left_qf;
    int qt = qe;
    int sf = pca->enzyme_soff;
    int st = pca->send;
    //HBN_LOG("le qf = %d, qt = %d, sf = %d, st = %d\n", qf, qt, sf, st);
    int q_d = qt - qf;
    int s_d = st - sf;
    if (s_d <= 0) return false;
    double ddf = fabs(1.0 - 1.0 * q_d / s_d);
    if (ddf > 0.2) return false;

    const u8* read = nullptr;
    int rl = qt - qf;
    if (pca->qdir == FWD) {
        read = fwd_read + qf;
    } else {
        read = rev_read + (pca->qsize - qt);
    }
    const u8* sbj = SeqReader_Seq(ref, pca->sid, FWD) + sf;
    int sl = st - sf;
    int tol = hbn_max(rl, sl) * 0.5;
    int r = edlib_nw(tbck_data->edlib, read, rl, sbj, sl, tol, &tbck_data->qabuf, &tbck_data->sabuf);
    if (!r) return false;

    const char* qas = ks_s(tbck_data->qabuf);
    const char* sas = ks_s(tbck_data->sabuf);
    int as_size = ks_size(tbck_data->qabuf);
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    if (ident < pca->pi - 5.0) return false;

    //HBN_LOG("Left fix"); dump_chain_pca(fprintf, stderr, *pca, -1);
    int xqb, xqe, xeqb, xeqe, xsb, xse, xesb, xese;
    if (pca->qdir == FWD) {
        xqb = qf;
        xqe = qt;
        xeqb = left_qf;
        xeqe = eqe;
    } else {
        xqb = pca->qsize - qt;
        xqe = pca->qsize - qf;
        xeqb = pca->qsize - eqe;
        xeqe = pca->qsize - left_qf;
        xqb -= enzyme_size;
        xeqb -= enzyme_size;
        xqe -= enzyme_size;
        xeqe -= enzyme_size;
    }
    xsb = sf;
    xse = st;
    xesb = pca->enzyme_soff;
    xese = pca->enzyme_send;
    //fprintf(stderr, "into [%d, %d] x [%d, %d], [%d, %d] x [%d, %d], %g\n", xqb, xqe, xsb, xse, xeqb, xeqe, xesb, xese, ident);

    pca->qoff = xqb;
    pca->qend = xqe;
    pca->enzyme_qoff = xeqb;
    pca->enzyme_qend = xeqe;
    pca->soff = xsb;
    pca->send = xse;
    pca->enzyme_soff = xesb;
    pca->enzyme_send = xese;
    set_pca_chain_offset(pca, enzyme_size);
    //dump_chain_pca(fprintf, stderr, *pca, -1);

    return true;
}

static bool 
s_gf_p2_right_extend(int qb, int qe, int eqb, int eqe, int right_qt,
    PoreCAlign* pca,
    SeqReader* ref,
    const u8* fwd_read,
    const u8* rev_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    int qf = qb;
    int qt = right_qt;
    int sf = pca->soff;
    int st = pca->enzyme_send;
    //HBN_LOG("re qf = %d, qt = %d, sf = %d, st = %d\n", qf, qt, sf, st);
    int q_d = qt - qf;
    int s_d = st - sf;
    if (s_d <= 0) return false;
    double ddf = fabs(1.0 - 1.0 * q_d / s_d);
//HBN_LOG("q_d = %d, s_d = %d, ddf = %g", q_d, s_d, ddf);
    if (ddf > 0.2) return false;

    const u8* read = nullptr;
    int rl = qt - qf;
    if (pca->qdir == FWD) {
        read = fwd_read + qf;
    } else {
        read = rev_read + (pca->qsize - qt);
    }
    const u8* sbj = SeqReader_Seq(ref, pca->sid, FWD) + sf;
    int sl = st - sf;
    int tol = hbn_max(rl, sl) * 0.5;
    int r = edlib_nw(tbck_data->edlib, read, rl, sbj, sl, tol, &tbck_data->qabuf, &tbck_data->sabuf);
    if (!r) return false;

    const char* qas = ks_s(tbck_data->qabuf);
    const char* sas = ks_s(tbck_data->sabuf);
    int as_size = ks_size(tbck_data->qabuf);
    double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
//fprintf(stderr, "pi: %g ---> %f\n", pca->pi, ident);
    if (ident < pca->pi - 5.0 && ident < 75.0) return false;

    //HBN_LOG("right fix"); dump_chain_pca(fprintf, stderr, *pca, -1);
    int xqb, xqe, xeqb, xeqe, xsb, xse, xesb, xese;
    if (pca->qdir == FWD) {
        xqb = qf;
        xqe = qt;
        xeqb = eqb;
        xeqe = right_qt;
    } else {
        xqb = pca->qsize - qt;
        xqe = pca->qsize - qf;
        xeqb = pca->qsize - right_qt;
        xeqe = pca->qsize - eqb;
        xqb -= enzyme_size;
        xeqb -= enzyme_size;
        xqe -= enzyme_size;
        xeqe -= enzyme_size;
    }
    xsb = sf;
    xse = st;
    xesb = pca->enzyme_soff;
    xese = pca->enzyme_send;
    //fprintf(stderr, "into [%d, %d] x [%d, %d], [%d, %d] x [%d, %d], %g\n", xqb, xqe, xsb, xse, xeqb, xeqe, xesb, xese, ident);

    pca->qoff = xqb;
    pca->qend = xqe;
    pca->enzyme_qoff = xeqb;
    pca->enzyme_qend = xeqe;
    pca->soff = xsb;
    pca->send = xse;
    pca->enzyme_soff = xesb;
    pca->enzyme_send = xese;
    set_pca_chain_offset(pca, enzyme_size);
    //dump_chain_pca(fprintf, stderr, *pca, -1);

    return true;
}

static bool 
s_gf_p2_left_most_extend(PoreCAlign* pca,
    SeqReader* ref,
    const u8* fwd_read,
    const u8* rev_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    vector<u8> qfrag, sfrag;
    const u8* read = (pca->qdir == FWD) ? fwd_read : rev_read;
    const u8* subject = SeqReader_Seq(ref, pca->sid, FWD);
    kstring_t* qabuf = &tbck_data->qabuf;
    kstring_t* sabuf = &tbck_data->sabuf;

    if (pca->qdir == FWD) {
        int qblk = min(pca->qoff, pca->soff + 50);
        int sblk = min(pca->qoff + 50, pca->soff);
        qfrag.clear(); for (int i = 1; i <= qblk; ++i) qfrag.push_back(read[pca->qoff-i]);
        sfrag.clear(); for (int i = 1; i <= sblk; ++i) sfrag.push_back(subject[pca->soff-i]);
        //fprintf(stderr, "lme 1 qf = %d, sf = %d, qblk = %d, sblk = %d\n", pca->qoff, pca->soff, qblk, sblk);
        if (qblk > 5000 || sblk > 5000) fprintf(stderr, "lme 1 qf = %d, sf = %d, qblk = %d, sblk = %d\n", pca->qoff, pca->soff, qblk, sblk);
        int qfae = 0, sfae = 0;
        edlib_shw(tbck_data->edlib, qfrag.data(), qblk, sfrag.data(), sblk, &qfae, &sfae, qabuf, sabuf);
        //HBN_LOG("LM fix"); dump_chain_pca(fprintf, stderr, *pca, -1);
        pca->qoff -= qfae;
        pca->chain_qoff = pca->qoff;
        pca->soff -= sfae;
        //fprintf(stderr, "into\t"); dump_chain_pca(fprintf, stderr, *pca, -1);
	return true;
    } else {
        int qt = pca->qend;
        int st = pca->send;
        int dq = pca->qsize - qt;
        int ds = pca->ssize - st;
        int qblk = min(ds + 50, dq);
        int sblk = min(dq + 50, ds);
        //HBN_LOG("lme2 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, pca->qsize, pca->ssize);
        if (qblk > 5000 || sblk > 5000)HBN_LOG("lme2 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, pca->qsize, pca->ssize);
        //HBN_LOG("RM fix"); dump_chain_pca(fprintf, stderr, *pca, -1);
        int qfae = 0, sfae = 0;
        edlib_shw(tbck_data->edlib, read + qt, qblk, subject + st, sblk, &qfae, &sfae, qabuf, sabuf);
        qt += qfae;
        st += sfae;
        hbn_assert(qt <= pca->qsize);
        hbn_assert(st <= pca->ssize);
        pca->qend = qt;
        pca->chain_qoff = pca->qsize - qt;
        pca->send = st;
        //fprintf(stderr, "into\t"); dump_chain_pca(fprintf, stderr, *pca, -1);
	return true;
    }
}

static bool 
s_gf_p2_right_most_extend(PoreCAlign* pca,
    SeqReader* ref,
    const u8* fwd_read,
    const u8* rev_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    vector<u8> qfrag, sfrag;
    const u8* read = (pca->qdir == FWD) ? fwd_read : rev_read;
    const u8* subject = SeqReader_Seq(ref, pca->sid, FWD);
    kstring_t* qabuf = &tbck_data->qabuf;
    kstring_t* sabuf = &tbck_data->sabuf;

    if (pca->qdir == FWD) {
        int qt = pca->qend;
        int st = pca->send;
        int dq = pca->qsize - qt;
        int ds = pca->ssize - st;
        int qblk = min(ds + 50, dq);
        int sblk = min(dq + 50, ds);
	if (qblk == 0 || sblk == 0) return false;
        //HBN_LOG("rme1 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, pca->qsize, pca->ssize);
        if (qblk > 5000 || sblk > 5000) HBN_LOG("rme1 qt = %d, qblk = %d, st = %d, sblk = %d, qsize = %d, ssize = %d", qt, qblk, st, sblk, pca->qsize, pca->ssize);
        //HBN_LOG("RM fix"); dump_chain_pca(fprintf, stderr, *pca, -1);
        int qfae = 0, sfae = 0;
        edlib_shw(tbck_data->edlib, read + qt, qblk, subject + st, sblk, &qfae, &sfae, qabuf, sabuf);
        qt += qfae;
        st += sfae;
        hbn_assert(qt <= pca->qsize);
        hbn_assert(st <= pca->ssize);
        pca->qend = qt;
        pca->chain_qend = qt;
        pca->send = st;
        //fprintf(stderr, "into\t"); dump_chain_pca(fprintf, stderr, *pca, -1);
	return true;
    } else {
        int qblk = min(pca->qoff, pca->soff + 50);
        int sblk = min(pca->qoff + 50, pca->soff);
	if (qblk == 0 || sblk == 0) return false;
        qfrag.clear(); for (int i = 1; i <= qblk; ++i) qfrag.push_back(read[pca->qoff-i]);
        sfrag.clear(); for (int i = 1; i <= sblk; ++i) sfrag.push_back(subject[pca->soff-i]);
        //fprintf(stderr, "rme2 qf = %d, sf = %d, qblk = %d, sblk = %d\n", pca->qoff, pca->soff, qblk, sblk);
        if (qblk > 5000 || sblk > 5000)fprintf(stderr, "rme2 qf = %d, sf = %d, qblk = %d, sblk = %d\n", pca->qoff, pca->soff, qblk, sblk);
        int qfae = 0, sfae = 0;
        edlib_shw(tbck_data->edlib, qfrag.data(), qblk, sfrag.data(), sblk, &qfae, &sfae, qabuf, sabuf);
        //HBN_LOG("RM fix"); dump_chain_pca(fprintf, stderr, *pca, -1);
        pca->qoff -= qfae;
        pca->chain_qend = pca->qsize - pca->qoff;
        pca->soff -= sfae;     
        //fprintf(stderr, "into\t"); dump_chain_pca(fprintf, stderr, *pca, -1);
	return true;
    }
}

void
s_gf_p2_set_fwd_read_offsets(PoreCAlign* pca, int& qb, int& qe, int& eqb, int& eqe, int enzyme_size)
{
    if (pca->qdir == FWD) {
        qb = pca->qoff;
        qe = pca->qend;
        eqb = pca->enzyme_qoff;
        eqe = pca->enzyme_qend;
    } else {
        qb = pca->qsize - pca->qend;
        qe = pca->qsize - pca->qoff;
        eqb = pca->qsize - pca->enzyme_qend;
        eqe = pca->qsize - pca->enzyme_qoff;
        if (qb >= enzyme_size) qb -= enzyme_size;
        qe -= enzyme_size;
        if (eqb >= enzyme_size) eqb -= enzyme_size;
        eqe -= enzyme_size;
    }
}

void
smooth_pca_list_pass2(std::vector<PoreCAlign>& chain,
    EChainType& chain_type,
    SeqReader* ref,
    const u8* fwd_read,
    const u8* rev_read,
    const int* vdfa,
    const int vdfc,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    if (chain_type == ePerfectChain) return;

    PoreCAlign* pcaa = chain.data();
    int pcac = chain.size();

    {
        int qb, qe, eqb, eqe;
        s_gf_p2_set_fwd_read_offsets(pcaa, qb, qe, eqb, eqe, enzyme_size);
        if (eqb == 0 && qb > 0) s_gf_p2_left_most_extend(pcaa, ref, fwd_read, rev_read, enzyme_size, tbck_data);
    }

    for (int i = 0; i < pcac - 1; ++i) {
        PoreCAlign* pi = pcaa + i;
        PoreCAlign* pj = pcaa + i + 1;
        if (pi->chain_qend >= pj->chain_qoff) continue;
        int iqb, iqe, ieqb, ieqe;
        s_gf_p2_set_fwd_read_offsets(pi, iqb, iqe, ieqb, ieqe, enzyme_size);
        int jqb, jqe, jeqb, jeqe;
        s_gf_p2_set_fwd_read_offsets(pj, jqb, jqe, jeqb, jeqe, enzyme_size);
        //HBN_LOG("filling gaps for");
        //fprintf(stderr, "[%d, %d] x [%d, %d], [%d, %d] x [%d, %d]\n", iqb, iqe, pi->soff, pi->send, ieqb, ieqe, pi->enzyme_soff, pi->enzyme_send);
        //fprintf(stderr, "[%d, %d] x [%d, %d], [%d, %d] x [%d, %d]\n", jqb, jqe, pj->soff, pj->send, jeqb, jeqe, pj->enzyme_soff, pj->enzyme_send);

        if (iqe == ieqe) {
            if (s_gf_p2_left_extend(jqb, jqe, jeqb, jeqe, iqe, pj, ref, fwd_read, rev_read, enzyme_size, tbck_data)) {
                continue;
            }
        }
        if (jqb == jeqb) {
            if (s_gf_p2_right_extend(iqb, iqe, ieqb, ieqe, jqb, pi, ref, fwd_read, rev_read, enzyme_size, tbck_data)) {
                continue;
            }
        }
    }

    {
        int qb, qe, eqb, eqe;
        s_gf_p2_set_fwd_read_offsets(pcaa + pcac - 1, qb, qe, eqb, eqe, enzyme_size);
        if (eqe == pcaa[0].qsize && qe != pcaa[0].qsize) s_gf_p2_right_most_extend(pcaa + pcac - 1, ref, fwd_read, rev_read, enzyme_size, tbck_data);
    }

    chain_type = pca_chain_type(vdfa, vdfc, chain.data(), chain.size());
    //HBN_LOG("chain_type = %s", get_chain_type_name(chain_type));
}

/////////////////

void
s_merge_adjacent_pca(std::vector<PoreCAlign>& pca_list,
    SeqReader* ref,
    const u8* fwd_read,
    const u8* rev_read,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    PoreCAlign* pca_a = pca_list.data();
    int pca_c = pca_list.size();
    int i = 0;
    while (i < pca_c) {
        PoreCAlign* pi = pca_a + i;
        if (pi->qid == -1) { ++i; continue; }
        int j = i + 1;
        while (j < pca_c) {
            PoreCAlign* pj = pca_a + j;
            if (pj->qid == -1) { ++j; continue; }
            if (pi->sid != pj->sid) break;
            if (pi->qdir != pj->qdir) break;
            int iqb = pi->qoff;
            int iqe = pi->qend;
            int isb = pi->soff;
            int ise = pi->send;
            int jqb = pj->qoff;
            int jqe = pj->qend;
            int jsb = pj->soff;
            int jse = pj->send;
            int q_d = jqe - iqb;
            int s_d = jse - isb;
            if (s_d <= 0) break;
            double ddf = fabs(1.0 - 1.0 * q_d / s_d);
            if (ddf > .2) break;
            const u8* sbj = SeqReader_Seq(ref, pi->sid, FWD);
            int tol = hbn_max(q_d, s_d) * 0.5;
            const u8* read = (pi->qdir == FWD) ? fwd_read : rev_read;
            int r = edlib_nw(tbck_data->edlib, read + iqb, q_d, sbj + isb, s_d, tol, &tbck_data->qabuf, &tbck_data->sabuf);
            if (!r) break;
            const char* qas = ks_s(tbck_data->qabuf);
            const char* sas = ks_s(tbck_data->sabuf);
            int as_size = ks_size(tbck_data->qabuf);
            double ident = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
	    double avg_pi = (pi->pi + pj->pi) / 2;
	    if (ident < avg_pi - 4.0 && ident < 80.0) { ++j; continue; }
            PoreCAlign pca = *pi;
            pca.qend = jqe;
            pca.send = jse;
            pca.enzyme_qend = pj->enzyme_qend;
            pca.enzyme_send = pj->enzyme_send;
            pca.pi = ident;
            set_pca_chain_offset(&pca, enzyme_size);
            //HBN_LOG("merge, pi = %g, avg_pi = %g", ident, avg_pi); dump_chain_pca(fprintf, stderr, *pi, i); dump_chain_pca(fprintf, stderr, *pj, j);
            //fprintf(stderr, "into"); dump_chain_pca(fprintf, stderr, pca, -1);
            *pi = pca;
            pj->qid = -1;
            ++j;
        }
        i = j;
    }
    int n = 0;
    for (i = 0; i < pca_c; ++i) if (pca_a[i].qid >= 0) pca_a[n++] = pca_a[i];
    pca_list.resize(n);
}

static void
x_update_chain(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name, int qdir)
{
    PoreCAlign* pa = chain.data();
    int pc = chain.size();
    for (auto& p : pca_list) {
        if (p.qdir != qdir) continue;
        //if (!p.lqm) continue;
        //if (!p.rqm) continue;
        //if (!p.is_perfect) continue;
        for (int i = 0; i < pc - 1; ++i) {
            if (pa[i].qid == -1 || pa[i+1].qid == -1) continue;
            if (p.qoff > pa[i].qoff) break;
            if (p.qoff == pa[i].qoff && p.qend == pa[i+1].qend) {
                double avg_pi = (pa[i].pi + pa[i+1].pi) / 2;
                bool r1 = (p.pi > avg_pi - 4.0) || p.pi >= 80.0;
                //bool r2 = (p.pi > avg_pi - 4.0) && ((!pa[i].lqm) || (!pa[i].rqm) || (!pa[i+1].lqm) || (!pa[i+1].rqm));
                if (r1) {
                    //HBN_LOG("%s replace", read_name); dump_chain_pca(fprintf, stderr, pa[i], i); dump_chain_pca(fprintf, stderr, pa[i+1], i+1);
                    //fprintf(stderr, "by\n"); dump_chain_pca(fprintf, stderr, p, -1); 
                    pa[i].qid = -1;
                    pa[i+1] = p;
                    break;
                }
            }
        }
    }
    int n = 0;
    for (int i = 0; i < pc; ++i) if (pa[i].qid >= 0) pa[n++] = pa[i];
    chain.resize(n);
}

static void
x_update_chain_p2(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    PoreCAlign* pa = chain.data();
    int pc = chain.size();
    for (auto& pca : pca_list) {
        int f = -1, t = -1;
        for (int i = 0; i < pc; ++i) {
            if (pa[i].chain_qoff == pca.chain_qoff) {
                f =i;
                break;
            }
        }
        if (f == -1) continue;
        for (int i = 0; i < pc; ++i) {
            if (pa[i].chain_qend == pca.chain_qend) {
                t = i;
                break;
            }
        }
        if (t == -1) continue;
        if (f >= t) continue;

        double avg_pi = 0;
        for (int i = f; i <= t; ++i) avg_pi += pa[i].pi;
        avg_pi /= (t - f + 1);
        if (pca.pi < avg_pi - 4.0 && pca.pi < 80.0) continue;
        //HBN_LOG("++++ merge, pi = %g, avg_pi = %g", pca.pi, avg_pi);
        //for (int i = f; i <=t ; ++i) dump_chain_pca(fprintf, stderr, pa[i], i);
        //fprintf(stderr, "into\n");
        //dump_chain_pca(fprintf, stderr, pca, -1);
        pa[f] = pca;
        for (int i = f + 1; i <= t; ++i) pa[i].qid = -1;
        int n = 0;
        for (int i = 0; i < pc; ++i) if (pa[i].qid >= 0) pa[n++] = pa[i];
        pc = n;
    }
    chain.resize(pc);
}

static void
x_update_chain_p3(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    PoreCAlign* pa = chain.data();
    int pc = chain.size();
    for (auto& pca : pca_list) {
        for (int i = 0; i < pc; ++i) {
            //if (pca.qoff == pa[i].qoff && pca.qend == pa[i].qend && pca.pi > pa[i].pi) {
	    if (pca.chain_qoff == pa[i].chain_qoff && pca.chain_qend == pa[i].chain_qend && pca.pi > pa[i].pi) {
                //HBN_LOG("**** %s better replace", read_name); dump_chain_pca(fprintf, stderr, pa[i], i);
                //dump_chain_pca(fprintf, stderr, pca, -1);
                pa[i] = pca;
                break;
            }
        }
    }
}

void
update_chain(std::vector<PoreCAlign>& chain, std::vector<PoreCAlign>& pca_list, const char* read_name)
{
    vector<PoreCAlign> fwd_pca_list, rev_pca_list;
    for (auto& pca : chain) {
        if (pca.qdir == FWD) {
            fwd_pca_list.push_back(pca);
        } else {
            rev_pca_list.push_back(pca);
        }
    }
    sort(fwd_pca_list.begin(), fwd_pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    sort(rev_pca_list.begin(), rev_pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    x_update_chain(fwd_pca_list, pca_list, read_name, FWD);
    x_update_chain(rev_pca_list, pca_list, read_name, REV);
    chain.clear();
    chain.insert(chain.end(), fwd_pca_list.begin(), fwd_pca_list.end());
    chain.insert(chain.end(), rev_pca_list.begin(), rev_pca_list.end());
    sort(chain.begin(), chain.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.chain_qoff < b.chain_qoff; });
    x_update_chain_p2(chain, pca_list, read_name);
    x_update_chain_p3(chain, pca_list, read_name);
}

void
smooth_pca_list(std::vector<PoreCAlign>& pca_list,
    EChainType& chain_type,
    std::vector<PoreCAlign>& all_pca_list,
    SeqReader* ref,
    const char* read_name,
    const u8* fwd_read,
    const u8* rev_read,
    const int* vdfa, 
    const int vdfc,
    const int enzyme_size,
    HbnTracebackData* tbck_data)
{
    vector<PoreCAlign> fwd_pca_list, rev_pca_list;
//fprintf(stderr, "======= merging pca list\n");
    for (auto& pca : pca_list) {
	//dump_chain_pca(fprintf, stderr, pca, -1);
        if (pca.qdir == FWD) {
            fwd_pca_list.push_back(pca);
        } else {
            rev_pca_list.push_back(pca);
        }
    }
    sort(fwd_pca_list.begin(), fwd_pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    sort(rev_pca_list.begin(), rev_pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.qoff < b.qoff; });
    s_merge_adjacent_pca(fwd_pca_list, ref, fwd_read, rev_read, enzyme_size, tbck_data);
    s_merge_adjacent_pca(rev_pca_list, ref, fwd_read, rev_read, enzyme_size, tbck_data);
    pca_list.clear();
    pca_list.insert(pca_list.end(), fwd_pca_list.begin(), fwd_pca_list.end());
    pca_list.insert(pca_list.end(), rev_pca_list.begin(), rev_pca_list.end());
    sort(pca_list.begin(), pca_list.end(), [](const PoreCAlign& a, const PoreCAlign& b) { return a.chain_qoff < b.chain_qoff; });

    update_chain(pca_list, all_pca_list, read_name);

    smooth_pca_list_pass2(pca_list, chain_type, ref, fwd_read, rev_read, vdfa, vdfc, enzyme_size, tbck_data);
}