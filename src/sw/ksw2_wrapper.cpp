#include "ksw2_wrapper.hpp"

#include "ksw2.h"
#include "hbn_traceback_aux.h"

#include <cstring>

void extract_sw_scoring_params(int* match_reward, int* mismatch_penalty, int* gap_open, int* gap_extend, int* gap_open1, int* gap_extend1)
{
#if 0
    data->reward = 2;
    data->penalty = 4;
    data->go = 4;
    data->ge = 2;
    data->go1 = 24;
    data->ge1 = 1;
#elif 0
    //-o 5 -O 56 -e 4 -E 1 -A 2 -B 5
    data->reward = 2;
    data->penalty = 5;
    data->go = 5;
    data->ge = 4;
    data->go1 = 56;
    data->ge1 = 1;
#endif 
    if (match_reward) *match_reward = 2;
    if (mismatch_penalty) *mismatch_penalty = 4;
    if (gap_open) *gap_open = 4;
    if (gap_extend) *gap_extend = 2;
    if (gap_open1) *gap_open1 = 24;
    if (gap_extend1) *gap_extend1 = 1;
}

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi)
{
	int i, j;
	a = a < 0? -a : a;
	b = b > 0? -b : b;
	sc_ambi = sc_ambi > 0? -sc_ambi : sc_ambi;
	for (i = 0; i < m - 1; ++i) {
		for (j = 0; j < m - 1; ++j)
			mat[i * m + j] = i == j? a : b;
		mat[i * m + m - 1] = sc_ambi;
	}
	for (j = 0; j < m; ++j)
		mat[(m - 1) * m + j] = sc_ambi;
}

void
ksw2_set_params(Ksw2Data* data,
    int reward,
    int penalty,
    int ambi,
    int go,
    int ge,
    int zdrop,
    int band_width)
{
    data->reward = reward;
    data->penalty = penalty;
    data->ambi_penalty = ambi;
    data->go = go;
    data->ge = ge;
    data->go1 = 0;
    data->ge1 = 0;
    data->zdrop = zdrop; // 400;
    data->band_width = band_width; // 2048;
    data->end_bonus = -1;
    ksw_gen_simple_mat(5, data->mat, reward, penalty, ambi);
    data->score_param_is_set = 1;
}

void
ksw2_extd2_set_params(Ksw2Data* data)
{
#if 1
    data->reward = 2;
    data->penalty = 4;
    data->go = 4;
    data->ge = 2;
    data->go1 = 24;
    data->ge1 = 1;
#else
    //-o 5 -O 56 -e 4 -E 1 -A 2 -B 5
    data->reward = 2;
    data->penalty = 5;
    data->go = 5;
    data->ge = 4;
    data->go1 = 56;
    data->ge1 = 1;
#endif 
    data->ambi_penalty = 1;
    data->zdrop = -1;
    data->end_bonus = -1;
    data->band_width = -1;
    ksw_gen_simple_mat(5, data->mat, data->reward, data->penalty, data->ambi_penalty);
}

int nw_ksw2_extd2(Ksw2Data* data,
        const int qid,
        const u8* query,
        const int qfrom,
        const int qto,
        const int qsize,
        const int sid,
        const u8* subject,
        const int sfrom,
        const int sto,
        const int ssize,
        const int min_align_size,
        const double min_ident_perc,
        int max_distance,
        int* qoff,
        int* qend,
        int* soff,
        int* send,
        double* ident_perc,
	    std::string& qaln,
	    std::string& saln)
{
	hbn_assert(qfrom < qto, "qid = %d, [%d, %d, %d] x [%d, %d, %d]", data->query_id, qfrom, qto, qsize, sfrom, sto, ssize);
	hbn_assert(qto <= qsize, "qid = %d, [%d, %d, %d] x [%d, %d, %d]", data->query_id, qfrom, qto, qsize, sfrom, sto, ssize);
	hbn_assert(sfrom < sto, "qid = %d, [%d, %d, %d] x [%d, %d, %d]", data->query_id, qfrom, qto, qsize, sfrom, sto, ssize);
	hbn_assert(sto <= ssize, "qid = %d, [%d, %d, %d] x [%d, %d, %d]", data->query_id, qfrom, qto, qsize, sfrom, sto, ssize);
    ksw_extz_t ez; memset(&ez, 0, sizeof(ksw_extz_t));
    int flag = 0;
    if (max_distance == 0) max_distance = hbn_max(qto - qfrom, sto - sfrom) * 0.1;
    //int max_band_width = hbn_min(max_distance, KSW_MAX_BW);
	int max_band_width = max_distance;
    int zdrop = -1;
    //HBN_LOG("distance = %d, band_width = %d", max_distance, max_band_width);
    const u8* qsubseq = query + qfrom;
    const int qsubseq_size = qto - qfrom;
    const u8* ssubseq = subject + sfrom;
    const int ssubseq_size = sto - sfrom;
    //ksw_extd2_sse(data->km, qsubseq_size, qsubseq, ssubseq_size, ssubseq, 5, data->mat,
    //    data->go, data->ge, data->go1, data->ge1, max_band_width, zdrop, data->end_bonus, flag, &ez);
    ksw_extz2_sse(data->km, qsubseq_size, qsubseq, ssubseq_size, ssubseq, 5, data->mat,
        data->go, data->ge, max_band_width, zdrop, data->end_bonus, flag, &ez);
    if (ez.n_cigar == 0) {
        //HBN_LOG("[%d, %d, %d, %d] x [%d, %d, %d, %d] %g%% sw fail, dist = %d", 
        //    qid, *qoff, *qend, qsize, sid, *soff, *send, ssize, *ident_perc, max_distance);
        if (ez.cigar) kfree(data->km, ez.cigar);
        return 0;
    }

    int qi = 0;
    int si = 0;
    qaln.clear();
    saln.clear();
    for (int k = 0; k < ez.n_cigar; ++k) {
        int op_num = ez.cigar[k]>>4;
        int op_type = "MIDN"[ez.cigar[k]&0xf];
        int c = 0;
        switch (op_type)
        {
        case 'M':
            for (int t = 0; t < op_num; ++t, ++qi, ++si) {
                c = qsubseq[qi];
                hbn_assert(c >= 0 && c < 4, "query_id: %d, qi: %d, si: %d, c: %d [%d, %d, %d] x [%d, %d, %d]", data->query_id, qi, si, c, qfrom, qto, qsize, sfrom, sto, ssize);
                c = DECODE_RESIDUE(c);
                qaln += c;
                c = ssubseq[si];
                hbn_assert(c >= 0 && c < 4, "query_id: %d, qi: %d, si: %d, c: %d [%d, %d, %d] x [%d, %d, %d]", data->query_id, qi, si, c, qfrom, qto, qsize, sfrom, sto, ssize);
                c = DECODE_RESIDUE(c);
                saln += c;
            }
            break;
        case 'I':
            for (int t = 0; t < op_num; ++t, ++qi) {
                c = qsubseq[qi];
                hbn_assert(c >= 0 && c < 4);
                c = DECODE_RESIDUE(c);
                qaln += c;
                saln += GAP_CHAR;
            }
            break;
        case 'D':
            for (int t = 0; t < op_num; ++t, ++si) {
                qaln += GAP_CHAR;
                c = ssubseq[si];
                hbn_assert(c >= 0 && c < 4);
                c = DECODE_RESIDUE(c);
                saln += c;
            }
            break;            
        default:
            HBN_LOG("invalid op_type: %d", op_type);
            break;
        }
    }
    if (qi < qsubseq_size || si < ssubseq_size) {
        //HBN_LOG("[%d, %d, %d, %d] x [%d, %d, %d, %d] %g%% sw fail, dist = %d, qi = %d, si = %d", 
        //    qid, *qoff, *qend, qsize, sid, *soff, *send, ssize, *ident_perc, max_distance, qi, si);
        if (ez.cigar) kfree(data->km, ez.cigar);
        return 0;
    }
    hbn_assert(qi == qsubseq_size);
    hbn_assert(si == ssubseq_size);
    kfree(data->km, ez.cigar);

    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, qsubseq, 0, qsubseq_size, qaln.c_str(),
        0, ssubseq, 0, ssubseq_size, saln.c_str(),
        qaln.size(), TRUE);
    *ident_perc = calc_ident_perc(qaln.c_str(), saln.c_str(), qaln.size(), NULL, NULL);
    if (*ident_perc < min_ident_perc) {
        HBN_LOG("identity is too low: %g", *ident_perc);
        return 0;
    }
    *qoff = qfrom;
    *qend = qto;
    *soff = sfrom;
    *send = sto;
    validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        0, query, *qoff, *qend, qaln.c_str(),
        0, subject, *soff, *send, saln.c_str(),
        qaln.size(), TRUE);

    return 1;    
}
