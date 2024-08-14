#ifndef KSW2_WRAPPER_H
#define KSW2_WRAPPER_H

#include "../corelib/hbn_aux.h"

#include <string>
#include <vector>

typedef struct {
    void* km;
    int reward, penalty, ambi_penalty;
    int go, ge;
    int go1, ge1;
    int zdrop;
    int band_width;
    int end_bonus;
    int8_t mat[25];
    int score_param_is_set;
    std::vector<u8> qfrag;
    std::vector<u8> tfrag;
    int query_id;
} Ksw2Data;

void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b, int8_t sc_ambi);

void extract_sw_scoring_params(int* match_reward, int* mismatch_penalty, int* gap_open, int* gap_extend, int* gap_open1, int* gap_extend1);

void
ksw2_set_params(Ksw2Data* data,
    int reward,
    int penalty,
    int ambi,
    int go,
    int ge,
    int zdrop,
    int band_width);

void
ksw2_extd2_set_params(Ksw2Data* data);

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
        std::string& saln);

#endif // KSW2_WRAPPER_H