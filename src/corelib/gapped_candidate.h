#ifndef __GAPPED_CANDIDATE_H
#define __GAPPED_CANDIDATE_H

#include "hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int qid;
    int qdir;
    int qbeg;
    int qend;
    int qsize;
    int sid;
    int sdir;
    size_t sbeg;
    size_t send;
    size_t ssize;
    int score;
    BOOL is_hom;
} HbnInitHit;

typedef kvec_t(HbnInitHit) vec_init_hit;

#define dump_init_hit(output_func, stream, hit) \
    output_func(stream, "[%d, %d, %d, %d] x [%d, %d, %d], %d\n", \
        (hit).qid, \
        (hit).qdir, \
        (hit).qbeg, \
        (hit).qend, \
        (hit).sid, \
        (hit).sbeg, \
        (hit).send, \
        (hit).score)

void ks_introsort_init_hit_score_gt(size_t n, HbnInitHit* a);
void ks_introsort_init_hit_qid_lt(size_t n, HbnInitHit* a);
void ks_introsort_init_hit_sid_lt(size_t n, HbnInitHit* a);

#ifdef __cplusplus
}
#endif

#endif // __GAPPED_CANDIDATE_H