#ifndef __EXTEND_HIT_LIST_HPP
#define __EXTEND_HIT_LIST_HPP

#include "../../corelib/gapped_candidate.h"
#include "../../sw/hbn_traceback.hpp"
#include "../../sw/hbn_traceback_aux.h"
#include "chain_align_list.hpp"
#include "hbn_options.hpp"

#include <string>
#include <vector>

struct frag_align_struct {
    BOOL is_valid;
    int qdir, qoff, qend, qsize;
    int sid, soff, send, ssize;
    int score, chain_score;
    double pi;
    int qas_offset;
    int sas_offset;
    int as_size;

    int bqoff, bqend;
    int bsoff, bsend;

    int sqoff, sqend;
    int ssoff, ssend;

    int left_align_id, right_align_id;
};

#define dump_frag_align(output_func__, stream__, align__, idx__) do { \
    if ((idx__) != -1) output_func__(stream__, "%d\t", idx__); \
    output_func__(stream__, "[%d, %d, %d, %d] x [%d, %d, %d], %d, %g\n", \
        (align__).qdir, \
        (align__).qoff, \
        (align__).qend, \
        (align__).qsize, \
        (align__).soff, \
        (align__).send, \
        (align__).ssize, \
        (align__).score, \
        (align__).pi); \
} while(0)

struct frag_align_list_struct {
    std::vector<frag_align_struct> frag_align_list;
    std::string align_strings;
};

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
    std::vector<PoreCAlign>& pca_list);

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
    std::vector<PoreCAlign>& align_list);

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
    HbnTracebackData* tbck_data);

#endif // __EXTEND_HIT_LIST_HPP