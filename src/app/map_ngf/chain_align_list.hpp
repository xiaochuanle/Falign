#ifndef __CHAIN_ALIGN_LIST_H_
#define __CHAIN_ALIGN_LIST_H_

#include "../../algo/hbn_traceback.h"
#include "../../corelib/restrict_enzyme_loci_list.h"
#include "../../corelib/small_object_alloc.h"

#include <vector>

#define PC_MAX_ALIGN_OVLP 20
#define kMaxMissedCovBase 200

typedef struct {
    int qid, qdir, qoff, qend, qsize;
    int qc;
    int chain_qoff, chain_qend;
    int enzyme_qoff, enzyme_qend;

    int sid, soff, send, ssize;
    int enzyme_soff, enzyme_send;

    bool lqm, rqm, lsm, rsm;
    bool lp, lpp, rp, rpp;
    bool lcp, lcpp, rcp, rcpp;
    bool is_perfect;
    bool is_pseudo_perfect;
    bool is_homologous;

    int map_q;
    double pi;
    int score;
    int chain_score;
} PoreCAlign;

#define dump_pca(output_func, stream, __pca, __idx) do { \
    if ((__idx) >= 0) output_func(stream, "%d\t", __idx); \
    output_func(stream, "[%d, %d, %d(%d), %d(%d)]", \
        (__pca).qid, (__pca).qdir, \
        (__pca).qoff, (__pca).lpp, \
        (__pca).qend, (__pca).rpp); \
    output_func(stream, " x "); \
    output_func(stream, "[%d, %d, %d]", \
        (__pca).sid, (__pca).soff, (__pca).send); \
    output_func(stream, " q_enzyme = [%d, %d]", \
        (__pca).enzyme_qoff, (__pca).enzyme_qend); \
    output_func(stream, " s_enzyme = [%d, %d]", \
        (__pca).enzyme_soff, (__pca).enzyme_send); \
    output_func(stream, ", score = %d", (__pca).score); \
    output_func(stream, ", %g\n", (__pca).pi); \
} while(0)

#define dump_chain_pca(output_func, stream, __pca, __idx) do { \
    if ((__idx) >= 0) output_func(stream, "%d\t", __idx); \
    output_func(stream, "[%d, %d, %d(%d), %d(%d)]", \
        (__pca).qid, (__pca).qdir, \
        (__pca).chain_qoff, (__pca).lcpp, \
        (__pca).chain_qend, (__pca).rcpp); \
    output_func(stream, " x "); \
    output_func(stream, "[%d, %d, %d]", \
        (__pca).sid, (__pca).soff, (__pca).send); \
    output_func(stream, " s_enzyme = [%d, %d]", \
        (__pca).enzyme_soff, (__pca).enzyme_send); \
    output_func(stream, ", score = %d", (__pca).score); \
    output_func(stream, ", %g\n", (__pca).pi); \
} while(0)

BOOL set_pca_chain_offset(PoreCAlign* pca, const int enzyme_size);

struct PoreCAlignChain {
    int is_valid;
    int sid;
    int q_cov;
    double avg_pi;
    int is_complete_chain;
    int is_pseudo_perfect_chain;
    int is_perfect_chain;
    std::vector<int> pca_a_idx_list;
    std::vector<PoreCAlign> chain;

    void clear() {
        is_valid = 1;
        sid = 0;
        q_cov = 0;
        avg_pi = 0.0;
        is_complete_chain = 0;
        is_perfect_chain = 0;
        is_pseudo_perfect_chain = 0;
        pca_a_idx_list.clear();
        chain.clear();
    }
};

struct PoreCAlignList{
    PoreCAlignChain best_chain;
    std::vector<PoreCAlign> best_pca_chain;
    std::vector<PoreCAlign> pca_list;
};

struct pca_chain_dp_point {
    int pca_a_idx;
    int q_cov;
    size_t s_dist;
    int score;
    void* parent;
};

struct pca_support {
    int has_precessor;
    int has_successor;
    int succ_offset;
    int succ_cnt;
    int prec_offset;
    int prec_cnt;
} ;

struct pca_supports {
    std::vector<pca_support> supports;
    std::vector<int> succ_idx_list;
} ;

struct PoreCAlignChainData {
    pca_supports supports;
    std::vector<pca_chain_dp_point> dp_stack;
    SmallObjectAlloc* dp_point_soa;
    std::vector<PoreCAlignList> subject_align_list_array;

    PoreCAlignChainData() {
        dp_point_soa = SmallObjectAllocNew(sizeof(pca_chain_dp_point));
    }

    ~PoreCAlignChainData() {
        if (dp_point_soa) SmallObjectAllocFree(dp_point_soa);
    }

    void release_soa() {
        dp_point_soa = SmallObjectAllocFree(dp_point_soa);
        dp_point_soa = SmallObjectAllocNew(sizeof(pca_chain_dp_point));
    }
};

typedef enum {
    ePerfectChain,
    ePseudoPerfectChain,
    eCompleteChain,
    eMaxCovChain
} EChainType;

const char* get_chain_type_name(EChainType type);

EChainType pca_chain_type(const int* vdfa, const int vdfc, PoreCAlign* pca_a, int pca_c);

void remove_duplicate_pca(PoreCAlign* pca_a, int* pca_c_);

void remove_contained_pca(PoreCAlign* pca_a, int* _pca_c);
void remove_contained_pca(std::vector<PoreCAlign>& pca_list);

double 
compute_pca_list_avg_pi(PoreCAlign* pca_a, int pca_c);

bool 
select_pca_chain(PoreCAlignChainData* pca_chain_data,
    const int* vdfa, 
    const int vdfc, 
    PoreCAlign* pca_array, 
    int pca_count, 
    const EChainType chain_type, 
    std::vector<PoreCAlign>& chain, 
    int* cov);

#endif // __CHAIN_ALIGN_LIST_H_