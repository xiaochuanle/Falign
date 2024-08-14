#ifndef __HBN_TRACE_BACK_HPP
#define __HBN_TRACE_BACK_HPP

#include "dalign.hpp"
#include "edlib_wrapper.hpp"
#include "ksw2.h"
#include "ksw2_wrapper.hpp"
#include "small_edlib_align.hpp"

#include <string>
#include <vector>

typedef struct {
    int qoff;
    int soff;
    int length;
} DDFS_Seed;

struct HbnTracebackData {
    int qoff, qend, qsize;
    int soff, send, ssize;
    int score;
    int dist;
    double ident_perc;

    std::vector<char> qabuf;
    std::vector<char> sabuf;
    std::string ext_qabuf;
    std::string ext_sabuf;
    std::vector<u8> qfrag;
    std::vector<u8> sfrag;
    std::vector<DDFS_Seed> trace_seeds;

    char* qas;
    char* qae;
    char* sas;
    char* sae;

    DalignData* dalign;
    EdlibAlignData* edlib;
    Ksw2Data* ksw;
    small_edlib_align_struct* small_edlib_data;

    void dump(FILE* out = stderr) {
        fprintf(out, "[%d, %d, %d] x [%d, %d, %d], %g\n",
            qoff, qend, qsize, soff, send, ssize, ident_perc);
    }

    HbnTracebackData() {
        dalign = new DalignData(0.35);
        edlib = new EdlibAlignData();
        ksw = new Ksw2Data();
        ksw2_extd2_set_params(ksw);
        const int kMaxOverHang = 3000;
        ksw->band_width = kMaxOverHang * 0.2;
        ksw->zdrop = 40;
        small_edlib_data = small_edlib_align_struct_new();
    }

    ~HbnTracebackData() {
        delete dalign;
        small_edlib_align_struct_free(small_edlib_data);
        delete edlib;
        delete ksw;
    }

    void init_ul(int qoff, int qsize, int soff, int ssize) {
        int Lb = 2 * std::min(qoff, soff) + 50000;
        int Rb = 2 * std::min(qsize - qoff, ssize - soff) + 50000;
        int Bs = Lb + Rb;
        hbn_assert(qoff >= 0 && qoff < qsize, "qoff = %d, qsize = %d, soff = %d, ssize = %d, Lb = %d, Rb = %d",
            qoff, qsize, soff, ssize, Lb, Rb);
        hbn_assert(soff >= 0 && soff < ssize, "qoff = %d, qsize = %d, soff = %d, ssize = %d, Lb = %d, Rb = %d",
            qoff, qsize, soff, ssize, Lb, Rb);
        qabuf.resize(Bs);
        sabuf.resize(Bs);
        qas = qae = qabuf.data() + Lb;
        sas = sae = sabuf.data() + Lb;

#if 0
        int wrk_l = 4 * hbn_max(qsize, ssize);
        qabuf.resize(wrk_l);
        sabuf.resize(wrk_l);
        int wrk_ll = 2 * hbn_max(qsize, ssize);
        qas = qae = qabuf.data() + wrk_ll;
        sas = sae = sabuf.data() + wrk_ll;
#endif

        *qae = '\0';
        *sae = '\0';
        qoff = qend = 0;
        soff = send = 0;
        dist = 0;
        ident_perc = 0.0;
    }

    void init(int qoff, int qsize, int soff, int ssize) {
        int Lb = 2 * std::min(qoff, soff) + 5000;
        int Rb = 2 * std::min(qsize - qoff, ssize - soff) + 5000;
        int Bs = Lb + Rb;
        hbn_assert(qoff >= 0 && qoff < qsize, "qoff = %d, qsize = %d, soff = %d, ssize = %d, Lb = %d, Rb = %d",
            qoff, qsize, soff, ssize, Lb, Rb);
        hbn_assert(soff >= 0 && soff < ssize, "qoff = %d, qsize = %d, soff = %d, ssize = %d, Lb = %d, Rb = %d",
            qoff, qsize, soff, ssize, Lb, Rb);
        qabuf.resize(Bs);
        sabuf.resize(Bs);
        qas = qae = qabuf.data() + Lb;
        sas = sae = sabuf.data() + Lb;

#if 0
        int wrk_l = 4 * hbn_max(qsize, ssize);
        qabuf.resize(wrk_l);
        sabuf.resize(wrk_l);
        int wrk_ll = 2 * hbn_max(qsize, ssize);
        qas = qae = qabuf.data() + wrk_ll;
        sas = sae = sabuf.data() + wrk_ll;
#endif

        *qae = '\0';
        *sae = '\0';
        qoff = qend = 0;
        soff = send = 0;
        dist = 0;
        ident_perc = 0.0;
    }
};

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
    const int process_over_hang);

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
    const char** sae_);

void
validate_mem(HBN_LOG_PARAMS_GENERIC,
    const u8* read, 
    const u8* subject,
    const DDFS_Seed* cdpsa,
    const int cdpsc);

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
    const double min_ident_perc);

#endif // __HBN_TRACE_BACK_HPP