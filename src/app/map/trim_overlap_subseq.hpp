#ifndef __TRIM_OVERLAP_SUBSEQ_H
#define __TRIM_OVERLAP_SUBSEQ_H

#include "hbn_options.hpp"
#include "chain_align_list.hpp"
#include "../../corelib/restrict_enzyme_loci_list.hpp"
#include "../../corelib/unpacked_seqdb.hpp"
#include "../../sw/hbn_traceback.hpp"

#include <string>
#include <vector>

 struct TrimPca {
    int is_valid;
    PoreCAlign pca;
    int qas_offset;
    int sas_offset;
    int as_size;
    int qdir;
    int qoff, qend;
    int fwd_qoff, fwd_qend;
    int e_qoff, e_qend;
    int sdir;
    int soff, send;
    int e_soff, e_send;
    int frag_id;
} ;

#define dump_tpca(output_func, stream, tpca) output_func(stream, "[%d, %d, %d, %d] x [%d, %d, %d, %d]\n", \
        tpca.pca.qid, tpca.qdir, tpca.qoff, tpca.qend,\
        tpca.pca.sid, tpca.sdir, tpca.soff, tpca.send)

static inline bool operator < (const TrimPca& lhs, const TrimPca& rhs)
{
        if (lhs.qoff >= rhs.qoff && lhs.qend <= rhs.qend) return true;
        if (lhs.pca.sid == rhs.pca.sid && lhs.soff >= rhs.soff && lhs.send <= rhs.send) return true;
        return false;
}

struct TrimPcaList
{
    std::vector<TrimPca> tpca_list;
    std::string align_strings;

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

    void clear() {
        tpca_list.clear();
        align_strings.clear();
    }
};

void
trim_overlap_subseqs(HbnTracebackData* tbck_data,
    RestrictEnzymeLociList* reloci_list,
    QueryVdfEndPointList* qvep_list,
    HbnUnpackedDatabase* subjects,
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    PoreCAlign* pca_a,
    int pca_c,
    const EChainType chain_type,
    TrimPcaList& trim_pca_list);

#endif // __TRIM_OVERLAP_SUBSEQ_H