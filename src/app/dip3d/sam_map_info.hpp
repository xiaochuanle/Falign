#ifndef __SAM_MAP_INFO_HPP
#define __SAM_MAP_INFO_HPP

#include "../../htslib/sam.h"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/seq_name2id_map.hpp"
#include "../../corelib/unpacked_seqdb.hpp"

#include <string>
#include <vector>

class SAM_MapInfo 
{
public:
    SAM_MapInfo() {}

    ~SAM_MapInfo() {}

    bool parse(sam_hdr_t* hdr, bam1_t* bam, HbnUnpackedDatabase* ref, SeqName2IdMap* ref_name2id);

public:

    const char* query_name() const {
        return M_query_name.c_str();
    }

    int query_strand() const {
        return M_query_strand;
    }

    int query_size() const {
        return M_query_size;
    }

    const char* fwd_seq() const {
        return M_fwd_seq.c_str();
    }

    const char* rev_seq() const {
        return M_rev_seq.c_str();
    }

    int fwd_pass() const {
        return M_fwd_pass;
    }

    int rev_pass() const {
        return M_rev_pass;
    }

    bool is_mapped() const {
        return M_is_mapped;
    }

    int qb() const {
        return M_qb;
    }

    int qe() const {
        return M_qe;
    }

    int sid() const {
        return M_sid;
    }

    const char* sname() const {
        return M_sname.c_str();
    }

    int sb() const {
        return M_sb;
    }

    int se() const {
        return M_se;
    }

    int ssize() const {
        return M_ss;
    }

    int mapQ() const {
        return M_mapQ;
    }

    double identity() const {
        return M_identity;
    }

    const char* qas() const {
        return M_qas.c_str();
    }

    const char* sas() const {
        return M_sas.c_str();
    }

    int as_size() const {
        return M_as_size;
    }

    const int* qas_pos_list() const {
        return M_qas_pos_list.data();
    }

    const int* sas_pos_list() const {
        return M_sas_pos_list.data();
    }

    void dump_map_info() const {
        if (!M_is_mapped) return;
        fprintf(stderr, "[%d, %d, %d, %d]", M_query_strand, M_qb, M_qe, M_query_size);
        fprintf(stderr, " x ");
        fprintf(stderr, "[%d, %d, %d, %d]", M_sid, M_sb, M_se, M_ss);
        fprintf(stderr, ", %d, %g", M_mapQ, M_identity);
        fprintf(stderr, ", %s", M_query_name.c_str());
        fprintf(stderr, ", %s", M_sname.c_str());
        fprintf(stderr, "\n");
    }

    void dump() const {
        fprintf(stderr, "query name:     %s\n", M_query_name.c_str());
        fprintf(stderr, "query strand:   %d\n", M_query_strand);
        fprintf(stderr, "query size:     %d\n", M_query_size);
        fprintf(stderr, "fwd-pass:       %d\n", M_fwd_pass);
        fprintf(stderr, "rev-pass:       %d\n", M_rev_pass);

        const char* query_seq = (M_query_strand == FWD) ? fwd_seq() : rev_seq();
        fprintf(stderr, "query:          ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%c", query_seq[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, "%c", query_seq[i]);
        fprintf(stderr, "\n");

        const char* rev_query_seq = (M_query_strand == FWD) ? rev_seq() : fwd_seq();
        fprintf(stderr, "query:          ");
        for (int i = 0; i < 10; ++i) fprintf(stderr, "%c", rev_query_seq[i]);
        fprintf(stderr, "...");
        for (int i = M_query_size - 10; i < M_query_size; ++i) fprintf(stderr, "%c", rev_query_seq[i]);
        fprintf(stderr, "\n");

        dump_map_info();
    }

private:

    void x_clear_map_info() {
        M_query_strand = F_R;
        M_query_size = 0;
        M_fwd_seq.clear();
        M_rev_seq.clear();
        M_fwd_pass = 0;
        M_rev_pass = 0;

        M_is_mapped = false;
        M_sid = -1;
        M_sname.clear();
        M_mapQ = 0;
        M_identity = 0.0;
        M_qas.clear();
        M_sas.clear();
        M_as_size = 0;
        M_qas_pos_list.clear();
        M_sas_pos_list.clear();
    }

private:
    std::vector<std::pair<char, int>>           M_cigar_op_list;

    std::string                                 M_query_name;
    int                                         M_query_strand;
    int                                         M_query_size;
    std::string                                 M_fwd_seq;
    std::string                                 M_rev_seq;
    int                                         M_fwd_pass;
    int                                         M_rev_pass;

    /// mapping info
    bool                                        M_is_mapped;
    int                                         M_qb;
    int                                         M_qe;
    int                                         M_sid;
    std::string                                 M_sname;
    int                                         M_sb;
    int                                         M_se;
    int                                         M_ss;
    int                                         M_mapQ;
    double                                      M_identity;
    std::string                                 M_qas;
    std::string                                 M_sas;
    int                                         M_as_size;
    std::vector<int>                            M_qas_pos_list;
    std::vector<int>                            M_sas_pos_list;
};

#endif // __SAM_MAP_INFO_HPP
