#include "sam_map_info.hpp"

#include "../../sw/hbn_traceback_aux.h"

using namespace std;

#if 0
CIGAR operation

OP      Description                                                 consume-query           consume-reference

M       alignment match (can be match or mismatch)                  yes                     yes
I       insertion to the reference                                  yes                     no 
D       deletion from the reference                                 no                      yes 
N       skipped region from the reference                           no                      yes
S       soft clipping (clipped sequences represent in SEQ)          yes                     no 
H       hard clipping (clipped sequences NOT represent in SEQ)      no                      no 
P       padding (silent deletion from padded reference)             no                      no 
=       sequence match                                              yes                     yes
X       sequence mismatch                                           yes                     yes
#endif 

static void cigar_to_alignment(const char* query,
    const int query_size,
    const u8* subject,
    const int subject_size,
    vector<pair<char, int>>& op_list,
    int& qb,
    int& qe,
    int& sb,
    int& se,
    std::string& qas,
    std::string& sas,
    std::vector<int>& q_pos_list,
    std::vector<int>& s_pos_list)
{
    pair<char, int>* opa = op_list.data();
    int opc = op_list.size();
    int opi = 0;
    qb = 0;
    sb = 0;
    if (opa[0].first == 'S') {
        qb = opa[0].second;
        opi = 1;
    } else if (opa[0].first == 'H') {
        opi = 1;
    }
    int qi = qb - 1;
    int si = sb - 1;

    qas.clear();
    sas.clear();
    q_pos_list.clear();
    s_pos_list.clear();
    for (; opi < opc; ++opi) {
        char op = opa[opi].first;
        int num = opa[opi].second;

        if (op == 'M') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'I') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                qas += query[qi];
                sas += GAP_CHAR;
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'D') {
            for (int i = 0; i < num; ++i) {
                ++si;
                qas += GAP_CHAR;
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else if (op == 'N') {
            for (int i = 0; i < num; ++i) {
                ++si;
                qas += GAP_CHAR;
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }            
        } else if (op == 'S') {
            continue;
        } else if (op == 'H') {
            continue;
        } else if (op == 'P') {
            continue;
        } else if (op == '=') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += DECODE_RESIDUE(subject[si]);
		        int qc = query[qi];
		        int sc = DECODE_RESIDUE(subject[si]);
		        if (qc != sc) {
			        fprintf(stderr, "qb = %d, sb = %d, qi = %d, si = %d, qc = %c, sc = %c\n", qb, sb, qi, si, qc, sc);
			        //cerr << NStr::CTempString(cigar, cigar_size) << '\n';
		        }
		        //hbn_assert(qc == sc);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }            
        } else if (op == 'X') {
            for (int i = 0; i < num; ++i) {
                ++qi;
                ++si;
                qas += query[qi];
                sas += DECODE_RESIDUE(subject[si]);
                q_pos_list.push_back(qi);
                s_pos_list.push_back(si);
            }
        } else {
            fprintf(stderr, "ERROR: Unrecognised CIGAR operation '%c' in\n", op);
            //cerr << NStr::CTempString(cigar, cigar_size) << '\n';
            exit(1);
        }
    }
    qe = qi;
    se = si;
}

static inline char s_decode_bam_query_base(int c)
{
    // 1 -> A
    // 2 -> C
    // 4 -> G
    // 8 -> T
    // 15 -> N
    switch (c) {
        case 1:
            return 'A';
        case 2:
            return 'C';
        case 4:
            return 'G';
        case 8:
            return 'T';
        case 15:
            return 'N';
        default:
            HBN_ERR("Illegal BAM base encoded value %s", c);
    }
}

static inline char completement_residue(const char c)
{
    switch (c) {
        case 'a':
        case 'A':
            return 'T';
        case 'c':
        case 'C':
            return 'G';
        case 'g':
        case 'G':
            return 'C';
        case 't':
        case 'T':
            return 'A';
        case 'n':
        case 'N':
            return 'N';
        default:
            return 'A';
    }
}

bool SAM_MapInfo::parse(sam_hdr_t* hdr, bam1_t* bam, HbnUnpackedDatabase* ref, SeqName2IdMap* ref_name2id)
{
    x_clear_map_info();
    if (ref == nullptr || ref_name2id == nullptr || bam->core.n_cigar == 0) return false;

    M_query_name = bam_get_qname(bam);
    M_query_strand = bam_is_rev(bam) ? REV : FWD;
    M_query_size = bam->core.l_qseq;

    if (M_query_strand == FWD) {
        uint8_t* seq = bam_get_seq(bam);
        uint8_t* qual = bam_get_qual(bam);
        for (int i = 0; i < M_query_size; ++i) {
            int c = bam_seqi(seq, i);
            c = s_decode_bam_query_base(c);
            int rc = completement_residue(c);
            M_fwd_seq += c;
            M_rev_seq += rc;
        }
        reverse(M_rev_seq.begin(), M_rev_seq.end());
    } else {
        uint8_t* seq = bam_get_seq(bam);
        uint8_t* qual = bam_get_qual(bam);
        for (int i = 0; i < M_query_size; ++i) {
            int rc = bam_seqi(seq, i);
            rc = s_decode_bam_query_base(rc);
            int c = completement_residue(rc);
            M_fwd_seq += c;
            M_rev_seq += rc;
        }        
        reverse(M_fwd_seq.begin(), M_fwd_seq.end());
    }

    M_is_mapped = true;
    uint32_t* cigar = bam_get_cigar(bam);
    M_cigar_op_list.clear();
    for (int i = 0; i < bam->core.n_cigar; ++i) {
        char op = bam_cigar_opchr(cigar[i]);
        int num = bam_cigar_oplen(cigar[i]);
        //fprintf(stderr, "%d\t%c\t%d\n", i, op, num);
        M_cigar_op_list.emplace_back(op, num);
    }
    //fprintf(stderr, "sid = %d\n", bam->core.tid);
    M_sname = hdr->target_name[bam->core.tid];
    M_sid = ref_name2id->GetIdFromNameSafe(M_sname);
    const u8* chr_seq = ref->GetSequence(M_sid);
    M_ss = ref->SeqSize(M_sid);
    const char* query_seq = (M_query_strand == FWD) ? fwd_seq() : rev_seq();
    cigar_to_alignment(query_seq, M_query_size, chr_seq + bam->core.pos, M_ss, M_cigar_op_list, M_qb, M_qe, M_sb, M_se, M_qas, M_sas, M_qas_pos_list, M_sas_pos_list);
    M_as_size = M_qas.size();
    M_sb += bam->core.pos;
    M_se += bam->core.pos;
    for (auto& p : M_sas_pos_list) p += bam->core.pos;
    ++M_qe;
    ++M_se;

    M_mapQ = bam->core.qual;
    M_identity = calc_effective_ident_perc(M_qas.c_str(), M_sas.c_str(), M_as_size);

    int qi = M_qb - 1, si = M_sb - 1;
    for (int i = 0; i < M_as_size; ++i) {
        if (M_qas[i] != GAP_CHAR) {
            int qc = M_qas[i];
            ++qi;
            hbn_assert(qi == M_qas_pos_list[i]);
            int qc1 = query_seq[qi];
            hbn_assert(qc == qc1);
        }
        if (M_sas[i] != GAP_CHAR) {
            int sc = M_sas[i];
            ++si;
            hbn_assert(si == M_sas_pos_list[i]);
            int sc1 = chr_seq[si];
            sc1 = DECODE_RESIDUE(sc1);
            hbn_assert(sc == sc1);
        }
    }

    return true;
}