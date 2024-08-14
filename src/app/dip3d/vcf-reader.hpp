#ifndef __VCF_READER_HPP
#define __VCF_READER_HPP

#include "../../corelib/seq_name2id_map.hpp"
#include "../../htslib/vcf.h"

#include <string>
#include <map>
#include <vector>

struct PhasedHetSNP
{
    int sid;
    int soff;
    char REF;
    char ALT;
    int phase_set;
    int gt[2];

    void dump() {
        fprintf(stderr, "sid = %d, soff = %d, REF = %c, ALT = %c, phase-set: %d, gt = [%d, %d]\n",
            sid, soff, REF, ALT, phase_set, gt[0], gt[1]);
    }
};

int get_phase_set_type(bcf_hdr_t* bcf_hdr);

void get_vcf_gt(bcf_hdr_t* bcf_hdr, bcf1_t* bcf, int gt[], int& is_phased);

void 
load_phased_het_snps_from_vcf_file(const char* vcf_path,
    const bool load_bcf_records,
    std::map<std::string, int>** p_phase_set_name2id,
    std::vector<std::string>** p_phase_set_names,
    bcf_hdr_t** p_bcf_hdr,
    std::vector<std::pair<PhasedHetSNP, bcf1_t*>>& phased_snps);

class VcfPhasedHetSnpList
{
public:
    VcfPhasedHetSnpList(const char* vcf_path);

    ~VcfPhasedHetSnpList() {
        bcf_hdr_destroy(M_bcf_hdr);
    }

    bool get_snp_idx(const char* chr_name, const int offset, PhasedHetSNP** snp_list, int* snp_list_idx, int* snp_list_size) {
        int sid = bcf_hdr_name2id(M_bcf_hdr, chr_name);
        if (sid == -1) return false;
        auto pos = M_chr_snps.find(sid);
        if (pos == M_chr_snps.end()) return false;
        PhasedHetSNP* na = M_snp_list.data() + pos->second.first;
        int nc = pos->second.second;

        if (na[0].soff >= offset) {
            if (snp_list) *snp_list = na;
            if (snp_list_idx) *snp_list_idx = 0;
            if (snp_list_size) *snp_list_size = nc;
            return true;
        }
        if (na[nc-1].soff < offset) return false;

        int L = 0, R = nc - 1;
        while (R - L > 1) {
            int M = (L + R) / 2;
            if (na[M].soff >= offset) {
                R = M;
            } else {
                L = M;
            }
        }
        hbn_assert(R - L == 1);
        hbn_assert(na[L].soff < offset && offset <= na[R].soff);
        if (snp_list) *snp_list = na;
        if (snp_list_idx) *snp_list_idx = R;
        if (snp_list_size) *snp_list_size = nc;
        return true;
    }

private:
    std::string                 M_vcf_path;
    bcf_hdr_t*                  M_bcf_hdr;
    std::map<int, std::pair<size_t, int>>
                                M_chr_snps;
    std::vector<PhasedHetSNP>   M_snp_list;

    int                         M_phase_set_value_type;
    std::map<std::string, int>  M_phase_set_name2id;
    std::vector<std::string>    M_phase_set_names;
};

#endif // __VCF_READER_HPP