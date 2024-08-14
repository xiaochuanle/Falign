#include "vcf-reader.hpp"

#include "../../corelib/pdqsort.h"

using namespace std;

int get_phase_set_type(bcf_hdr_t* bcf_hdr)
{
    bcf_hrec_t* ps_hrec = bcf_hdr_get_hrec(bcf_hdr, BCF_HL_FMT, "ID", "PS", nullptr);
    if (!ps_hrec) HBN_ERR("Phase set format '##FORMAT=<ID=PS, ...>' is not present in VCF header");

    int type_key_idx = -1;
    for (int i = 0; i < ps_hrec->nkeys; ++i) {
        if (strcmp(ps_hrec->keys[i], "Type") == 0) {
            type_key_idx = i;
            break;
        }
    }
    if (type_key_idx == -1) HBN_ERR("Type information is missing in the phase set format '##FORMAT=<ID=PS, ...>' of VCF header");
    
    int phase_set_type = -1;
    if (strcmp(ps_hrec->vals[type_key_idx], "String") == 0) {
        phase_set_type = BCF_HT_STR;
    } else if (strcmp(ps_hrec->vals[type_key_idx], "Integer") == 0) {
        phase_set_type = BCF_HT_INT;
    } else {
        HBN_ERR("Unsupported phase set (PS) data type '%s'", ps_hrec->vals[type_key_idx]);
    }
    HBN_LOG("Phase set (PS) data type: %s", ps_hrec->vals[type_key_idx]);

    return phase_set_type;
}

void get_vcf_gt(bcf_hdr_t* bcf_hdr, bcf1_t* bcf, int gt[], int& is_phased)
{
    gt[0] = gt[1] = -1;
    is_phased = 0;

    const int num_samples = bcf_hdr_nsamples(bcf_hdr);
    int32_t* gt_arr = nullptr;
    int ngt_arr = 0;
    int ngt = bcf_get_genotypes(bcf_hdr, bcf, &gt_arr, &ngt_arr);
    if ( ngt<=0 ) return; // GT not present

    int max_ploidy = ngt / num_samples;
    int32_t* ptr = gt_arr + 0 * max_ploidy; // use 0-th sample
    for (int j = 0; j < 2 && j < max_ploidy; ++j) {
        // if true, the sample has smaller ploidy
        if ( ptr[j]==bcf_int32_vector_end ) break;

        // missing allele
        if ( bcf_gt_is_missing(ptr[j]) ) continue;

        // the VCF 0-based allele index
        gt[j] = bcf_gt_allele(ptr[j]);

        // is phased?
        if (bcf_gt_is_phased(ptr[j])) is_phased = 1;
    }
}

bool get_vcf_phase_set(bcf_hdr_t* bcf_hdr, bcf1_t* bcf,
    const int phase_set_value_type,
    std::string& s_phase_set,
    int& i_phase_set)
{
    bool status = false;
    if (phase_set_value_type == BCF_HT_INT) {
        int* arr = nullptr;
        int n = 0;
        if (bcf_get_format_int32(bcf_hdr, bcf, "PS", &arr, &n) > 0) {
            i_phase_set = arr[0];
            status = true;
        }
        free(arr);
    } else if (phase_set_value_type == BCF_HT_STR) {
        int ndst = 0; 
        char **dst = NULL;
        if (bcf_get_format_string(bcf_hdr, bcf, "PS", &dst, &ndst) > 0) {
            s_phase_set = dst[0];
            status = true;
        }
        free(dst[0]); free(dst);
    }

    return status;
}

bool parse_phased_snp_from_one_vcf(bcf_hdr_t* bcf_hdr,
    bcf1_t* bcf,
    const int phase_set_value_type,
    std::string& s_phase_set,
    int& i_phase_set,
    PhasedHetSNP& snp)
{
    bcf_unpack(bcf, BCF_UN_ALL);
    if (!bcf_is_snp(bcf)) return false;

    int REF = 'N';
    int ALT = 'N';
    if (bcf->n_allele == 2) {
        if (bcf->d.allele[0][1] == 0) {
            REF = bcf->d.allele[0][0];
        }
        if (bcf->d.allele[1][1] == 0) {
            ALT = bcf->d.allele[1][0];
        }
    }
    int r = (nst_nt4_table[REF] < 4) && (nst_nt4_table[ALT] < 4);
    if (!r) return false;

    int gt[2], is_phased = 0;
    get_vcf_gt(bcf_hdr, bcf, gt, is_phased);
    r = is_phased;
    if (r) r = (gt[0] == 0 && gt[1] == 1) || (gt[0] == 1 && gt[1] == 0);
    if (!r) return false;

    if (!get_vcf_phase_set(bcf_hdr, bcf, phase_set_value_type, s_phase_set, i_phase_set)) {
        return false;
    }

    snp.sid = bcf->rid;
    snp.soff = bcf->pos;
    snp.REF = REF;
    snp.ALT = ALT;
    snp.gt[0] = gt[0];
    snp.gt[1] = gt[1];
    if (phase_set_value_type == BCF_HT_INT) {
        snp.phase_set = i_phase_set;
    } else {
        snp.phase_set = -1;
    }
    return true;
}

void 
load_phased_het_snps_from_vcf_file(const char* vcf_path,
    const bool load_bcf_records,
    std::map<std::string, int>** p_phase_set_name2id,
    std::vector<std::string>** p_phase_set_names,
    bcf_hdr_t** p_bcf_hdr,
    std::vector<std::pair<PhasedHetSNP, bcf1_t*>>& phased_snps)
{
    vcfFile* vcfin = vcf_open(vcf_path, "r");
    hts_set_threads(vcfin, 8);
    bcf_hdr_t* vcfhdr = bcf_hdr_read(vcfin);
    if (!vcfhdr) HBN_ERR("FAIL at reading VCF header from %s", vcf_path);
    int phase_set_type = get_phase_set_type(vcfhdr);

    const int num_samples = bcf_hdr_nsamples(vcfhdr);
    if (num_samples== 0) {
        HBN_ERR("No sample is present in the VCF file");
    } 
    fprintf(stderr, "Found %d samples: %s", num_samples, vcfhdr->samples[0]);
    for (int i = 1; i < num_samples; ++i) fprintf(stderr, ", %s", vcfhdr->samples[i]);
    fprintf(stderr, "\n");
    if (num_samples > 1) fprintf(stderr, "We only use FORMAT of the first sample (%s)\n", vcfhdr->samples[0]);

    map<string, int>* phase_set_name2id = nullptr;
    vector<string>* phase_set_names = nullptr;
    if (phase_set_type == BCF_HT_STR) {
        phase_set_name2id = new map<string, int>();
        phase_set_names = new vector<string>();
        *p_phase_set_name2id = phase_set_name2id;
        *p_phase_set_names = phase_set_names;
    }
    string s_phase_set;
    int i_phase_set;

    while (1) {
        bcf1_t* bcf = bcf_init1();
        int r = vcf_read(vcfin, vcfhdr, bcf);
        if (r == -1) {
            bcf_destroy1(bcf);
            break;
        }
        if (r < 0) HBN_ERR("FAIL at reading VCF reacord");
        bcf_unpack(bcf, BCF_UN_ALL);

        PhasedHetSNP snp;
        if (parse_phased_snp_from_one_vcf(vcfhdr, bcf, phase_set_type, s_phase_set, i_phase_set, snp)) {
            if (phase_set_type == BCF_HT_STR) {
                auto pos = phase_set_name2id->find(s_phase_set);
                if (pos != phase_set_name2id->end()) {
                    snp.phase_set = pos->second;
                } else {
                    snp.phase_set = phase_set_name2id->size();
                    phase_set_name2id->insert(pair<string, int>(s_phase_set, snp.phase_set));
                    phase_set_names->push_back(s_phase_set);
                }
            } 
            if (load_bcf_records) {
                phased_snps.emplace_back(snp, bcf);
            } else {
                phased_snps.emplace_back(snp, nullptr);
                bcf_destroy1(bcf);
            }
        } else {
            bcf_destroy1(bcf);
        }
    }

    if (p_bcf_hdr) *p_bcf_hdr = vcfhdr; else bcf_hdr_destroy(vcfhdr);
    vcf_close(vcfin);

    if (p_phase_set_name2id) *p_phase_set_name2id = nullptr;
    if (p_phase_set_names) *p_phase_set_names = nullptr;
    if (phase_set_type == BCF_HT_STR) {
        if (p_phase_set_name2id) {
            *p_phase_set_name2id = phase_set_name2id;
        } else {
            delete phase_set_name2id;
        }
        if (p_phase_set_names) {
            *p_phase_set_names = phase_set_names;
        } else {
            delete phase_set_names;
        }
    }

    HBN_LOG("Load %zu phased het SNPs", phased_snps.size());
}

////////////////////////////

VcfPhasedHetSnpList::VcfPhasedHetSnpList(const char* vcf_path)
{
    M_vcf_path = vcf_path;

    vcfFile* vcfin = vcf_open(vcf_path, "r");
    hts_set_threads(vcfin, 8);
    M_bcf_hdr = bcf_hdr_read(vcfin);
    if (!M_bcf_hdr) HBN_LOG("FAIL at reading VCF header from %s", vcf_path);
    M_phase_set_value_type = get_phase_set_type(M_bcf_hdr);

    const int num_samples = bcf_hdr_nsamples(M_bcf_hdr);
    if (num_samples== 0) HBN_ERR("No sample is present in VCF file %s", vcf_path);
    fprintf(stderr, "Found %d samples: %s", num_samples, M_bcf_hdr->samples[0]);
    for (int i = 1; i < num_samples; ++i) fprintf(stderr, ", %s", M_bcf_hdr->samples[i]);
    fprintf(stderr, "\n");
    if (num_samples > 1) fprintf(stderr, "We will use only FORMAT of the first sample (%s)\n", M_bcf_hdr->samples[0]);

    string s_phase_set;
    int i_phase_set;
    while (1) {
        bcf1_t* bcf = bcf_init1();
        int r = vcf_read(vcfin, M_bcf_hdr, bcf);
        if (r == -1) {
            bcf_destroy1(bcf);
            break;
        }
        if (r < 0) HBN_ERR("FAIL at reading VCF reacord");
        bcf_unpack(bcf, BCF_UN_ALL);

        PhasedHetSNP snp;
        if (parse_phased_snp_from_one_vcf(M_bcf_hdr, bcf, M_phase_set_value_type, s_phase_set, i_phase_set, snp)) {
            if (M_phase_set_value_type == BCF_HT_STR) {
                auto pos = M_phase_set_name2id.find(s_phase_set);
                if (pos != M_phase_set_name2id.end()) {
                    snp.phase_set = pos->second;
                } else {
                    snp.phase_set = M_phase_set_name2id.size();
                    M_phase_set_name2id.insert(pair<string, int>(s_phase_set, snp.phase_set));
                    M_phase_set_names.push_back(s_phase_set);
                }
            }
            M_snp_list.push_back(snp);
        } 

        bcf_destroy1(bcf);
    }
    vcf_close(vcfin);

    pdqsort(M_snp_list.begin(), M_snp_list.end(),
        [](const PhasedHetSNP& x, const PhasedHetSNP& y) {
            return (x.sid < y.sid) || (x.sid == y.sid && x.soff < y.soff);
        });
    PhasedHetSNP* na = M_snp_list.data();
    size_t nc = M_snp_list.size();
    size_t i = 0;
    while (i < nc) {
        size_t j = i + 1;
        while (j < nc && na[i].sid == na[j].sid) ++j;

        M_chr_snps[na[i].sid] = pair<size_t, int>(i, j - i);

        i = j;
    }

    HBN_LOG("Load %zu phased het snps from %s", M_snp_list.size(), vcf_path);
}