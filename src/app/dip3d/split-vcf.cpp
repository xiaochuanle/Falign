#include "../../corelib/hbn_aux.h"
#include "../../htslib/vcf.h"

#include <map>
#include <string>
#include <vector>

using namespace std;

static constexpr const bool kExtractSnpOnly = false;
static bool extract_snp_only = kExtractSnpOnly;
static constexpr const char* _extract_snp_only = "-v";

static constexpr const bool kExtractPassVariantsOnly = false;
static bool extract_pass_variants_only = kExtractPassVariantsOnly;
static constexpr const char* _extract_pass_variants_only = "-p";

static void
dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] vcf output-dir chr-1 [...chr-n]\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");

    fprintf(stderr, "  %s\n", _extract_snp_only);
    fprintf(stderr, "    Only extract SNPs and ignore other types of variants\n");
    fprintf(stderr, "  %s\n", _extract_pass_variants_only);
    fprintf(stderr, "    Only extract variants having 'PASS' FILETER\n");
}

static bool
parse_arguments(int argc, char* argv[], const char** vcf_path, const char** wrk_dir,
    map<string, int>& chrs)
{
    int i = 2;

    while (i < argc) {
        bool r =(argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (strcmp(argv[i], _extract_snp_only) == 0) {
            extract_snp_only = true;
            i += 1;
            continue;
        }

        if (strcmp(argv[i], _extract_pass_variants_only) == 0) {
            extract_pass_variants_only = true;
            i += 1;
            continue;
        }

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    *vcf_path = argv[i];
    ++i;

    if (i >= argc) return false;
    *wrk_dir = argv[i];
    ++i;

    if (i >= argc) return false;
    for (; i < argc; ++i) chrs[argv[i]] = 0;

    return true;
}

static void
dump_chr_bcf_list(bcf_hdr_t* bcf_hdr, map<string, int>& chrs, const char* wrk_dir, vector<bcf1_t*>& bcf_list)
{
    if (bcf_list.empty()) return;

    string chrname = bcf_seqname(bcf_hdr, bcf_list[0]);
    if (chrs.find(chrname) != chrs.end()) {
        char path[HBN_MAX_PATH_LEN];
        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s", wrk_dir, chrname.c_str());
        create_directory(path);
        if (extract_snp_only) {
            snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/%s.snp.vcf", wrk_dir, chrname.c_str(), chrname.c_str());
        } else {
            snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/%s.vcf", wrk_dir, chrname.c_str(), chrname.c_str());
        }
        fprintf(stderr, "%s\n", path);

        vcfFile* out = vcf_open(path, "w");
        hts_set_threads(out, 8);
        if (vcf_hdr_write(out, bcf_hdr)) HBN_ERR("FAIL at writing VCF header to %s", path);
        for (auto bcf : bcf_list) {
            if (extract_snp_only && (!bcf_is_snp(bcf))) continue;
            if (extract_pass_variants_only && (!bcf_has_filter(bcf_hdr, bcf, "PASS"))) continue;
            if (vcf_write(out, bcf_hdr, bcf)) HBN_ERR("FAIL at writing VCF record to %s", path);
        }
        vcf_close(out);
        chrs[chrname] = 1;
    }

    for (auto bcf : bcf_list) bcf_destroy1(bcf);
    bcf_list.clear();
}

int split_vcf_main(int argc, char* argv[])
{
    const char* vcf_path = nullptr;
    const char* wrk_dir = nullptr;
    map<string, int> chrs;
    if (!parse_arguments(argc, argv, &vcf_path, &wrk_dir, chrs)) {
        dump_usage(argc, argv);
        exit (1);
    }
    create_directory(wrk_dir);

    vcfFile* vcfin = vcf_open(vcf_path, "r");
    hts_set_threads(vcfin, 8);
    bcf_hdr_t* bcfhdr = vcf_hdr_read(vcfin);
    if (!bcfhdr) HBN_ERR("FAIL at reading VCF header from %s", vcf_path);
    vector<bcf1_t*> bcf_list;
    int last_rid = -1;

    while (1) {
        bcf1_t* bcf = bcf_init1();
        int r = vcf_read(vcfin, bcfhdr, bcf);
        if (r == -1) {
            bcf_destroy1(bcf);
            break;
        }
        if (r < 0) HBN_ERR("FAIL at reading VCF record from %s", vcf_path);
        if (bcf->rid != last_rid) {
            dump_chr_bcf_list(bcfhdr, chrs, wrk_dir, bcf_list);
            last_rid = bcf->rid;
        }
        bcf_list.push_back(bcf);
    }
    dump_chr_bcf_list(bcfhdr, chrs, wrk_dir, bcf_list);

    char path[HBN_MAX_PATH_LEN];
    for (auto& chr : chrs) {
        if (chr.second) continue;
        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s", wrk_dir, chr.first.c_str());
        create_directory(path);
        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/%s.vcf", wrk_dir, chr.first.c_str(), chr.first.c_str());
        fprintf(stderr, "%s\n", path);
        vcfFile* out = vcf_open(path, "w");
        if (vcf_hdr_write(out, bcfhdr)) HBN_ERR("FAIL at writing VCF header to %s", path);
        vcf_close(out);
    }

    bcf_hdr_destroy(bcfhdr);
    vcf_close(vcfin);

    return 0;
}