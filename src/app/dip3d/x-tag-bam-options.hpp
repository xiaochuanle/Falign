#ifndef __X_TAG_BAM_OPTIONS_HPP
#define __X_TAG_BAM_OPTIONS_HPP

#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"

#include <cstring>

struct BamTagOptions
{
    int snp_tag_min_mapQ { 5 };
    int snp_tag_match_base { 3 };
    double snp_tag_identity { 85.0 };
    
    int impute_adj_frag_dist { 29500000 };
    int impute_non_adj_frag_dist { 16500000 };

    int num_threads { 1 };

    const char* wrk_dir { nullptr };
    const char* reference_path { nullptr };
    const char* vcf_path { nullptr };
    const char* bam_path { nullptr };

    void dump() {
        fprintf(stderr, "\n");
        fprintf(stderr, "====> Parameters:\n");
        fprintf(stderr, "SNP-tag-min-mapQ: %d\n", snp_tag_min_mapQ);
        fprintf(stderr, "SNP-tag-match-bases: %d\n", snp_tag_match_base);
        fprintf(stderr, "SNP-tag-min-identity: %g\n", snp_tag_identity);
        fprintf(stderr, "Impute-adj-dist: %d\n", impute_adj_frag_dist);
        fprintf(stderr, "Impute-non-adj-dist: %d\n", impute_non_adj_frag_dist);
        fprintf(stderr, "CPU-threads: %d\n", num_threads);
        fprintf(stderr, "Wrk-dir: %s\n", wrk_dir);
        fprintf(stderr, "reference: %s\n", reference_path);
        fprintf(stderr, "VCF: %s\n", vcf_path);
        fprintf(stderr, "BAM: %s\n", bam_path);
        fprintf(stderr, "\n");
    }

    void dump_usage(int argc, char* argv[]) {
        fprintf(stderr, "\n");
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s [OPTIONS] wrk-dir reference vcf input-bam\n", argv[0], argv[1]);

        fprintf(stderr, "\n");
        fprintf(stderr, "OPTIONAL ARGUMENTS\n");
        
        fprintf(stderr, "  *** SNP-tagging options\n");
        fprintf(stderr, "  -q <Integer>\n");
        fprintf(stderr, "    Minimum mapping quality\n");
        fprintf(stderr, "    Default = '%d'\n", snp_tag_min_mapQ);
        fprintf(stderr, "  -m <Integer>\n");
        fprintf(stderr, "    Number of match bases surrounding overlapped SNP\n");
        fprintf(stderr, "    Default = '%d'\n", snp_tag_match_base);
        fprintf(stderr, "  -i <Percentage>\n");
        fprintf(stderr, "    Minimum alignment identity\n");
        fprintf(stderr, "    Default '%g'\n", snp_tag_identity);

        fprintf(stderr, "  *** Imputation options\n");
        fprintf(stderr, "  -a <Integer>\n");
        fprintf(stderr, "    Maximum genomic distance between adjacent fragment pair\n");
        fprintf(stderr, "    Default ='%d'\n", impute_adj_frag_dist);
        fprintf(stderr, "  -d <Integer>\n");
        fprintf(stderr, "    Maximum genomic distance between non-adjacent fragment pair\n");
        fprintf(stderr, "    Default = '%d'\n", impute_non_adj_frag_dist);

        fprintf(stderr, "  *** General running options\n");
        fprintf(stderr, "  -t <Integer>\n");
        fprintf(stderr, "    Number of CPU threads used\n");
        fprintf(stderr, "    Default = '%d'\n", num_threads);
    }

    bool parse_arguments(int argc, char* argv[]) {
        int i = 2;
        while (i < argc) {
            bool r = (argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
            if (r) break;

            if (parse_int_arg_value(argc, argv, i, "-q", snp_tag_min_mapQ)) continue;
            if (parse_int_arg_value(argc, argv, i, "-m", snp_tag_match_base)) continue;
            if (parse_real_arg_value(argc, argv, i, "-i", snp_tag_identity)) continue;

            if (parse_int_arg_value(argc, argv, i, "-a", impute_adj_frag_dist)) continue;
            if (parse_int_arg_value(argc, argv, i, "-d", impute_non_adj_frag_dist)) continue;

            if (parse_int_arg_value(argc, argv, i, "-t", num_threads)) continue;

            fprintf(stderr, "ERROR: Unrecognised option %s\n", argv[i]);
            return false;
        }

        if (i >= argc) return false;
        wrk_dir = argv[i];
        ++i;

        if (i >= argc) return false;
        reference_path = argv[i];
        ++i;

        if (i >= argc) return false;
        vcf_path = argv[i];
        ++i;

        if (i >= argc) return false;
        bam_path = argv[i];
        ++i;

        return true;
    }
};

#endif // __X_TAG_BAM_OPTIONS_HPP