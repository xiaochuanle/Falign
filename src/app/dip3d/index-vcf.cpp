#include "../../htslib/bgzf.h"
#include "../../htslib/tbx.h"
#include "../../htslib/vcf.h"
#include "../../corelib/arg_parse.hpp"

static void
dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s vcf.gz\n", argv[0], argv[1]);
}

int index_vcf_main(int argc, char* argv[])
{
    if (argc != 3) {
        dump_usage(argc, argv);
        exit (1);
    }
    const char* gz_vcf_path = argv[2];
    char idx_path[HBN_MAX_PATH_LEN];
    snprintf(idx_path, HBN_MAX_PATH_LEN, "%s.tbi", gz_vcf_path);

    int r = tbx_index_build3(gz_vcf_path, idx_path, 0, 8, &tbx_conf_vcf);
    if (r) HBN_ERR("FAIL at indexing %s", gz_vcf_path);
    
    return 0;
}