#include "../../htslib/bgzf.h"
#include "../../htslib/vcf.h"
#include "../../corelib/arg_parse.hpp"

static void
dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s vcf\n", argv[0], argv[1]);
}

int compress_vcf_main(int argc, char* argv[])
{
    if (argc != 3) {
        dump_usage(argc, argv);
        exit (1);
    }
    const char* vcf_path = argv[2];

    BGZF* bzin = bgzf_open(vcf_path, "r");
    int r = bgzf_compression(bzin);
    if (r) {
        fprintf(stderr, "%s is already in compressed format. No compression is needed.\n", vcf_path);
        exit (0);
    }

    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s.gz", vcf_path);

    vcfFile* vcfin = vcf_open(vcf_path, "r");
    hts_set_threads(vcfin, 8);
    bcf_hdr_t* bcfhdr = vcf_hdr_read(vcfin);
    if (!bcfhdr) HBN_ERR("FAIL at reading VCF header from %s", vcf_path);

    vcfFile* vcfout = vcf_open(path, "wb");
    hts_set_threads(vcfout, 8);
    if (vcf_hdr_write(vcfout, bcfhdr)) HBN_ERR("FAIL at writing VCF header to %s", path);
    
    bcf1_t* bcf = bcf_init1();
    while (1) {
        r = vcf_read(vcfin, bcfhdr, bcf);
        if (r == -1) break;
        if (r < 0) HBN_ERR("FAIL reading VCF reacord from %s", vcf_path);
        if (bcf_write(vcfout, bcfhdr, bcf)) HBN_ERR("FAIL at writing VCF record to %s", path);
    }
    bcf_destroy1(bcf);

    bcf_hdr_destroy(bcfhdr);
    vcf_close(vcfout);
    vcf_close(vcfin);

    return 0;
}