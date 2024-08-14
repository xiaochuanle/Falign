#include "../../htslib/bgzf.h"
#include "../../htslib/faidx.h"
#include "../../htslib/vcf.h"
#include "../../corelib/arg_parse.hpp"

static void
dump_usage(int argc, char* argv[])
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s fasta\n", argv[0], argv[1]);
}

int index_fasta_main(int argc, char* argv[])
{
    if (argc != 3) {
        dump_usage(argc, argv);
        exit (1);
    }
    const char* fasta_path = argv[2];

    int r = fai_build(fasta_path);
    if (r) HBN_ERR("FAIL at indexing %s", fasta_path);
    
    return 0;
}