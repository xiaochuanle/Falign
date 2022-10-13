#include "../../corelib/hbn_package_version.h"

#include <cstdlib>
#include <cstring>

static void
s_dump_usage(const char* pn)
{
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "  %s <command> [OPTIONS]\n", pn);

    fprintf(stderr, "\n");
    fprintf(stderr, "DESCRIPTION\n");
    fprintf(stderr, " Utility for falign\n");

    fprintf(stderr, "\n");
    fprintf(stderr, "VERSION\n");
    fprintf(stderr, "  %s\n", HBN_PACKAGE_VERSION);

    fprintf(stderr, "\n");
    fprintf(stderr, "COMMANDS\n");
    fprintf(stderr, "  sam2salsa2      Transfer SAM results to pairwise contacts for SALSA2\n");
    fprintf(stderr, "  sam23ddna       Transfer SAM results to pairwise contacts for 3DDNA\n");
    fprintf(stderr, "  paf2salsa2      Transfer PAF results to pairwise contacts for SALSA2\n");
    fprintf(stderr, "  paf23ddna       Transfer PAF results ot pairwise contacts for 3DDNA\n");
    fprintf(stderr, "  sam2frag-sam    Transfer SAM results output by falign to fragment SAM mapping results\n");
}

int make_frag_sam_main(int argc, char* argv[]);
int paf_to_pairwise_contact_bed_main(int argc, char* argv[], const char* target);
int sam_to_pairwise_contact_bed_main(int argc, char* argv[], const char* target);

int main(int argc, char* argv[])
{
    if (argc < 2) {
        s_dump_usage(argv[0]);
        return 1;
    }

    if (strcmp(argv[1], "sam2salsa2") == 0) {
        return sam_to_pairwise_contact_bed_main(argc, argv, "salsa2");
    }
    
    if (strcmp(argv[1], "sam23ddna") == 0) {
        return sam_to_pairwise_contact_bed_main(argc, argv, "3ddna");
    }

    if (strcmp(argv[1], "paf2salsa2") == 0) {
        return paf_to_pairwise_contact_bed_main(argc, argv, "salsa2");
    }

    if (strcmp(argv[1], "paf23ddna") == 0) {
        return paf_to_pairwise_contact_bed_main(argc, argv, "3ddna");
    } 

    if (strcmp(argv[1], "sam2frag-sam") == 0) {
        return make_frag_sam_main(argc, argv);
    }

    s_dump_usage(argv[0]);
    return 1;
}