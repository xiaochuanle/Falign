#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "../../htslib/sam.h"

#include <algorithm>
#include <vector>

using namespace std;

static constexpr const int kNumThreads = 16;

static int num_threads = kNumThreads;

static const char* _num_threads = "-t";

static const char* sorted_bam_path = nullptr;

static void
dump_usage(int argc, char* argv[])
{
    string size;
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] sorted-bam\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  %s <integer>\n", _num_threads);
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "    Default = '%d'\n", kNumThreads);
}

static bool 
parse_arguments(int argc, char* argv[])
{
    int i = 2;
    while (i < argc) {
        bool r =(argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (parse_int_arg_value(argc, argv, i, _num_threads, num_threads)) continue;

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    sorted_bam_path = argv[i];
    ++i;
    return true;
}

static void
dump_parameters()
{
    string size;
    fprintf(stderr, "====> Parameters:\n");
    fprintf(stderr, "CPU-threads: %d\n", num_threads);
    fprintf(stderr, "Sorted-BAM: %s\n", sorted_bam_path);
    fprintf(stderr, "\n");
}

int index_bam_main(int argc, char* argv[])
{
    if (!parse_arguments(argc, argv)) {
        dump_usage(argc, argv);
        exit (EXIT_FAILURE);
    }
    dump_parameters();

    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s.bai", sorted_bam_path);

    int r = sam_index_build3(sorted_bam_path, path, 0, num_threads);
    if (r) HBN_ERR("FAIL at creating index for %s", sorted_bam_path);

    return 0;
}