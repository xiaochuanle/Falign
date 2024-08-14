#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/pdqsort.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "../../htslib/sam.h"
#include "bam-writer.hpp"

#include <algorithm>
#include <vector>

#include <cmath>

using namespace std;

static constexpr const char* _num_threads = "-t";
static constexpr const int kNumThreads = 16;
static int num_threads = kNumThreads;

static constexpr const char* _output_bam_path = "-o";
static constexpr const char* kOutputBamPath = "-";
static const char* output_bam_path = kOutputBamPath;

const char* input_bam_dir = nullptr;

static void
dump_usage(int argc, char* argv[])
{
    string size;
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] input-bam-dir chr-1 [...chr-n]\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");

    fprintf(stderr, "  %s <integer>\n", _num_threads);
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "    Default = '%d'\n", kNumThreads);

    fprintf(stderr, "  %s <integer>\n", _output_bam_path);
    fprintf(stderr, "    Output BAM path\n");
    fprintf(stderr, "    Default = '%s'\n", kOutputBamPath);
}

struct ChrInfo
{
    const char* name;
    int id;

    ChrInfo(const char* _name, int _id): name(_name), id(_id) {}
};

static bool 
parse_arguments(int argc, char* argv[], vector<ChrInfo>& chr_list)
{
    int i = 2;
    while (i < argc) {
        bool r =(argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (parse_int_arg_value(argc, argv, i, _num_threads, num_threads)) continue;
        if (parse_string_arg_value(argc, argv, i, _output_bam_path, output_bam_path)) continue;

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    input_bam_dir = argv[i];
    ++i;

    for (; i < argc; ++i) chr_list.emplace_back(argv[i], -1);

    return !chr_list.empty();
}

static void
dump_parameters(vector<ChrInfo>& chr_list)
{
    string size;
    fprintf(stderr, "====> Parameters:\n");
    fprintf(stderr, "CPU-threads: %d\n", num_threads);
    fprintf(stderr, "Input-BAM-dir: %s\n", input_bam_dir);
    fprintf(stderr, "Output-BAM: %s\n", output_bam_path);
    fprintf(stderr, "Chromosomes:"); for (auto& chr : chr_list) fprintf(stderr, " %s", chr.name); fprintf(stderr, "\n");
    fprintf(stderr, "\n");
}

static inline void
s_make_chr_bam_path(const char* chr_name, char path[])
{
    snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/sorted.snp.bam", input_bam_dir, chr_name);
}

static sam_hdr_t*
s_extract_chr_id_and_length(vector<ChrInfo>& chr_list)
{
    char path[HBN_MAX_PATH_LEN];
    s_make_chr_bam_path(chr_list[0].name, path);
    samFile* in = sam_open(path, "rb");
    if (!in) exit(1);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", path);
    
    for (auto& chr : chr_list) {
        const char* chr_name = chr.name;
        int chr_id = bam_name2id(hdr, chr_name);
        if (chr_id == -1) HBN_ERR("Chromosome %s is not present in BAM file %s", chr_name, path);
        if (chr_id == -2) HBN_ERR("BAM header could not be parsed from file %s", path);
        chr.id = chr_id;
    }

    pdqsort(chr_list.begin(), chr_list.end(), [](const ChrInfo& x, const ChrInfo& y) { return x.id < y.id; });
    return hdr;
}

static void
s_merge_one_bam(const char* chr_name, samFile* out, sam_hdr_t* out_hdr) 
{
    char path[HBN_MAX_PATH_LEN];
    s_make_chr_bam_path(chr_name, path);
    HBN_LOG("Merge %s", path);
    samFile* in = sam_open(path, "rb");
    hts_set_threads(in, num_threads);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", path);

    bam1_t* bam = bam_init1();
    while (1) {
        int r = sam_read1(in, hdr, bam);
        if (r == -1) {
            bam_destroy1(bam);
            break;
        }
        if (r < 0) HBN_ERR("FAIL at reading BAM record from %s", path);

        if (sam_write1(out, out_hdr, bam) == -1) HBN_ERR("FAIL at writing BAM reocrd to %s", output_bam_path);
    }

    sam_hdr_destroy(hdr);
    sam_close(in);
    HBN_LOG("Done");
}

int merge_snp_bam_main(int argc, char* argv[])
{
    vector<ChrInfo> chr_list;
    if (!parse_arguments(argc, argv, chr_list)) {
        dump_usage(argc, argv);
        exit (EXIT_FAILURE);
    }
    dump_parameters(chr_list);

    sam_hdr_t* hdr = s_extract_chr_id_and_length(chr_list);
    samFile* out = sam_open(output_bam_path, "wb");
    if (!out) abort();
    hts_set_threads(out, num_threads);
    int r = sam_hdr_write(out, hdr);
    if (r) HBN_ERR("FAIL at writing BAM header to '%s'", output_bam_path);
    for (auto& chr : chr_list) s_merge_one_bam(chr.name, out, hdr);
    sam_hdr_destroy(hdr);
    sam_close(out);

    return 0;
}
