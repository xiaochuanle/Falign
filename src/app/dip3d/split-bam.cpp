#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/pdqsort.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "../../htslib/sam.h"
#include "bam-writer.hpp"

#include <algorithm>
#include <map>
#include <string>
#include <vector>

using namespace std;

constexpr const size_t kInputBamChunkSize = 10000000000;
static constexpr const int kNumThreads = 8;

static const char* _input_bam_chunk_size = "-w";
static const char* _num_threads = "-t";

static size_t input_bam_chunk_size = kInputBamChunkSize;
static int num_threads = kNumThreads;

static const char* input_bam_path = nullptr;
static const char* output_dir = nullptr;

static void
dump_usage(int argc, char* argv[])
{
    string size;
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] input-bam ouput-dir [chr1 ...chr-n]\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");

    fprintf(stderr, "  %s <DataSize>\n", _input_bam_chunk_size);
    fprintf(stderr, "    Each time process this number of bps in RAM\n");
    size = NStr::UInt8ToString_DataSize(kInputBamChunkSize);
    fprintf(stderr, "    Default = '%s'\n", size.c_str());

    fprintf(stderr, "  %s <integer>\n", _num_threads);
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "    Default = '%d'\n", kNumThreads);
}

static bool 
parse_arguments(int argc, char* argv[], vector<string>& chr_list)
{
    int i = 2;
    while (i < argc) {
        bool r =(argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (parse_data_size_arg_value(argc, argv, i, _input_bam_chunk_size, input_bam_chunk_size)) continue;
        if (parse_int_arg_value(argc, argv, i, "-t", num_threads)) continue;

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    input_bam_path = argv[i];
    ++i;

    if (i >= argc) return false;
    output_dir = argv[i];
    ++i;

    for (; i < argc; ++i) chr_list.push_back(argv[i]);

    return true;
}

static void
dump_parameters(vector<string>& chr_list)
{
    string size;
    fprintf(stderr, "====> Parameters:\n");
    size = NStr::UInt8ToString_DataSize(input_bam_chunk_size);
    fprintf(stderr, "Input-bam-chunk-size: %s\n", size.c_str());
    fprintf(stderr, "CPU-threads: %d\n", num_threads);
    fprintf(stderr, "Input-BAM: %s\n", input_bam_path);
    fprintf(stderr, "Output-dir: %s\n", output_dir);
    fprintf(stderr, "Chr:"); for (auto &chr : chr_list) fprintf(stderr, " %s", chr.c_str()); fprintf(stderr, "\n");
    fprintf(stderr, "\n");
}

static void
create_chr_bam_out(vector<string>& chr_list, sam_hdr_t* hdrin, map<int, samFile*>& outs)
{
    create_directory(output_dir);
    char path[HBN_MAX_PATH_LEN];
    for (auto& chr : chr_list) {
        int chr_id = sam_hdr_name2tid(hdrin, chr.c_str());
        if (chr_id == -1) {
            continue;
        } else if (chr_id == -2) {
            HBN_ERR("FAIL at parsing BAM header from %s", input_bam_path);
        }
        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s", output_dir, chr.c_str());
        create_directory(path);
        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/%s.bam", output_dir, chr.c_str(), chr.c_str());
        fprintf(stderr, "%s\n", path);
        samFile* out = sam_open(path, "wb");
        hts_set_threads(out, num_threads);
        int r = sam_hdr_write(out, hdrin);
        if (r) HBN_ERR("FAIL at writing BAM header to '%s'", path);
        outs[chr_id] = out;
    }
}

static void 
s_dump_bam_list(vector<bam1_t*>& bam_list, sam_hdr_t* hdrin, map<int, samFile*>& outs)
{
    if (bam_list.empty()) return;
    pdqsort(bam_list.begin(), bam_list.end(), [](bam1_t* x, bam1_t* y) { return x->core.tid < y->core.tid; });
    bam1_t** ba = bam_list.data();
    size_t bc = bam_list.size();
    size_t i = 0;
    while (i < bc) {
        int tid = ba[i]->core.tid;
        size_t j = i + 1;
        while (j < bc && ba[j]->core.tid == tid) ++j;

        auto pos = outs.find(tid);
        if (pos != outs.end()) {
            samFile* out = pos->second;
            for (size_t k = i; k < j; ++k) {
                if (sam_write1(out, hdrin, ba[k]) == -1) HBN_ERR("FAIL at writing BAM reocrd");
            }
        }

        i = j;
    }
}

int split_bam_main(int argc, char* argv[])
{
    vector<string> chr_list;
    if (!parse_arguments(argc, argv, chr_list)) {
        dump_usage(argc, argv);
        exit (EXIT_FAILURE);
    }

    samFile* in = sam_open(input_bam_path, "rb");
    hts_set_threads(in, num_threads);
    sam_hdr_t* hdrin = sam_hdr_read(in);

    if (chr_list.empty()) {
        int num_chr = sam_hdr_nref(hdrin);
        for (int i = 0; i < num_chr; ++i) chr_list.push_back(sam_hdr_tid2name(hdrin, i));
    }
    dump_parameters(chr_list);

    map<int, samFile*> outs;
    create_chr_bam_out(chr_list, hdrin, outs);
    vector<bam1_t*> bam_list;

    size_t total_frags = 0, total_bases = 0;
    string size;
    while (1) {
        bool eof = false;
        size_t num_frags = 0, num_bases = 0;
        while (1) {
            bam1_t* bam = bam_init1();
            int r = sam_read1(in, hdrin, bam);
            if (r == -1) {
                eof = true;
                bam_destroy1(bam);
                break;
            }
            if (r < 0) HBN_ERR("FAIL at reading BAM record");

            bam_list.push_back(bam);
            ++num_frags;   
            num_bases += bam->core.l_qseq;

            if (num_bases >= input_bam_chunk_size) break;     
        }
        if (bam_list.empty()) break;

        size = NStr::UInt8ToString_DataSize(num_bases);
        HBN_LOG("Save %zu BAM records (%s)", num_frags, size.c_str());
        s_dump_bam_list(bam_list, hdrin, outs);
        for (auto bam : bam_list) bam_destroy1(bam);
        bam_list.clear();

        total_frags += num_frags;
        total_bases += num_bases;
        size = NStr::UInt8ToString_DataSize(total_bases);
        HBN_LOG("%zu BAM records (%s) dumpped", total_frags, size.c_str());

        if (eof) break;
    }
    for (auto xout : outs) sam_close(xout.second);

    char path[HBN_MAX_PATH_LEN];
    for (auto& chr : chr_list) {
        int chr_id = sam_hdr_name2tid(hdrin, chr.c_str());
        if (chr_id == -2) {
            HBN_ERR("FAIL at parsing BAM header from %s", input_bam_path);
        }
        if (chr_id != -1) continue;

        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s", output_dir, chr.c_str());
        create_directory(path);
        snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/%s.vcf", output_dir, chr.c_str(), chr.c_str());
        fprintf(stderr, "%s\n", path);
        samFile* out = sam_open(path, "wb");
        hts_set_threads(out, num_threads);
        int r = sam_hdr_write(out, hdrin);
        if (r) HBN_ERR("FAIL at writing BAM header to '%s'", path);
        sam_close(out);
    }

    return 0;
}
