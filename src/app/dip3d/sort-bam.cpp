#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/pdqsort.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "../../htslib/sam.h"
#include "bam-writer.hpp"

#include <algorithm>
#include <vector>

using namespace std;

constexpr const size_t kInputBamChunkSize = 10000000000;
constexpr const size_t kOutputBamChunkSize = 1000000000;
constexpr const int kNumThreads = 8;
constexpr const char* kSortedBamPath = "-";

static size_t input_bam_chunk_size = kInputBamChunkSize;
static size_t output_bam_chunk_size = kOutputBamChunkSize;
static int num_threads = kNumThreads;
static const char* sorted_bam_path = kSortedBamPath;

static const char* _input_bam_chunk_size = "-w";
static const char* _output_bam_chunk_size = "-u";
static const char* _num_threads = "-t";
static const char* _sorted_bam_path = "-o";

static void
dump_usage(int argc, char* argv[])
{
    string size;
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] bam-1 [...bam-n]\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  %s <DataSize>\n", _input_bam_chunk_size);
    fprintf(stderr, "    Each time process this number of bps in RAM\n");
    size = NStr::UInt8ToString_DataSize(kInputBamChunkSize);
    fprintf(stderr, "    Default = '%s'\n", size.c_str());
    fprintf(stderr, "  %s <DataSize>\n", _output_bam_chunk_size);
    fprintf(stderr, "    Output buffer size\n");
    size = NStr::UInt8ToString_DataSize(kOutputBamChunkSize);
    fprintf(stderr, "    Default = '%s'\n", size.c_str());
    fprintf(stderr, "  %s <integer>\n", _num_threads);
    fprintf(stderr, "    Number of CPU threads used\n");
    fprintf(stderr, "    Default = '%d'\n", kNumThreads);
    fprintf(stderr, "  %s <path>\n", _sorted_bam_path);
    fprintf(stderr, "    Path to save sorted BAM file\n");
    fprintf(stderr, "    Default = '%s'\n", kSortedBamPath);
}

static bool 
parse_arguments(int argc, char* argv[], vector<string>& input_bam_list)
{
    int i = 2;
    while (i < argc) {
        bool r =(argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (parse_data_size_arg_value(argc, argv, i, _input_bam_chunk_size, input_bam_chunk_size)) continue;
        if (parse_data_size_arg_value(argc, argv, i, _output_bam_chunk_size, output_bam_chunk_size)) continue;
        if (parse_int_arg_value(argc, argv, i, _num_threads, num_threads)) continue;
        if (parse_string_arg_value(argc, argv, i, _sorted_bam_path, sorted_bam_path)) continue;

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    for (; i < argc; ++i) input_bam_list.push_back(argv[i]);

    return !input_bam_list.empty();
}

static void
dump_parameters(vector<string>& input_bam_list)
{
    string size;
    fprintf(stderr, "====> Parameters:\n");
    size = NStr::UInt8ToString_DataSize(input_bam_chunk_size);
    fprintf(stderr, "Input-bam-chunk-size: %s\n", size.c_str());
    size = NStr::UInt8ToString_DataSize(output_bam_chunk_size);
    fprintf(stderr, "Output-bam-chunk-size: %s\n", size.c_str());
    fprintf(stderr, "CPU-threads: %d\n", num_threads);
    fprintf(stderr, "Sorted-BAM: %s\n", sorted_bam_path);
    fprintf(stderr, "Input-bams:"); for (auto& s : input_bam_list) fprintf(stderr, " %s", s.c_str()); fprintf(stderr, "\n");
    fprintf(stderr, "\n");
}

int sort_bam_main(int argc, char* argv[])
{
    vector<string> input_bam_list;
    if (!parse_arguments(argc, argv, input_bam_list)) {
        dump_usage(argc, argv);
        exit (EXIT_FAILURE);
    }
    dump_parameters(input_bam_list);

    vector<BamInfo> chunk_bam_list;
    char path[HBN_MAX_PATH_LEN];
    int chunk_id = 0;
    size_t num_frags = 0, num_bases = 0;
    size_t chunk_frags = 0, chunk_bases = 0;
    string size;
    sam_hdr_t* outhdr = nullptr;
    for (auto& s : input_bam_list) {
        HBN_LOG("Process %s", s.c_str());
        samFile* in = sam_open(s.c_str(), "rb");
        hts_set_threads(in, num_threads);
        sam_hdr_t* hdrin = sam_hdr_read(in);
        if (!outhdr) outhdr = sam_hdr_dup(hdrin);
        while (1) {
            bam1_t* bam = bam_init1();
            int r = sam_read1(in, hdrin, bam);
            if (r == -1) {
                bam_destroy1(bam);
                break;
            }
            if (r < 0) HBN_ERR("FAIL at reading BAM record from %s", s.c_str());
            chunk_bam_list.push_back(BamInfo(bam->core.tid, bam->core.pos, bam));
            chunk_bases += bam->core.l_qseq;
            ++chunk_frags;
            if (chunk_bases >= kInputBamChunkSize) {
                make_bam_chunk_path(sorted_bam_path, chunk_id, path);
                size = NStr::UInt8ToString_DataSize(chunk_bases);
                HBN_LOG("Save %zu BAM records (%s) to %s", chunk_frags, size.c_str(), path);
                samFile* out = sam_open(path, "wb");
	            hts_set_threads(out, num_threads);
                r = sam_hdr_write(out, outhdr);
                if (r) HBN_ERR("FAIL at writing BAM header");
                pdqsort(chunk_bam_list.begin(), chunk_bam_list.end());
                for (auto& bi : chunk_bam_list) {
                    r = sam_write1(out, outhdr, bi.bam);
                    if (r < 0) HBN_ERR("FAIL at writing BAM record");
                }
                for (auto& bi : chunk_bam_list) {
                    bam_destroy1(bi.bam);
                }
                sam_close(out);

                num_frags += chunk_frags;
                num_bases += chunk_bases;
                chunk_frags = 0;
                chunk_bases = 0;
                chunk_bam_list.clear();
                ++chunk_id;
            }           
        }
        sam_hdr_destroy(hdrin);
        sam_close(in);
    }
    if (chunk_frags) {
        make_bam_chunk_path(sorted_bam_path, chunk_id, path);
        size = NStr::UInt8ToString_DataSize(chunk_bases);
        HBN_LOG("Save %zu BAM records (%s) to %s", chunk_frags, size.c_str(), path);
        samFile* out = sam_open(path, "wb");
	    hts_set_threads(out, num_threads);
        int r = sam_hdr_write(out, outhdr);
        if (r) HBN_ERR("FAIL at writing BAM header");
        pdqsort(chunk_bam_list.begin(), chunk_bam_list.end());
        for (auto& bi : chunk_bam_list) {
            r = sam_write1(out, outhdr, bi.bam);
            if (r < 0) HBN_ERR("FAIL at writing BAM record");
        }
        for (auto& bi : chunk_bam_list) {
            bam_destroy1(bi.bam);
        }
        sam_close(out);

        num_frags += chunk_frags;
        num_bases += chunk_bases;
        chunk_frags = 0;
        chunk_bases = 0;
        chunk_bam_list.clear();
        ++chunk_id;        
    }
    const int num_chunks = chunk_id;
    size = NStr::UInt8ToString_DataSize(num_bases);
    HBN_LOG("Split %zu BAM records (%s) into %d chunks", num_frags, size.c_str(), num_chunks);

    HBN_LOG("Merging %d BAM chunks into %s", num_chunks, sorted_bam_path);
    BamChunkReader** chunk_list = new BamChunkReader*[num_chunks];
    for (int i = 0; i < num_chunks; ++i) {
        make_bam_chunk_path(sorted_bam_path, i, path);
        chunk_list[i] = new BamChunkReader(path, num_threads);
    }
    BamWriter* out = new BamWriter(sorted_bam_path, outhdr, output_bam_chunk_size, num_threads);

    while (1) {
        bam1_t* bam = select_next_bam(chunk_list, num_chunks);
        if (!bam) break;
        out->save_one_bam(bam);
    }

    for (int i = 0; i < num_chunks; ++i) delete chunk_list[i];
    delete[] chunk_list;
    sam_hdr_destroy(outhdr);
    delete out;

    for (int i = 0; i < num_chunks; ++i) {
        make_bam_chunk_path(sorted_bam_path, i, path);
        fprintf(stderr, "Remove file %s\n", path);
        remove(path);
    }

    HBN_LOG("Create index for %s", sorted_bam_path);
    snprintf(path, HBN_MAX_PATH_LEN, "%s.bai", sorted_bam_path);
    int r = sam_index_build3(sorted_bam_path, path, 0, num_threads);
    if (r) HBN_ERR("FAIL at creating index for %s", sorted_bam_path);
    HBN_LOG("Done");

    return 0;
}