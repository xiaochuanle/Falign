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

static int enzyme_dist = 50;

static constexpr const int kMinFragSize = 100;
static constexpr const int kMinMapQ = 5;
static constexpr const double kMinIdentity = 90.0;
constexpr const size_t kInputBamChunkSize = 10000000000;
static constexpr const size_t kOutputBamChunkSize = 1000000000;
static constexpr const int kNumThreads = 8;

static const char* _min_frag_size = "-l";
static const char* _min_mapQ = "-q";
static const char* _min_identity = "-i";
static const char* _input_bam_chunk_size = "-w";
static const char* _output_bam_chunk_size = "-u";
static const char* _num_threads = "-t";

static int min_frag_size = kMinFragSize;
static int min_mapQ = kMinMapQ;
static double min_identity = kMinIdentity;
static size_t input_bam_chunk_size = kInputBamChunkSize;
static size_t output_bam_chunk_size = kOutputBamChunkSize;
static int num_threads = kNumThreads;

static const char* input_bam_path = nullptr;
static size_t snp_data_size = 0;
static const char* snp_bam_path = nullptr;

static void
dump_usage(int argc, char* argv[])
{
    string size;
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] input-bam snp-data-size snp-bam\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");

    fprintf(stderr, "  %s <integer>\n", _min_frag_size);
    fprintf(stderr, "    Minimum fragment length (in bp)\n");
    fprintf(stderr, "    Default = '%d'\n", kMinFragSize);

    fprintf(stderr, "  %s <integer>\n", _min_mapQ);
    fprintf(stderr, "    Minimum mapping quality\n");
    fprintf(stderr, "    Default = '%d'\n", kMinMapQ);

    fprintf(stderr, "  %s <percentage>\n", _min_identity);
    fprintf(stderr, "    Minimum alignment identity\n");
    fprintf(stderr, "    Default = '%g'\n", kMinIdentity);

    fprintf(stderr, "  %s <DataSize>\n", _input_bam_chunk_size);
    fprintf(stderr, "    Each time process this number of bps in RAM\n");
    size = NStr::UInt8ToString_DataSize(kInputBamChunkSize);
    fprintf(stderr, "    Default = '%s'\n", size.c_str());

    fprintf(stderr, "  %s <DataSize>\n", _output_bam_chunk_size);
    fprintf(stderr, "    Output BAM buffer size\n");
    size = NStr::UInt8ToString_DataSize(kOutputBamChunkSize);
    fprintf(stderr, "    Default = '%s'\n", size.c_str());

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

        if (parse_int_arg_value(argc, argv, i, _min_frag_size, min_frag_size)) continue;
        if (parse_int_arg_value(argc, argv, i, _min_mapQ, min_mapQ)) continue;
        if (parse_real_arg_value(argc, argv, i, _min_identity, min_identity)) continue;
        if (parse_data_size_arg_value(argc, argv, i, _input_bam_chunk_size, input_bam_chunk_size)) continue;
        if (parse_data_size_arg_value(argc, argv, i, _output_bam_chunk_size, output_bam_chunk_size)) continue;
        if (parse_int_arg_value(argc, argv, i, "-t", num_threads)) continue;

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    input_bam_path = argv[i];
    ++i;

    if (i >= argc) return false;
    snp_data_size = NStr::StringToUInt8_DataSize(argv[i]);
    ++i;

    if (i >= argc) return false;
    snp_bam_path = argv[i];
    ++i;

    return true;
}

static void
dump_parameters()
{
    string size;
    fprintf(stderr, "====> Parameters:\n");
    fprintf(stderr, "min-frag-size: %dbp\n", min_frag_size);
    fprintf(stderr, "min-mapQ: %d\n", min_mapQ);
    fprintf(stderr, "min-identity: %g\n", min_identity);
    size = NStr::UInt8ToString_DataSize(input_bam_chunk_size);
    fprintf(stderr, "Input-bam-chunk-size: %s\n", size.c_str());
    size = NStr::UInt8ToString_DataSize(output_bam_chunk_size);
    fprintf(stderr, "Output-bam-chunk-size: %s\n", size.c_str());
    fprintf(stderr, "CPU-threads: %d\n", num_threads);
    fprintf(stderr, "Input-BAM: %s\n", input_bam_path);
    size = NStr::UInt8ToString_DataSize(snp_data_size);
    fprintf(stderr, "SNP-data-size: %s\n", size.c_str());
    fprintf(stderr, "SNP-BAM: %s\n", snp_bam_path);
    fprintf(stderr, "\n");
}

static int
s_select_and_split_snp_bams(samFile* in, sam_hdr_t* hdrin, const char* snp_bam_path)
{
    char path[HBN_MAX_PATH_LEN];
    int chunk_id = 0;
    char pitagname[2] = { 'p', 'i' };
    vector<BamInfo> bam_list;
    size_t total_frags = 0, total_bases = 0;
    string size;

    while (1) {
        bool eof = false;
        size_t num_frags = 0, num_bases = 0;
        bam_list.clear();
        while (1) {
            bam1_t* bam = bam_init1();
            int r = sam_read1(in, hdrin, bam);
            if (r == -1) {
                eof = true;
                bam_destroy1(bam);
                break;
            }
            if (r < 0) HBN_ERR("FAIL at reading BAM record");
            int frag_size = bam->core.l_qseq;
            int mapQ = bam->core.qual;
	    pitagname[0] = 'p';
	    pitagname[1] = 'i';
            uint8_t* pitag = bam_aux_get(bam, pitagname);
            if (!pitag) HBN_ERR("FAIL at retrieving alignment identity from %c%c BAM tag", pitagname[0], pitagname[1]);
            double pi = bam_aux2f(pitag);         


            r = (frag_size >= min_frag_size) && (mapQ >= min_mapQ) && (pi >= min_identity);
            if (!r) {
                bam_destroy1(bam);
		        continue;    
	        }

	    pitagname[0] = 'q';
	    pitagname[1] = 's';
	    uint8_t* qstag = bam_aux_get(bam, pitagname);
	    int qb = bam_aux2i(qstag);

	    pitagname[0] = 'q';
	    pitagname[1] = 'S';
	    uint8_t* qStag = bam_aux_get(bam, pitagname);
	    int qS = bam_aux2i(qStag);

	    pitagname[0] = 'q';
	    pitagname[1] = 'e';
	    uint8_t* qetag = bam_aux_get(bam, pitagname);
	    int qe = bam_aux2i(qetag);

	    pitagname[0] = 'q';
	    pitagname[1] = 'E';
	    uint8_t* qEtag = bam_aux_get(bam, pitagname);
	    int qE = bam_aux2i(qEtag);

	    pitagname[0] = 's';
	    pitagname[1] = 's';
	    uint8_t* sstag = bam_aux_get(bam, pitagname);
	    int sb = bam_aux2i(sstag);

	    pitagname[0] = 'v';
	    pitagname[1] = 'S';
	    uint8_t* sStag = bam_aux_get(bam, pitagname);
	    int sS = bam_aux2i(sStag);

	    pitagname[0] = 's';
	    pitagname[1] = 'e';
	    uint8_t* setag = bam_aux_get(bam, pitagname);
	    int se = bam_aux2i(setag);

	    pitagname[0] = 'v';
	    pitagname[1] = 'E';
	    uint8_t* sEtag = bam_aux_get(bam, pitagname);
	    int sE = bam_aux2i(sEtag);

	    //r = (abs(qb - qS) <= enzyme_dist) || (abs(qe - qE) <= enzyme_dist) || (abs(sb - sS) <= enzyme_dist) || (abs(se - sE) <= enzyme_dist);
	    //if (!r) { bam_destroy1(bam); continue; }

            bam_list.push_back(BamInfo(bam->core.tid, bam->core.pos, bam));
            ++num_frags;   
            ++total_frags;
            num_bases += frag_size;   
            total_bases += frag_size;

            if (num_bases >= input_bam_chunk_size || total_bases >= snp_data_size) break;     
        }

        if (!bam_list.empty()) {
            make_bam_chunk_path(snp_bam_path, chunk_id, path);
            ++chunk_id;
            size = NStr::UInt8ToString_DataSize(num_bases);
            HBN_LOG("Save %zu BAM records (%s) into %s", num_frags, size.c_str(), path);
            samFile* out = sam_open(path, "wb");
	        hts_set_threads(out, num_threads);
            int r = sam_hdr_write(out, hdrin);
            if (r) HBN_ERR("FAIL at writing BAM header");
            pdqsort(bam_list.begin(), bam_list.end());
            for (auto& bi : bam_list) {
                r = sam_write1(out, hdrin, bi.bam);
                if (r < 0) HBN_ERR("FAIL at writing BAM record");
            }
            for (auto& bi : bam_list) bam_destroy1(bi.bam);
            bam_list.clear();
            sam_close(out);
            HBN_LOG("Done");
        }

        if (eof || num_bases >= snp_data_size) break;
    }

    size = NStr::UInt8ToString_DataSize(total_bases);
    HBN_LOG("Extract %zu BAM recordds (%s) and split into %d chunks", total_frags, size.c_str(), chunk_id);
    return chunk_id;
}

int select_snp_bam_main(int argc, char* argv[])
{
    if (!parse_arguments(argc, argv)) {
        dump_usage(argc, argv);
        exit (1);
    }
    dump_parameters();

    samFile* in = sam_open(input_bam_path, "rb");
    hts_set_threads(in, num_threads);
    sam_hdr_t* hdrin = sam_hdr_read(in);
    const int num_chunks = s_select_and_split_snp_bams(in, hdrin, snp_bam_path);

    HBN_LOG("Merging %d BAM chunks into %s", num_chunks, snp_bam_path);
    BamChunkReader** chunk_list = new BamChunkReader*[num_chunks];
    char path[HBN_MAX_PATH_LEN];
    for (int i = 0; i < num_chunks; ++i) {
        make_bam_chunk_path(snp_bam_path, i, path);
        chunk_list[i] = new BamChunkReader(path, num_threads);
    }
    BamWriter* out = new BamWriter(snp_bam_path, hdrin, output_bam_chunk_size, num_threads);

    while (1) {
        bam1_t* bam = select_next_bam(chunk_list, num_chunks);
        if (!bam) break;
        out->save_one_bam(bam);
    }

    for (int i = 0; i < num_chunks; ++i) delete chunk_list[i];
    delete[] chunk_list;
    sam_hdr_destroy(hdrin);
    sam_close(in);
    delete out;

    for (int i = 0; i < num_chunks; ++i) {
        make_bam_chunk_path(snp_bam_path, i, path);
        fprintf(stderr, "Remove file %s\n", path);
        remove(path);
    }

    HBN_LOG("Create index for %s", snp_bam_path);
    snprintf(path, HBN_MAX_PATH_LEN, "%s.bai", snp_bam_path);
    int r = sam_index_build3(snp_bam_path, path, 0, num_threads);
    if (r) HBN_ERR("FAIL at creating index for %s", snp_bam_path);
    HBN_LOG("Done");

    return 0;
}
