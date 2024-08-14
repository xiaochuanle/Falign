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

static constexpr const int kMinFragSize = 100;
static constexpr const int kMinMapQ = 5;
static constexpr const double kMinIdentity = 90.0;
constexpr const size_t kInputBamChunkSize = 10000000000;
static constexpr const size_t kOutputBamChunkSize = 1000000000;
static constexpr const int kNumThreads = 16;

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

static const char* chr_bam_dir = nullptr;
static int snp_frag_coverage = 0;
static const char* snp_bam_path = nullptr;

///////////////////

static u8* build_cnt_table()
{
    const u8 kMaxCov = 200;
    u8* cnt_table = new u8[256];
    for (u8 i = 0; i <= kMaxCov; ++i) cnt_table[i] = i + 1;
    cnt_table[kMaxCov] = kMaxCov;

    return cnt_table;
}

static bool
s_bam_cov_is_full(bam1_t* bam, const u8* cnt_table, const int max_cov, u8* cov_stats)
{
    uint8_t* sstag = bam_aux_get(bam, "ss");
    hbn_assert(sstag);
    int ss = bam_aux2i(sstag);
    uint8_t* setag = bam_aux_get(bam, "se");
    hbn_assert(setag);
    int se = bam_aux2i(setag);

    int full_sites = 0;
    for (int i = ss; i < se; ++i) if (cov_stats[i] >= max_cov) ++full_sites;
    if (se - ss - full_sites >= 20) {
        for (int i = ss; i < se; ++i) cov_stats[i] = cnt_table[ cov_stats[i] ];
        return false;
    }

    return true;
}

struct SnpChrInfo
{
    const char* name;
    int id;
    int length;

    SnpChrInfo(const char* _name, int _id, int _length): name(_name), id(_id), length(_length) {}
};

struct SnpBamInfo
{
    int chr_id;
    int chr_pos;
    bam1_t* bam;

    SnpBamInfo(int _id, int _pos, bam1_t* _bam): chr_id(_id), chr_pos(_pos), bam(_bam) {}

    bool operator<(const SnpBamInfo& rhs) const {
        return (this->chr_id < rhs.chr_id) || (this->chr_id == rhs.chr_id && this->chr_pos < rhs.chr_pos);
    }
};

static inline void
s_make_chr_bam_path(const char* chr_name, char path[])
{
    snprintf(path, HBN_MAX_PATH_LEN, "%s/%s/%s.bam", chr_bam_dir, chr_name, chr_name);
}

static void
dump_usage(int argc, char* argv[])
{
    string size;
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] chr-bam-dir coverage snp-bam chr-1 [...chr-n]\n", argv[0], argv[1]);

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
parse_arguments(int argc, char* argv[], vector<SnpChrInfo>& chr_list)
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
    chr_bam_dir = argv[i];
    ++i;

    if (i >= argc) return false;
    snp_frag_coverage = NStr::StringToInt(argv[i]);
    ++i;

    if (i >= argc) return false;
    snp_bam_path = argv[i];
    ++i;

    for (; i < argc; ++i) chr_list.emplace_back(argv[i], -1, -1);

    return !chr_list.empty();
}

static void
dump_parameters(vector<SnpChrInfo>& chr_list)
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
    fprintf(stderr, "Chr-BAM-dir: %s\n", chr_bam_dir);
    fprintf(stderr, "Coverage: %d\n", snp_frag_coverage);
    fprintf(stderr, "SNP-BAM: %s\n", snp_bam_path);
    fprintf(stderr, "Chromosomes:"); for (auto& chr : chr_list) fprintf(stderr, " %s", chr.name); fprintf(stderr, "\n");
    fprintf(stderr, "\n");
}

static sam_hdr_t*
s_extract_chr_id_and_length(vector<SnpChrInfo>& chr_list)
{
    char path[HBN_MAX_PATH_LEN];
    s_make_chr_bam_path(chr_list[0].name, path);
    samFile* in = sam_open(path, "rb");
    if (!in) HBN_ERR("Could not open file %s for reading", path);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", path);

    for (auto& chr : chr_list) {
        int chr_id = sam_hdr_name2tid(hdr, chr.name);
        if (chr_id == -1) HBN_ERR("Chromosome %s is not present in BAM file %s", chr.name, path);
        if (chr_id == -2) HBN_ERR("BAM header could not be parsed from file %s", path);
        int chr_length = sam_hdr_tid2len(hdr, chr_id);
        if (chr_length == 0) HBN_ERR("Could not extract chromosome length from BAM header of %s", path);
        string size = NStr::UInt8ToString_DataSize(chr_length);
        fprintf(stderr, "Chromosome: %s, id: %d, length: %s\n", chr.name, chr_id, size.c_str());
        chr.id = chr_id;
        chr.length = chr_length;
    }

    sam_close(in);
    pdqsort(chr_list.begin(), chr_list.end(), [](const SnpChrInfo& x, const SnpChrInfo& y) { return x.id < y.id; });
    return hdr;
}

static void
s_select_snp_bam_for_one_chr(SnpChrInfo& chr, const u8* cnt_table, vector<SnpBamInfo>& bam_list)
{
    HBN_LOG("Extract SNP BAM records for %d:%s:%d", chr.id, chr.name, chr.length);
    char path[HBN_MAX_PATH_LEN];
    s_make_chr_bam_path(chr.name, path);
    samFile* in = sam_open(path, "rb");
    if (!in) HBN_ERR("Could not open file %s for reading", path);
    hts_set_threads(in, num_threads);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", path);

    u8* cov_stats = new u8[chr.length]; fill(cov_stats, cov_stats + chr.length, 0);

    char pitagname[2] = { 'p', 'i' };
    string size;
    size_t chr_frags = 0, chr_bases = 0;
    size_t snp_data_size = chr.length; snp_data_size *= snp_frag_coverage;

    bam1_t* bam = nullptr;
    while (chr_bases < snp_data_size) {
        if (!bam) bam = bam_init1();
        int r = sam_read1(in, hdr, bam);
        if (r == -1) break;
        if (r < 0) HBN_ERR("FAIL at reading BAM record from %s", path);

        int frag_size = bam->core.l_qseq;
        int mapQ = bam->core.qual;
        uint8_t* pitag = bam_aux_get(bam, pitagname);
        if (!pitag) HBN_ERR("Could not retrieve alignment identity from %c%c BAM tag", pitagname[0], pitagname[1]);
        double pi = bam_aux2f(pitag);         
        r = (frag_size >= min_frag_size) && (mapQ >= min_mapQ) && (pi >= min_identity);
        if (!r) continue;   

        if (s_bam_cov_is_full(bam, cnt_table, snp_frag_coverage, cov_stats)) continue;

        bam_list.emplace_back(bam->core.tid, bam->core.pos, bam);
        bam = nullptr;
        ++chr_frags;  
        chr_bases += frag_size;

        if (chr_bases >= snp_data_size) break;       
    }
    if (bam) bam_destroy1(bam); 
    bam = nullptr;

    delete[] cov_stats;

    sam_hdr_destroy(hdr);
    sam_close(in);
    size = NStr::UInt8ToString_DataSize(chr_bases);
    const double cov = 1.0 * chr_bases / chr.length;
    HBN_LOG("Extract %zu BAM records (%s, %gX)", chr_frags, size.c_str(), cov);
}

int select_chr_snp_bam_main(int argc, char* argv[])
{
    vector<SnpChrInfo> chr_list;
    if (!parse_arguments(argc, argv, chr_list)) {
        dump_usage(argc, argv);
        exit (EXIT_FAILURE);
    }
    dump_parameters(chr_list);

    sam_hdr_t* hdr = s_extract_chr_id_and_length(chr_list);
    samFile* out = sam_open(snp_bam_path, "wb");
    if (!out) HBN_ERR("Could not open file %s for writing", snp_bam_path);
    hts_set_threads(out, num_threads);
    if (sam_hdr_write(out, hdr)) HBN_ERR("Could not write BAM header to %s", snp_bam_path);

    u8* cnt_table = build_cnt_table();
    vector<SnpBamInfo> bam_list;
    for (auto& chr : chr_list) {
        s_select_snp_bam_for_one_chr(chr, cnt_table, bam_list);
        pdqsort(bam_list.begin(), bam_list.end());
        for (auto& bi : bam_list) {
            if (sam_write1(out, hdr, bi.bam) < 0) HBN_ERR("Could not write BAM record to %s", snp_bam_path);
        }
        for (auto& bi : bam_list) bam_destroy1(bi.bam);
        bam_list.clear();
    }
    delete[] cnt_table;
    sam_hdr_destroy(hdr);
    sam_close(out);

    HBN_LOG("Create index for %s", snp_bam_path);
    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s.bai", snp_bam_path);
    int r = sam_index_build3(snp_bam_path, path, 0, num_threads);
    if (r) HBN_ERR("FAIL at creating index for %s", snp_bam_path);
    HBN_LOG("Done");

    return 0;
}
