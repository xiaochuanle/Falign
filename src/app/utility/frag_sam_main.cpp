#include "../utility/sam_reader.hpp"
#include "../../algo/hbn_traceback_aux.h"
#include "../../corelib/cstr_util.h"
#include "frag_id.hpp"
#include <cstdio>

using namespace std;

static void
s_copy_sam_hdr(const char* input_sam_path, FILE* out)
{
    HbnLineReader* line_reader = HbnLineReaderNew(input_sam_path);
    while (!HbnLineReaderAtEof(line_reader)) {
        HbnLineReaderReadOneLine(line_reader);
        const char* s = ks_s(line_reader->line);
        const int sl = ks_size(line_reader->line);
        if (s[0] != '@') break;
        hbn_fwrite(s, 1, sl, out);
        ::fprintf(out, "\n");
    }
    HbnLineReaderFree(line_reader);
}

static void
s_dump_one_frag_sam(map<string, int>& subject_name2id, FragmentAlignment& align, string& sam, string& align_strings, int read_id, int frag_id, FILE* out)
{
    ks_dinit(cigar);
    string subseq;
    const char* qas = align_strings.c_str() + align.qas_offset;
    const char* sas = align_strings.c_str() + align.sas_offset;
    const int as_size = align.as_size;
    int xi = 0;
    while (xi < as_size) {
        int num = 0;
        char op = 'N';
        int xj = xi + 1;
        if (qas[xi] != GAP_CHAR && sas[xi] != GAP_CHAR) {
                op = 'M';
                while (xj < as_size && qas[xj] != GAP_CHAR && sas[xj] != GAP_CHAR) ++xj;
        } else if (qas[xi] == GAP_CHAR) {
                op = 'D';
                while (xj < as_size && qas[xj] == GAP_CHAR) ++xj;
        } else {
                hbn_assert(sas[xi] == GAP_CHAR);
                op = 'I';
                while (xj < as_size && sas[xj] == GAP_CHAR) ++xj;
        }
        num = xj - xi;
        xi = xj;
        ksprintf(&cigar, "%d%c", num, op);
    }
    for (xi = 0; xi < as_size; ++xi) {
	if (qas[xi] != GAP_CHAR) {
		subseq += qas[xi];
		hbn_assert(qas[xi] == 'A' || qas[xi] == 'C' || qas[xi] == 'G' || qas[xi] == 'T', "xi = %d, qas[xi] = %d, sas[xi] = %d as_size = %d", xi, qas[xi], sas[xi], as_size);
	}
    }

    char id_buf[FRAD_ID_SIZE];
    auto pos = subject_name2id.find(align.sname);
    hbn_assert(pos != subject_name2id.end());
    frag_id_to_string(read_id, frag_id, pos->second, align.soff, id_buf);
    /// 1) query name
    fprintf(out, "%s_%s", align.qname.c_str(), id_buf);
    fprintf(out, "\t");
    /// 2) flag
    int sam_flag = (align.qdir == FWD) ? 0 : 0x10;
    fprintf(out, "%d", sam_flag);
    fprintf(out, "\t");
    /// 3) sname
    fprintf(out, "%s", align.sname.c_str());
    fprintf(out, "\t");
    /// 4) left most subject position (1-based)
    fprintf(out, "%d", align.soff + 1);
    fprintf(out, "\t");
    /// 5) map_q
    fprintf(out, "%d", align.map_q);
    fprintf(out, "\t");
    /// 6) cigar
    fprintf(out, "%s", ks_s(cigar));
    fprintf(out, "\t");
    /// 7) reference name of the mate/next read
    fprintf(out, "*");
    fprintf(out, "\t");
    /// 8) position of the mate/next read
    fprintf(out, "0");
    fprintf(out, "\t");
    /// 9) observed template LENgth
    fprintf(out, "0");
    fprintf(out, "\t");
    /// 10) segment SEQuence
    fprintf(out, "%s", subseq.c_str());
    fprintf(out, "\t");
    /// 11) ASCII of Phred-scaled base QUALity+33
    fprintf(out, "*");
    fprintf(out, "\t");

    const char* md = strstr(sam.c_str(), "MD:Z:");
    for (int i = 0; !isspace(md[i]); ++i) fprintf(out, "%c", md[i]);
    fprintf(out, "\n");

    ks_destroy(cigar);
}

static void
s_dump_one_align_list(vector<FragmentAlignment>& align_list,
    vector<string>& sam_list,
    string& align_strings,
    map<string, int>& sname2id_map,
    const int read_id,
    const int map_q,
    const double pi,
    const int frag_size,
    FILE* sam_out)
{
    int n_frag = align_list.size();
    /// frag sam
    for (int i = 0; i < n_frag; ++i) {
        FragmentAlignment& align = align_list[i];
        if (align.map_q < map_q) continue;
        if (align.pi < pi) continue;
        if (align.qend - align.qoff < frag_size) continue;
        s_dump_one_frag_sam(sname2id_map, align, sam_list[i], align_strings, read_id, i, sam_out);
    }
}

static void 
s_dump_usage(int argc, char* argv[])
{
    hbn_assert(argc >= 2);
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] sam-path frag-sam-path\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  -q <Integer>      Minimum mapping quality\n");
    fprintf(stderr, "  -p <Real>         Minimum percentage of identity of alignments\n");
    fprintf(stderr, "  -l <Integer>      Minimum length of fragment\n");
    fprintf(stderr, "  -g <Integer>      Genome size\n");
    fprintf(stderr, "  -c <Integer>      Coverage\n");
}

static bool 
s_parse_arguments(int argc, char* argv[], 
    int* map_q, 
    double* pi, 
    int* frag_size,
    size_t* genome_size,
    int* coverage,
    const char** sam_path)
{
    *map_q = 0;
    *pi = 0.0;
    *frag_size = 0;
    *genome_size = 0;
    *coverage = 0;
    *sam_path = nullptr;

    hbn_assert(argc >= 2);
    int i = 2;
    while (i < argc) {
        if (argv[i][0] != '-') break;

        if (strcmp(argv[i], "-q") == 0) {
            if (i + 1 >= argc) HBN_ERR("Argument to option '%s' is missing", argv[i]);
            *map_q = atoi(argv[i+1]);
            i += 2;
            continue;
        }     

        if (strcmp(argv[i], "-p") == 0) {
            if (i + 1 >= argc) HBN_ERR("Argument to option '%s' is missing", argv[i]);
            *pi = atof(argv[i+1]);
            i += 2;
            continue;
        }     

        if (strcmp(argv[i], "-l") == 0) {
            if (i + 1 >= argc) HBN_ERR("Argument to option '%s' is missing", argv[i]);
            *frag_size = atoi(argv[i+1]);
            i += 2;
            continue;
        }                

        if (strcmp(argv[i], "-g") == 0) {
            if (i + 1 >= argc) HBN_ERR("Argument to option '%s' is missing", argv[i]);
            *genome_size = atoll(argv[i+1]);
            i += 2;
            continue;            
        }

        if (strcmp(argv[i], "-c") == 0) {
            if (i + 1 >= argc) HBN_ERR("Argument to option '%s' is missing", argv[i]);
            *coverage = atoi(argv[i+1]);
            i += 2;
            continue;             
        }

        HBN_ERR("Unrecognised option '%s'", argv[i]);
    }

    if (i + 1 != argc) return false;
    *sam_path = argv[i];
    return true;
}

int make_frag_sam_main(int argc, char* argv[])
{
    const char* map_sam_path = nullptr;
    int map_q = 0;
    double pi = 0.0;
    int frag_size = 0;
    size_t genome_size = 0;
    int cov = 0;
    if (!s_parse_arguments(argc, argv, &map_q, &pi, &frag_size, &genome_size, &cov, &map_sam_path)) {
        s_dump_usage(argc, argv);
        return 1;
    }
    size_t target_size = (genome_size == 0 || cov == 0) ? U64_MAX : genome_size * cov;

    fprintf(stderr, "map_q: %d\n", map_q);
    fprintf(stderr, "identity: %g\n", pi);
    fprintf(stderr, "frag size: %d\n", frag_size);
    fprintf(stderr, "genome size: %zu\n", genome_size);
    fprintf(stderr, "coverage: %d\n", cov);

    FILE* sam_out = stdout;
    s_copy_sam_hdr(map_sam_path, sam_out);
    map<string, int> sname2id_map;
    load_subject_id_from_sam(map_sam_path, sname2id_map);

    vector<FragmentAlignment> align_list;
    vector<string> sam_list;
    string align_strings;
    bool is_complete_map = false;
    SamLoader sam(map_sam_path);
    int read_id = 0;
    size_t dumpped_size = 0;
    while (sam.load_next_align_list(align_list, sam_list, align_strings, is_complete_map)) {
        s_dump_one_align_list(align_list, sam_list, align_strings, sname2id_map, read_id, map_q, pi, frag_size, sam_out);
        ++read_id;
        dumpped_size += align_list.front().qsize;
        if (dumpped_size >= target_size) break;
        if ((read_id % 100000) == 0) fprintf(stderr, "%10d queries processed\n", read_id);
    }
    char buf[64];
    u64_to_string_datasize(dumpped_size, buf);
    HBN_LOG("dump frag sam records for %d reads (%s)", read_id, buf);
    return 0;
}