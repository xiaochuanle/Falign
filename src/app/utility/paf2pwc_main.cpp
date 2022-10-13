#include "paf_reader.hpp"

using namespace std;

static bool 
s_parse_arguments(int argc, char* argv[], 
    int* map_q, double* pi, 
    int* frag_size,
    const char** paf_path)
{
    hbn_assert(argc >= 2);
    *paf_path = nullptr;
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

        HBN_ERR("Unrecognised option '%s'", argv[i]);
    }

    if (i + 1 != argc) return false;
    *paf_path = argv[i];
    return true;
}

static void 
s_dump_usage(int argc, char* argv[])
{
    hbn_assert(argc >= 2);
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] paf_path\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  -q <Integer>      Minimum mapping quality\n");
    fprintf(stderr, "  -p <Real>         Minimum percentage of identity of alignments\n");
    fprintf(stderr, "  -l <Integer>      Minimum length of fragment\n");
}

static void 
s_dump_pairwise_contact_salsa2(vector<PAF>& frag_align_list, 
    const char* read_name, vector<string>& ref_name_list,
    int mapping_quality, double pi, int frag_size, kstring_t* out)
{
    size_t n = frag_align_list.size();
    int contact_id = 0;
    char strand;
    for (size_t i = 0; i < n; ++i) {
        PAF& fai = frag_align_list[i];
        if (fai.map_q < mapping_quality) continue;
        if (fai.qend - fai.qoff < frag_size) continue;
        const char* ref_i_name = ref_name_list[fai.sname_idx].c_str();
        for (size_t j = i + 1; j < n; ++j) {
            PAF& faj = frag_align_list[j];
            if (faj.map_q < mapping_quality) continue;
            if (faj.qend - faj.qoff < frag_size) continue;
            const char* ref_j_name = ref_name_list[faj.sname_idx].c_str();
            //if (fai.sname == faj.sname) continue;

            ksprintf(out, "%s", ref_i_name);
            ksprintf(out, "\t");
            ksprintf(out, "%d", fai.soff);
            ksprintf(out, "\t");
            ksprintf(out, "%d", fai.send);
            ksprintf(out, "\t");
            ksprintf(out, "read%s_%d/1", read_name, contact_id);
            ksprintf(out, "\t");
            ksprintf(out, "%d", fai.map_q);
            strand = '+';//(fai.qdir == FWD) ? '+' : '-';
            ksprintf(out, "\t");
            ksprintf(out, "%c", strand);
            ksprintf(out, "\n");

            ksprintf(out, "%s", ref_j_name);
            ksprintf(out, "\t");
            ksprintf(out, "%d", faj.soff);
            ksprintf(out, "\t");
            ksprintf(out, "%d", faj.send);
            ksprintf(out, "\t");
            ksprintf(out, "read%s_%d/2", read_name, contact_id);
            ksprintf(out, "\t");
            ksprintf(out, "%d", faj.map_q);
            strand = '-';//(faj.qdir == FWD) ? '+' : '-';
            ksprintf(out, "\t");
            ksprintf(out, "%c", strand);
            ksprintf(out, "\n");
            ++contact_id;
        }
    }
}

static void 
s_dump_pairwise_contact_3ddna(vector<PAF>& frag_align_list, 
    const char* read_name, vector<string>& ref_name_list,
    int mapping_quality, double pi, int frag_size, kstring_t* out)
{
    size_t n = frag_align_list.size();
    for (size_t i = 0; i < n; ++i) {
        PAF& fai = frag_align_list[i];
        if (fai.map_q < mapping_quality) continue;
        if (fai.qend - fai.qoff < frag_size) continue;
        const char* ref_i_name = ref_name_list[fai.sname_idx].c_str();
        for (size_t j = i + 1; j < n; ++j) {
            PAF& faj = frag_align_list[j];
            if (faj.map_q < mapping_quality) continue;
            if (faj.qend - faj.qoff < frag_size) continue;
            const char* ref_j_name = ref_name_list[faj.sname_idx].c_str();

            ksprintf(out, "0"); /// strand 1
            kputc(' ', out);
            ksprintf(out, "%s", ref_i_name); /// chr1
            kputc(' ', out);
            ksprintf(out, "%d", fai.soff); /// chr1 pos
            kputc(' ', out);
            kputc('0', out); /// frag 1
            kputc(' ', out);

            ksprintf(out, "-1"); /// strand 2
            kputc(' ', out); 
            ksprintf(out, "%s", ref_j_name); /// chr2
            kputc(' ', out);
            ksprintf(out, "%d", faj.soff); /// chr2 pos
            kputc(' ', out);
            kputc('1', out); /// frag2
            kputc(' ', out);

            ksprintf(out, "%d", fai.map_q); /// map_q 1
            kputc(' ', out);
            kputc('.', out); /// cigar 1
            kputc(' ', out);
            kputc('.', out); /// frag 1 sequence
            kputc(' ', out);

            ksprintf(out, "%d", faj.map_q); /// map_q 2
            kputc(' ', out);
            kputc('.', out); // cigar 2
            kputc(' ', out);
            kputc('.', out); /// frag 2 sequence
            kputc(' ', out);

            ksprintf(out, "%s_1", read_name); /// frag1 name
            kputc(' ', out);
            ksprintf(out, "%s_2", read_name); /// frag2 name
            kputc('\n', out);
        }
    }
}

typedef void (*pwc_dump_func)(vector<PAF>& frag_align_list, 
    const char* read_name, vector<string>& ref_name_list,
    int mapping_quality, double pi, int frag_size, kstring_t* out);

int paf_to_pairwise_contact_bed_main(int argc, char* argv[], const char* target)
{
    int mapq = 0;
    double pi = 0.0;
    int frag_size = 0;
    const char* paf_path = nullptr;
    if (!s_parse_arguments(argc, argv, &mapq, &pi, &frag_size, &paf_path)) {
        s_dump_usage(argc, argv);
        return 1;
    }
    vector<PAF> paf_list;
    vector<string> ref_name_list;
    string read_name;
    PAFReader paf_in(paf_path);
    int read_idx = 0;
    ks_dinit(out_buf);
    pwc_dump_func output = (strcmp(target, "salsa2") == 0) ? s_dump_pairwise_contact_salsa2 : s_dump_pairwise_contact_3ddna;
    while (paf_in.read_next_paf_list(paf_list, read_name, ref_name_list)) {
        ks_clear(out_buf);
        output(paf_list, read_name.c_str(), ref_name_list, mapq, pi, frag_size, &out_buf);
        hbn_fwrite(ks_s(out_buf), 1, ks_size(out_buf), stdout);
        ++read_idx;
        if ((read_idx % 100000) == 0) fprintf(stderr, "%10d queries processed\n", read_idx);
    }
    ks_destroy(out_buf);

    return 0;
}