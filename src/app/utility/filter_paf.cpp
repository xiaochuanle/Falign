#include "paf_reader.hpp"

using namespace std;

struct FragMapQStats
{
    int cnt_0_5;
    int cnt_5_10;
    int cnt_10_20;
    int cnt_20_30;
    int cnt_30_40;
    int cnt_40_50;
    int cnt_50_60;
    int cnt;

    FragMapQStats() {
        cnt_0_5 = 0;
        cnt_5_10 = 0;
        cnt_10_20 = 0;
        cnt_20_30 = 0;
        cnt_30_40 = 0;
        cnt_40_50 = 0;
        cnt_50_60 = 0;
        cnt = 0;
    }

    void add(int mapq) {
        if (mapq < 5) {
            ++cnt_0_5;
        } else if (mapq < 10) {
            ++cnt_5_10;
        } else if (mapq < 20) {
            ++cnt_10_20;
        } else if (mapq < 30) {
            ++cnt_20_30;
        } else if (mapq < 40) {
            ++cnt_30_40;
        } else if (mapq < 50) {
            ++cnt_50_60;
        } else {
            ++cnt_50_60;
        }
        ++cnt;
    }

    void stats() {
        double p = 1.0 * cnt_0_5 / cnt;
        fprintf(stderr, "[0, 5): %g\n", p);

        p = 1.0 * cnt_5_10 / cnt;
        fprintf(stderr, "[5, 10): %g\n", p);

        p = 1.0 * cnt_10_20 / cnt;
        fprintf(stderr, "[10, 20): %g\n", p);

        p = 1.0 * cnt_20_30 / cnt;
        fprintf(stderr, "[20, 30): %g\n", p);

        p = 1.0 * cnt_30_40 / cnt;
        fprintf(stderr, "[30, 40): %g\n", p);

        p = 1.0 * cnt_40_50 / cnt;
        fprintf(stderr, "[40, 50): %g\n", p);

        p = 1.0 * cnt_50_60 / cnt;
        fprintf(stderr, "[50, 60]: %g\n", p);
    }
};

static void 
s_dump_usage(int argc, char* argv[])
{
    hbn_assert(argc >= 2);
    fprintf(stderr, "USAGE:\n");
    fprintf(stderr, "%s %s [OPTIONS] paf-path\n", argv[0], argv[1]);

    fprintf(stderr, "\n");
    fprintf(stderr, "OPTIONAL ARGUMENTS\n");
    fprintf(stderr, "  -q <Integer>      Minimum mapping quality\n");
    fprintf(stderr, "  -p <Real>         Minimum percentage of identity of alignments\n");
    fprintf(stderr, "  -l <Integer>      Minimum length of fragment\n");
}

static bool 
s_parse_arguments(int argc, char* argv[], 
    int* map_q, 
    double* pi, 
    int* frag_size,
    const char** paf_path)
{
    *map_q = 0;
    *pi = 0.0;
    *frag_size = 0;
    *paf_path = nullptr;

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

        HBN_ERR("Unrecognised option '%s'", argv[i]);
    }

    if (i + 1 != argc) return false;
    *paf_path = argv[i];
    return true;
}

int filter_paf_main(int argc, char* argv[])
{
    const char* paf_path = nullptr;
    int map_q = 0;
    double pi = 0.0;
    int frag_size = 0;
    if (!s_parse_arguments(argc, argv, &map_q, &pi, &frag_size, &paf_path)) {
        s_dump_usage(argc, argv);
        return 1;
    }

    fprintf(stderr, "map_q: %d\n", map_q);
    fprintf(stderr, "identity: %g\n", pi);
    fprintf(stderr, "frag size: %d\n", frag_size);

    PAFReader paf = PAFReader(paf_path);
    vector<PAF> paf_list;
    string read_name;
    vector<string> ref_name_list;
    vector<string> line_list;
    FragMapQStats stats;
    while (paf.read_next_paf_list(paf_list, read_name, ref_name_list, line_list)) {
        int n = paf_list.size();
        for (int i = 0; i < n; ++i) {
            PAF& p = paf_list[i];
            stats.add(p.map_q);
            if (p.map_q < map_q) continue;
            if (p.qend - p.qoff < frag_size) continue;
            fprintf(stdout, "%s\n", line_list[i].c_str());
        }
    }
    stats.stats();

    return 0;
}