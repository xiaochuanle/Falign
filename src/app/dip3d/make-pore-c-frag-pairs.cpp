#include "../../corelib/frag_id.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "split_string_by_char.hpp"

#include <algorithm>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

using namespace std;

static const char kDelim = 32;

struct LocalFragAlign
{
    int read_id;
    int frag_id;
    size_t frag_offset;
    int frag_size;
    int subject_id;
    int subject_offset;
};

static char* g_frag_list = nullptr;
static vector<LocalFragAlign> g_align_list;
static vector<pair<int, int>> g_read_frag_info_list;
static int g_single_frags = 0;
static int g_frag_pairs = 0;

static inline bool
s_is_end_of_line(int c)
{
    return (c == '\r') || (c == '\n');
}

static void
s_load_frags(const char* frag_file_path)
{
    const size_t file_size = hbn_file_size(frag_file_path);
    g_frag_list = new char[file_size];
    hbn_dfopen(in, frag_file_path, "r");
    hbn_fread(g_frag_list, 1, file_size, in);
    hbn_fclose(in);

    vector<pair<const char*, int>> cols;
    size_t i = 0;
    int num_reads = 0;
    while (i < file_size) {
        size_t j = i + 1;
        while (j < file_size && (!s_is_end_of_line(g_frag_list[j]))) ++j;
        cols.clear();
        split_string_by_char(g_frag_list + i, j - i, kDelim, cols);
        int read_id = 0;
        int frag_id = 0;
        int subject_id = 0;
        int subject_offset = 0;
        //fprintf(stderr, "i = %zu, j = %zu, read name: %s, frag_id = %d\n", i, j, read_name.c_str(), frag_id);
        extract_frag_id_from_name(cols[1].first, cols[1].second, &read_id, &frag_id, &subject_id, &subject_offset);
        
        num_reads = max(num_reads, read_id);
        LocalFragAlign lfa;
        lfa.read_id = read_id;
        lfa.frag_id = frag_id;
        lfa.frag_offset = i;
        lfa.frag_size = j - i;
        lfa.subject_id = subject_id;
        lfa.subject_offset = subject_offset;
        g_align_list.push_back(lfa);
        i = j;
        while(i < file_size && isspace(g_frag_list[i])) ++i;
    }
    sort(g_align_list.begin(), g_align_list.end(),
	    [](const LocalFragAlign& a, const LocalFragAlign& b) { return (a.read_id < b.read_id) || (a.read_id == b.read_id && a.frag_id < b.frag_id); });

    ++num_reads;
    g_read_frag_info_list.resize(num_reads);
    for (auto& x : g_read_frag_info_list) { x.first = 0; x.second = 0; }
    i = 0;
    const size_t n_align = g_align_list.size();
    while (i < n_align) {
        size_t j = i + 1;
        while (j < n_align && g_align_list[i].read_id == g_align_list[j].read_id) ++j;
        g_read_frag_info_list[g_align_list[i].read_id].first = i;
        g_read_frag_info_list[g_align_list[i].read_id].second = j - i;
        i = j;
    }
    fprintf(stderr, "Load %zu fragment info from %s\n", g_align_list.size(), frag_file_path);
}

static void
s_make_and_dump_single_frag_hic(LocalFragAlign* align, FILE* out)
{
    ++g_single_frags;
    const char* frag = g_frag_list + align->frag_offset;
    vector<pair<const char*, int>> cols;
    split_string_by_char(frag, align->frag_size, kDelim, cols);
    int n_cols = cols.size();
    int date_type = 1; // 0 for normal, 1 for hic
    int mate_idx = -1; // -1 if no mate 2
    int insert_size = -1; // -1 if no mate 2

    hbn_fwrite(cols[0].first, 1, cols[0].second, out);
    fprintf(out, " ");
    hbn_fwrite(cols[1].first, 1, cols[1].second, out);
    fprintf(out, " ");
    fprintf(out, "%d", date_type);
    fprintf(out, " ");
    fprintf(out, "%d", mate_idx);
    fprintf(out, " ");
    fprintf(out, "%d", insert_size);
    fprintf(out, " ");
    for (int i = 2; i < n_cols; ++i) {
        hbn_fwrite(cols[i].first, 1, cols[i].second, out);
        if (i < n_cols - 1) fprintf(out, " ");
    }
    fprintf(out, "\n");
}

static void
s_extract_var_id_info(vector<pair<const char*, int>>& cols, int& f_varid, int& l_varid)
{
    const int n_blk = atoi(cols[0].first);
    f_varid = atoi(cols[2].first);
    int blk_i = 0, col_i = 2;
    while (blk_i < n_blk) {
        l_varid = atoi(cols[col_i].first);
        if (cols[col_i+1].second > 1) l_varid += (cols[col_i+1].second - 1);
        ++blk_i;
        col_i += 2;
    }
}

static void
s_validate_varid(const char* i_frag, const int i_frag_size, vector<pair<const char*, int>>& i_cols, const char* j_frag, const int j_frag_size, vector<pair<const char*, int>>& j_cols)
{
    int i_f_varid = -1, i_l_varid = -1;
    s_extract_var_id_info(i_cols, i_f_varid, i_l_varid);
    int j_f_varid = -1, j_l_varid = -1;
    s_extract_var_id_info(j_cols, j_f_varid, j_l_varid);
    if (i_l_varid >= j_f_varid) {
        fprintf(stderr, "inconsistent frag pair: i_l_varid = %d, j_f_varid = %d\n", i_l_varid, j_f_varid);
        hbn_fwrite(i_frag, 1, i_frag_size, stderr);
        fprintf(stderr, "\n");
        hbn_fwrite(j_frag, 1, j_frag_size, stderr);
        fprintf(stderr, "\n");
        abort();
    }
}

static void
s_make_and_dump_multi_frag_hic(LocalFragAlign* align_list, int align_cnt, FILE* pair_frag_out)
{
    vector<pair<const char*, int>> i_cols;
    vector<pair<const char*, int>> j_cols;
    ostringstream os;
    string oss;
    for (int i = 0; i < align_cnt; ++i) {
        LocalFragAlign* ai = align_list + i;
        const char* fragi = g_frag_list + ai->frag_offset;
        i_cols.clear();
        split_string_by_char(fragi, ai->frag_size, kDelim, i_cols);
        int i_sites_num = atoi(i_cols[0].first);
        int i_cols_cnt = i_cols.size();
        for (int j = i + 1; j < align_cnt; ++j) {
            LocalFragAlign* aj = align_list + j;
            if (ai->subject_id != aj->subject_id) continue;
            const char* fragj = g_frag_list + aj->frag_offset;
            j_cols.clear();
            split_string_by_char(fragj, aj->frag_size, kDelim, j_cols);
            int j_sites_mum = atoi(j_cols[0].first);
            int j_cols_cnt = j_cols.size();

            int sites_num = i_sites_num + j_sites_mum;
            int date_type = 1;
            int insert_size = abs(ai->subject_offset - aj->subject_offset);
            os.str("");
            if (ai->subject_offset < aj->subject_offset) { // ai -> aj
                s_validate_varid(fragi, ai->frag_size, i_cols, fragj, aj->frag_size, j_cols);

                os << sites_num;
                os << kDelim << NStr::CTempString(i_cols[1].first, i_cols[1].second);
                const char* j_suffix = j_cols[1].first + j_cols[1].second - 3;
                os << '-' << NStr::CTempString(j_suffix, 3);
                os << kDelim << date_type;
                os << kDelim << NStr::CTempString(j_cols[2].first, j_cols[2].second);
                os << kDelim << insert_size;

                os << kDelim;
                for (int t = 2; t < i_cols_cnt - 1; ++t) {
                    os << NStr::CTempString(i_cols[t].first, i_cols[t].second);
                    os << kDelim;
                }
                for (int t = 2; t < j_cols_cnt - 1; ++t) {
                    os << NStr::CTempString(j_cols[t].first, j_cols[t].second);
                    if (t < j_cols_cnt - 2) os << kDelim;
                }

                os << kDelim;
                os << NStr::CTempString(i_cols.back().first, i_cols.back().second);
                os << NStr::CTempString(j_cols.back().first, j_cols.back().second);
                os << '\n';
            } else { // aj -> ai
                s_validate_varid(fragj, aj->frag_size, j_cols, fragi, ai->frag_size, i_cols);

                os << sites_num;
                os << kDelim << NStr::CTempString(i_cols[1].first, i_cols[1].second);
                const char* j_suffix = j_cols[1].first + j_cols[1].second - 3;
                os << '-' << NStr::CTempString(j_suffix, 3);
                os << kDelim << date_type;
                os << kDelim << NStr::CTempString(i_cols[2].first, i_cols[2].second);
                os << kDelim << insert_size;

                os << kDelim;
                for (int t = 2; t < j_cols_cnt - 1; ++t) {
                    os << NStr::CTempString(j_cols[t].first, j_cols[t].second);
                    os << kDelim;
                }
                for (int t = 2; t < i_cols_cnt - 1; ++t) {
                    os << NStr::CTempString(i_cols[t].first, i_cols[t].second);
                    if (t < i_cols_cnt - 2) os << kDelim;
                }

                os << kDelim;
                os << NStr::CTempString(j_cols.back().first, j_cols.back().second);
                os << NStr::CTempString(i_cols.back().first, i_cols.back().second);
                os << '\n';         
            }
            oss = os.str();
            hbn_fwrite(oss.c_str(), 1, oss.size(), pair_frag_out);
            ++g_frag_pairs;
        }
    }
}

static void
s_make_and_dump_pair_frag(const char* input, FILE* out)
{
    const size_t n_read = g_read_frag_info_list.size();
    fprintf(stderr, "number of reads: %zu\n", n_read);
    for (size_t i = 0; i < n_read; ++i) {
        LocalFragAlign* a = g_align_list.data() + g_read_frag_info_list[i].first;
        int c = g_read_frag_info_list[i].second;
        sort(a, a + c, [](const LocalFragAlign& a, const LocalFragAlign& b) { return a.subject_id < b.subject_id; });
        int x = 0;
        while (x < c) {
            int y = x + 1;
            while (y < c && a[x].subject_id == a[y].subject_id) ++y;
            if (y - x == 1) {
                s_make_and_dump_single_frag_hic(a + x, out);
            } else {
                s_make_and_dump_multi_frag_hic(a + x, y - x, out);
            }
            x = y;
        }
    }
    HBN_LOG("%s", input);
    HBN_LOG("single fragments: %d", g_single_frags);
    HBN_LOG("fragment pairs: %d", g_frag_pairs);
}

int make_pore_c_frag_pair_main(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s single_frag_path pair_frag_path\n", argv[0], argv[1]);
        return 1;
    }
    const char* frag_file_path = argv[2];
    const char* pair_frag_path = argv[3];
    s_load_frags(frag_file_path);
    hbn_dfopen(out, pair_frag_path, "w");
    s_make_and_dump_pair_frag(frag_file_path, out);
    hbn_fclose(out);
    delete[] g_frag_list;
    return 0;
}