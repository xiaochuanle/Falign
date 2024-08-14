#include "split_string_by_char.hpp"

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#include <cstring>

using namespace std;

static vector<int> frag_length_list;
static vector<int> read_frag_list;
size_t num_reads = 0;
size_t num_contacts = 0;

static void
s_add_one_frag_list(vector<int>& length_list)
{
    frag_length_list.insert(frag_length_list.end(), length_list.begin(), length_list.end());
    read_frag_list.push_back(length_list.size());
    ++num_reads;
    int n = length_list.size();
    n = n * (n-1)/2;
    num_contacts += n;
}

static void
s_compute_mean_and_median(vector<int>& xlist, double& mean, int& median)
{
    size_t sum = 0;
    sort(xlist.begin(), xlist.end());
    for (auto x : xlist) sum += x;
    size_t N = xlist.size();
    mean = 1.0 * sum / N;
    median = xlist[N/2];
}

int paf_frag_stats_main(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s paf\n", argv[0], argv[1]);
        exit (1);
    }

    size_t num_paf = 0;
    vector<pair<const char*, int>> cols;
    ifstream in(argv[2]);
    string line;
    string last_qname;
    vector<int> length_list;
    while (getline(in, line)) {
	    if (!strstr(line.c_str(), "PAH58618")) continue;
        cols.clear();
        split_string_by_char(line.c_str(), line.size(), '\t', cols);
        if (last_qname.size() != cols[0].second || strncmp(last_qname.c_str(), cols[0].first, cols[0].second)) {
            if (!length_list.empty()) s_add_one_frag_list(length_list);
            length_list.clear();
            last_qname.assign(cols[0].first, cols[0].second);
        }
        int from = atoi(cols[2].first);
        int to = atoi(cols[3].first);
        length_list.push_back(to - from);
        ++num_paf;
        if ((num_paf % 1000000) == 0) fprintf(stderr, "%zu paf\n", num_paf);
    }

    fprintf(stderr, "reads: %zu\n", num_reads);
    fprintf(stderr, "contacts: %zu\n", num_contacts);

    double mean;
    int median;
    s_compute_mean_and_median(frag_length_list, mean, median);
    fprintf(stderr, "frag length, mean: %g, median: %d\n", mean, median);

    s_compute_mean_and_median(read_frag_list, mean, median);
    fprintf(stderr, "read frag, mean: %g, median: %d\n", mean, median);

    return 0;
}
