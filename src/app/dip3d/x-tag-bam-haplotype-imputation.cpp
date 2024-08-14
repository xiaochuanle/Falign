#include "x-tag-bam-haplotype-imputation.hpp"

#include "../../corelib/pdqsort.h"

#include <map>
#include <mutex>
#include <set>
#include <vector>

using namespace std;

class FragHapImputationThreadWorkdata
{
public:
    FragHapImputationThreadWorkdata(BamTagOptions* options, FragHapInfo* fhia, size_t fhic) {
        M_options = options;

        size_t i = 0;
        while (i < fhic) {
            size_t j = i + 1;
            while (j < fhic && fhia[i].read_id == fhia[j].read_id) ++j;

            M_read_frags.emplace_back(fhia + i, j - i);

            i = j;
        }

        M_read_idx = 0;
    }

    BamTagOptions* options() {
        return M_options;
    }

    bool get_next_read(FragHapInfo** a, int* c) {
        lock_guard<mutex> _(M_read_mutex);
        if (M_read_idx >= M_read_frags.size()) return false;
        *a = M_read_frags[M_read_idx].first;
        *c = M_read_frags[M_read_idx].second;
        ++M_read_idx;
        return true;
    }

private:
    BamTagOptions*                  M_options;
    vector<pair<FragHapInfo*, int>> M_read_frags;
    size_t                          M_read_idx;
    mutex                           M_read_mutex;
};

static void
s_merge_subfrag_list(vector<map<int, FragHapInfo*>>& subfrag_list)
{
    int n = subfrag_list.size();
    for (int i = 1; i < n; ++i) {
        if (subfrag_list[i].empty()) continue;
        for (int j = 0; j < i; ++j) {
            if (subfrag_list[j].empty()) continue;
            bool has_overlapped_fid = false;
            for (auto& jx : subfrag_list[j]) {
                if (subfrag_list[i].find(jx.first) != subfrag_list[i].end()) {
                    has_overlapped_fid = true;
                    break;
                }
            }
            if (!has_overlapped_fid) continue;
            for (auto& jx : subfrag_list[j]) {
                if (subfrag_list[i].find(jx.first) != subfrag_list[i].end()) continue;
                subfrag_list[i].insert(jx);
            }
            subfrag_list[j].clear();
        }
    }

    for (int i = 1; i < n; ++i) {
        if (subfrag_list[i].empty()) continue;
        for (int j = 0; j < i; ++j) {
            if (subfrag_list[j].empty()) continue;
            bool has_bridged_fid = false;
            for (auto& ix : subfrag_list[i]) {
                if (subfrag_list[j].find(ix.first+1) != subfrag_list[j].end() && subfrag_list[j].find(ix.first-1) != subfrag_list[j].end()) {
                    has_bridged_fid = true;
                    break;
                }
            }
            if (!has_bridged_fid) continue;
            for (auto& jx : subfrag_list[j]) {
                if (subfrag_list[i].find(jx.first) != subfrag_list[i].end()) continue;
                subfrag_list[i].insert(jx);
            }
            subfrag_list[j].clear();
        }
    }

    int m = 0;
    for (int i = 0; i < n; ++i) {
        if (subfrag_list[i].empty()) continue;
        if (i > m) {
            subfrag_list[m] = subfrag_list[i];
        }
        ++m;
    }
    subfrag_list.resize(m);
}

static void
s_cut_subfrag_list(vector<map<int, FragHapInfo*>>& h1_subfrag_list, vector<map<int, FragHapInfo*>>& h2_subfrag_list)
{
    map<int, FragHapInfo*> h1_frags, h2_frags;
    for (auto& v1 : h1_subfrag_list) {
        for (auto& x1 : v1) {
            if (h1_frags.find(x1.first) != h1_frags.end()) continue;
            h1_frags.insert(x1);
        }
    }
    for (auto& v2 : h2_subfrag_list) {
        for (auto& x2 : v2) {
            if (h2_frags.find(x2.first) != h2_frags.end()) continue;
            h2_frags.insert(x2);
        }
    }

    map<int, FragHapInfo*> com_frags, com1_frags, com2_frags;
    for (auto& x1 : h1_frags) {
        if (h2_frags.find(x1.first) == h2_frags.end()) continue;
        com_frags.insert(x1);
        if (x1.second->hp == 1) com1_frags.insert(x1);
        if (x1.second->hp == 2) com2_frags.insert(x1);
    }

    // filter frag in both h1 and h2 subread, with identified hp
    map<int, FragHapInfo*> frags;
    if (!com1_frags.empty()) {
        for (auto& v2 : h2_subfrag_list) {
            frags.clear();
            for (auto& x2 : v2) {
                if (com1_frags.find(x2.first) != com1_frags.end()) continue;
                frags.insert(x2);
            }
            v2 = frags;
        }
    }
    if (!com2_frags.empty()) {
        for (auto& v1 : h1_subfrag_list) {
            frags.clear();
            for (auto& x1 : v1) {
                if (com2_frags.find(x1.first) != com2_frags.end()) continue;
                frags.insert(x1);
            }
            v1 = frags;
        }
    }

    h1_frags.clear();
    for (auto& v1 : h1_subfrag_list) {
        for (auto& x1 : v1) {
            if (h1_frags.find(x1.first) != h1_frags.end()) continue;
            h1_frags.insert(x1);
        }
    }
    h2_frags.clear();
    for (auto& v2 : h2_subfrag_list) {
        for (auto& x2 : v2) {
            if (h2_frags.find(x2.first) != h2_frags.end()) continue;
            h2_frags.insert(x2);
        }
    }
    com_frags.clear();
    for (auto& x1 : h1_frags) {
        if (h2_frags.find(x1.first) == h2_frags.end()) continue;
        com_frags.insert(x1);
    }
    // filter frag in both h1 and h2 subread, unphased
    if (!com_frags.empty()) {
        for (auto& v1 : h1_subfrag_list) {
            frags.clear();
            for (auto& x1 : v1) {
                if (com_frags.find(x1.first) != com_frags.end()) continue;
                frags.insert(x1);
            }
            v1 = frags;
        }
        for (auto& v2 : h2_subfrag_list) {
            frags.clear();
            for (auto& x2 : v2) {
                if (com_frags.find(x2.first) != com_frags.end()) continue;
                frags.insert(x2);
            }
            v2 = frags;
        }
    }
}

static void
s_impute_one_list_large_dist(FragHapInfo* fhia, int fhic)
{
    int hp0 = 0, hp1 = 0, hp2 = 0;
    for (int i = 0; i < fhic; ++i) {
        if (fhia[i].hp == 1) ++hp1;
        if (fhia[i].hp == 2) ++hp2;
        if (fhia[i].hp == 0) ++hp0;
    }
    if (hp0 == 0) return;

    if (hp1 > 0 && hp2 == 0) {
	    for (int i = 0; i < fhic; ++i) if (fhia[i].hp == 0) fhia[i].hp = 1;
	    return;
    }

    if (hp1 == 0 && hp2 > 0) {
	    for (int i = 0; i < fhic; ++i) if (fhia[i].hp == 0) fhia[i].hp = 2;
	    return;
    }

    double p1 = 1.0 * hp1 / fhic;
    if (p1 >= 0.8) {
        for (int i = 0; i < fhic; ++i) {
            if (fhia[i].hp == 0) fhia[i].hp = 1;
        }
    }

    double p2 = 1.0 * hp2 / fhic;
    if (p2 >= 0.8) {
        for (int i = 0; i < fhic; ++i) {
            if (fhia[i].hp == 0) fhia[i].hp = 2;
        }
    }
}

static void
s_impute_one_list(FragHapInfo* fhia, int fhic, const int res_dist, const int adj_dist)
{
    vector<map<int, FragHapInfo*>> h1_subfrag_list;
    vector<map<int, FragHapInfo*>> h2_subfrag_list;
    map<int, FragHapInfo*> subfrag_list;
    map<int, FragHapInfo*> cf, bf;

    for (int i = 0; i < fhic; ++i) {
        if (fhia[i].hp == 0) continue;
        cf.clear();
        bf.clear();
        for (int j = 0; j < fhic; ++j) {
            int max_dist = (fhia[i].frag_id == fhia[j].frag_id + 1 || fhia[i].frag_id == fhia[j].frag_id - 1) ? adj_dist : res_dist;
            if (abs(fhia[i].subject_offset - fhia[j].subject_offset) <= max_dist) cf.insert(pair<int, FragHapInfo*>(fhia[j].frag_id, &fhia[j]));            
        }
        for (int j = 0; j < fhic; ++j) {
            if (cf.find(fhia[j].frag_id + 1) != cf.end() && cf.find(fhia[j].frag_id - 1) != cf.end()) bf.insert(pair<int, FragHapInfo*>(fhia[j].frag_id, &fhia[j]));
        }
        if (cf.empty()) continue;

        subfrag_list.clear();
        subfrag_list.insert(cf.begin(), cf.end());
        subfrag_list.insert(bf.begin(), bf.end());
        if (fhia[i].hp == 1) {
            h1_subfrag_list.push_back(subfrag_list);
        } else {
            hbn_assert(fhia[i].hp == 2);
            h2_subfrag_list.push_back(subfrag_list);
        }
    }
    
    if (!h1_subfrag_list.empty()) s_merge_subfrag_list(h1_subfrag_list);
    if (!h2_subfrag_list.empty()) s_merge_subfrag_list(h2_subfrag_list);
    s_cut_subfrag_list(h1_subfrag_list, h2_subfrag_list);

    for (int i = 0; i < fhic; ++i) fhia[i].hp = 0;
    for (auto& v1 : h1_subfrag_list) {
        for (auto& x1 : v1) {
            x1.second->hp = 1;
        }
    }
    for (auto& v2 : h2_subfrag_list) {
        for (auto& x2 : v2) {
            x2.second->hp = 2;
        }
    }

    s_impute_one_list_large_dist(fhia, fhic);
}

static void*
haplotype_imputation_thread(void* params)
{
    FragHapImputationThreadWorkdata* data = (FragHapImputationThreadWorkdata*)(params);
    BamTagOptions* options = data->options();
    FragHapInfo* a;
    int c;
    while (data->get_next_read(&a, &c)) {
        pdqsort(a, a + c, [](const FragHapInfo& x, const FragHapInfo& y) {
            return (x.subject_id < y.subject_id) || (x.subject_id == y.subject_id && x.subject_offset < y.subject_offset);
        });
        int i = 0;
        while (i < c) {
            int j = i + 1;
            while (j < c && a[i].subject_id == a[j].subject_id) ++j;
            s_impute_one_list(a + i, j - i, options->impute_non_adj_frag_dist, options->impute_adj_frag_dist);
            i = j;
        }
    }

    return nullptr;
}

void haplotype_imputation_mt(BamTagOptions* options, FragHapInfo* fhia, size_t fhic)
{
    FragHapImputationThreadWorkdata data(options, fhia, fhic);
    pthread_t jobs[options->num_threads];
    for (int i = 0; i < options->num_threads; ++i) {
        pthread_create(jobs + i, nullptr, haplotype_imputation_thread, &data);
    }
    for (int i = 0; i < options->num_threads; ++i) {
        pthread_join(jobs[i], NULL);
    }
}