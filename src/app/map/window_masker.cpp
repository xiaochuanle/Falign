#include "window_masker.hpp"

#include "../../corelib/parasort.h"
#include "../../corelib/pdqsort.h"
#include "../dip3d/split_string_by_char.hpp"

#include <algorithm>
#include <limits>
#include <set>

using namespace std;

static int
s_resolve_kmer_size(HbnUnpackedDatabase* reference)
{
    size_t dbsize = 0;
    const int num_chr = reference->NumSeqs();
    for (int i = 0; i < num_chr; ++i) dbsize += reference->SeqSize(i);
    HBN_LOG("DB size: %zu", dbsize);
    int kmer = 1;
    int x = 4;
    const double kFactor = 5.0;
    while (1) {
        double f = 1.0 * dbsize / x;
        if (f < kFactor) break;
        ++kmer;
        x *= 4;
    }
    fprintf(stderr, "Chroms: %d, bases: %zu, kmer: %d\n", num_chr, dbsize, kmer);
    return kmer;
}

static void
s_compute_hash_counts_for_one_chr(int* hash_counts,
    const u8* chr, 
    const int chr_size,
    const int kmer)
{
    if (chr_size < kmer) return;
    const u64 kMaxHashValue = U64_ONE << (kmer << 1); 
    const u64 kHashMask = kMaxHashValue - 1;

    u64 fhash = 0, rhash = 0;
    for (int i = 0; i < kmer - 1; ++i) {
        u8 c = chr[i];
        fhash = (fhash << 2) | c;
        fhash &= kHashMask;
        
        c = 3 - chr[chr_size-1-i];
        rhash = (rhash << 2) | c;
        rhash &= kHashMask;
    }

    for (int p = kmer - 1; p < chr_size; ++p) {
        u8 c = chr[p];
        fhash = (fhash << 2) | c;
        fhash &= kHashMask;
        ++hash_counts[fhash];

        c = 3 - chr[chr_size-1-p];
        rhash = (rhash << 2) | c;
        rhash &= kHashMask;
        ++hash_counts[rhash];
    }
}

static void
s_compute_kmer_score_and_threshold_value(HbnUnpackedDatabase* reference, 
    int kmer_size, 
    int* kmer_scores,
    int& T_threshold,
    int& T_extend,
    int& T_low,
    int& T_high)
{
    const int num_chr = reference->NumSeqs();
    for (int i = 0; i < num_chr; ++i) {
        const char* name = reference->SeqName(i);
        const u8* chr = reference->GetSequence(i);
        const int chr_size = reference->SeqSize(i);
        //HBN_LOG("%d:%s:%d", i, name, chr_size);
        s_compute_hash_counts_for_one_chr(kmer_scores, chr, chr_size, kmer_size);
    }

    vector<int> occ_list;
    const u64 kMaxHash = U64_ONE << (2 * kmer_size);
    for (u64 i = 0; i < kMaxHash; ++i) if (kmer_scores[i]) occ_list.push_back(kmer_scores[i]);

    sort(occ_list.begin(), occ_list.end(), greater<int>());
    size_t NK = occ_list.size();
    HBN_LOG("Distinct kmers: %zu, max-occ: %d, min-occ: %d", NK, occ_list[0], occ_list.back());

    size_t N = NK * (1.0 - 0.995);
    T_threshold = occ_list[N];
    N = NK * (1.0 - 0.99);
    T_extend = occ_list[N];
    N = NK * (1.0 - 0.9);
    T_low = occ_list[N];
    N = NK * (1.0 - 0.998);
    T_high = occ_list[N];

    for (u64 i = 0; i < kMaxHash; ++i) {
        if (!kmer_scores[i]) continue;
        int s = kmer_scores[i];
        if (s >= T_high) {
            kmer_scores[i] = T_high;
        } else if (s < T_low) {
            kmer_scores[i] = T_low / 2;
        } else {
            kmer_scores[i] = s;
        }
    }
}

int*
s_compute_window_score_for_one_chr(const int* kmer_scores, const u8* chr, const int chr_size, const int kmer, const int win_ext)
{
    if (chr_size < kmer) return nullptr;
    const u64 kMaxHashValue = U64_ONE << (kmer << 1); 
    const u64 kHashMask = kMaxHashValue - 1;

    int* scores = new int[chr_size];
    fill(scores, scores + chr_size, 0);

    u64 hash = 0;
    for (int i = 0; i < kmer - 1; ++i) {
        u8 c = chr[i];
        hash = (hash << 2) | c;
        hash &= kHashMask;
    }

    for (int i = 0, p = kmer - 1; p < chr_size; ++i, ++p) {
        u8 c = chr[p];
        hash = (hash << 2) | c;
        hash &= kHashMask;
        scores[i] = kmer_scores[hash];
    }

    int win_size = kmer + win_ext;
    for (int i = 0; i < chr_size - win_ext; ++i) {
        int score = 0;
        for (int j = 0; j <= win_ext; ++j) score += scores[i+j];
        score /= (win_ext + 1);
        scores[i] = score;
    }

    return scores;
}

static void
s_mask_one_chr(HbnUnpackedDatabase* reference,
    const int chr_id,
    const int kmer,
    const int* kmer_scores,
    const int T_threshold,
    const int T_extend,
    const int T_low,
    const int T_high,
    const int min_masked_intv_size,
    vector<ChrIntv>& masked_intv_list)
{
    const char* chr_name = reference->SeqName(chr_id);
    const u8* chr = reference->GetSequence(chr_id);
    const int chr_size = reference->SeqSize(chr_id);
    //fprintf(stderr, "%s:%d\n", chr_name, chr_size);
    u8* mask_list = new u8[chr_size];
    fill(mask_list, mask_list + chr_size, 0);
    const int kWinExt = 4;
    const int kWinSize = kmer + kWinExt;

    int* window_scores = s_compute_window_score_for_one_chr(kmer_scores, chr, chr_size, kmer, kWinExt);
    int* x_scores = new int[chr_size];
    fill(x_scores, x_scores + chr_size, 0);
    for (int i = 0; i < chr_size - kWinSize; ++i) {
        if (window_scores[i] >= T_threshold) {
            for (int j = 0; j < kWinSize; ++j) mask_list[i+j] = 1;
        }
        for (int j = 0; j < kWinSize; ++j) x_scores[i+j] = max(x_scores[i+j], window_scores[i]);
    }

    int x = 0;
    while (x < chr_size) {
        while (x < chr_size && mask_list[x]) ++x;
        if (x >= chr_size) break;
        int y = x + 1;
        while (y < chr_size && mask_list[y] == 0) ++y;
        bool mask_this_intv = true;
        for (int z = x; z < y; ++z) {
            if (x_scores[z] < T_extend) {
                mask_this_intv = false;
                break;
            }
        }
        if (mask_this_intv) for (int z = x; z < y; ++z) mask_list[z] = 1;
        x = y;
    }

    x = 0;
    ChrIntv intv;
    while (x < chr_size) {
        while (x < chr_size && mask_list[x] == 0) ++x;
        if (x >= chr_size) break;
        int y = x + 1;
        while (y < chr_size && mask_list[y]) ++y;
        if (y - x >= min_masked_intv_size) {
            intv.chr_id = chr_id;
            intv.from = x;
            intv.to = y;
            masked_intv_list.push_back(intv);
        }
        x = y;
    }

    delete[] window_scores;
    delete[] x_scores;
    delete[] mask_list;
}

static void 
s_merge_adj_intvs(vector<ChrIntv>& intv_list)
{
    pdqsort(intv_list.begin(), intv_list.end(), [](const ChrIntv& x, const ChrIntv& y) {
        return (x.chr_id < y.chr_id) || (x.chr_id == y.chr_id && x.from < y.from);
    });

    int N = intv_list.size();
    int i = 0;
    while (i < N) {
        int j = i + 1;
        int last_to = intv_list[i].to;
        while (j < N && intv_list[i].chr_id == intv_list[j].chr_id) {
            int from = intv_list[j].from;
            int to = intv_list[j].to;
            if (from - last_to > 10) break;
            last_to = max(last_to, to);
            ++j;
        }
        if (j - i == 1) { i = j; continue; }

        for (int k = i + 1; k < j; ++k) intv_list[k].chr_id = -1;
        intv_list[i].to = last_to;
        i = j;
    }
    int n = 0;
    for (i = 0; i < N; ++i) if (intv_list[i].chr_id >= 0) intv_list[n++] = intv_list[i];
    intv_list.resize(n);
}

void compute_repeat_regions(HbnUnpackedDatabase* reference,
    const int min_masked_intv_size,
    std::vector<ChrIntv>& repeat_intv_list)
{
    int kmer = s_resolve_kmer_size(reference);
    int T_threshold = 0, T_extend = 0, T_low = 0, T_high = 0;
    const u64 kMaxHashValue = U64_ONE << (2 * kmer);
    int* kmer_scores = new int[kMaxHashValue];
    fill(kmer_scores, kmer_scores + kMaxHashValue, 0);
    s_compute_kmer_score_and_threshold_value(reference, kmer, kmer_scores, T_threshold, T_extend, T_low, T_high);
    fprintf(stderr, "T_threshold: %d, T_extend: %d, T_low: %d, T_high: %d\n", T_threshold, T_extend, T_low, T_high);

    const int num_chr = reference->NumSeqs();
    for (int i = 0; i < num_chr; ++i) s_mask_one_chr(reference, i, kmer, kmer_scores, T_threshold, T_extend, T_low, T_high, min_masked_intv_size, repeat_intv_list);
    
    //s_merge_adj_intvs(repeat_intv_list);

    delete[] kmer_scores;
}

//////////////////

void load_bed_list(const char* bed_path, SeqName2IdMap& refname2id, std::vector<ChrIntv>& intv_list)
{
    ChrIntv intv;
    vector<pair<const char*, int>> cols;
    HbnLineReader in(bed_path);
    string last_chr_name;
    int last_chr_id = -1;
    while (in.ReadOneLine()) {
        NStr::CTempString line = *in;
        cols.clear();
        split_string_by_char(line.data(), line.size(), '\t', cols);
        if (cols.size() < 3) HBN_ERR("Corrupted file %s", bed_path);
        int n1 = last_chr_name.size();
        int n2 = cols[0].second;
        int n0 = min(n1, n2);
        if (n1 != n2 || strncmp(last_chr_name.c_str(), cols[0].first, n0)) {
            last_chr_name.assign(cols[0].first, cols[0].second);
            last_chr_id = refname2id.GetIdFromNameSafe(last_chr_name);
        }
        intv.chr_id = last_chr_id;
        intv.from = atoi(cols[1].first);
        intv.to = atoi(cols[2].first);
        intv_list.push_back(intv);
    }

    pdqsort(intv_list.begin(), intv_list.end(), [](const ChrIntv& x, const ChrIntv& y) {
        return (x.chr_id < y.chr_id) || (x.chr_id == y.chr_id && x.from < y.from);
    });

    int N = intv_list.size();
    int i = 0;
    while (i < N) {
        int j = i + 1;
        int last_to = intv_list[i].to;
        while (j < N && intv_list[i].chr_id == intv_list[j].chr_id) {
            int from = intv_list[j].from;
            int to = intv_list[j].to;
            if (from - last_to > 0) break;
            last_to = max(last_to, to);
            ++j;
        }
        if (j - i == 1) { i = j; continue; }

        for (int k = i + 1; k < j; ++k) intv_list[k].chr_id = -1;
        intv_list[i].to = last_to;
        i = j;
    }
    int n = 0;
    for (i = 0; i < N; ++i) if (intv_list[i].chr_id >= 0) intv_list[n++] = intv_list[i];
    intv_list.resize(n);

    size_t num_bases = 0;
    int min_intv = numeric_limits<int>::max();
    int max_intv = numeric_limits<int>::min();
    for (i = 0; i < n; ++i) {
        int L = (intv_list[i].to - intv_list[i].from);
        num_bases += L;
        min_intv = min(min_intv, L);
        max_intv = max(max_intv, L);
    }
    string size = NStr::UInt8ToString_DataSize(num_bases);
    HBN_LOG("Load %zu bed segments (%s) from %s", intv_list.size(), size.c_str(), bed_path);
    HBN_LOG("min-intv: %d, max-intv: %d", min_intv, max_intv);
}

static void
s_compute_complement_intv_list(ChrIntv* intva, int intvc, int chr_size, vector<ChrIntv>& complement_intv_list)
{
    ChrIntv intv;
    intv.chr_id = intva[0].chr_id;
    if (intva[0].from) {
        intv.from = 0;
        intv.to = intva[0].from;
        complement_intv_list.push_back(intv);
    }
    for (int i = 0; i < intvc - 1; ++i) {
        intv.from = intva[i].to;
        intv.to = intva[i+1].from;
        hbn_assert(intv.from < intv.to);
        complement_intv_list.push_back(intv);
    }
    if (intva[intvc-1].to < chr_size) {
        intv.from = intva[intvc-1].to;
        intv.to = chr_size;
        complement_intv_list.push_back(intv);
    }
}

void compute_complement_intv_list(ChrIntv* a, size_t c, HbnUnpackedDatabase& reference, std::vector<ChrIntv>& complement_intv_list)
{
    set<int> added_chrs;
    size_t i = 0;
    while (i < c) {
        int chr_id = a[i].chr_id;
        int chr_size = reference.SeqSize(chr_id);
        size_t j = i + 1;
        while (j < c && a[i].chr_id == a[j].chr_id) ++j;
        s_compute_complement_intv_list(a + i, j - i, chr_size, complement_intv_list);
        added_chrs.insert(chr_id);
        i = j;
    }

    const int num_chr = reference.NumSeqs();
    for (int p = 0; p < num_chr; ++p) {
        if (added_chrs.find(p) != added_chrs.end()) continue;
        int chr_size = reference.SeqSize(p);
        ChrIntv intv;
        intv.chr_id = p;
        intv.from = 0;
        intv.to = chr_size;
        complement_intv_list.push_back(intv);
    }

    int min_intv = numeric_limits<int>::max();
    int max_intv = numeric_limits<int>::min();
    size_t num_bases = 0;
    for (auto& intv : complement_intv_list) {
        int L = intv.to - intv.from;
        min_intv = min(min_intv, L);
        max_intv = max(max_intv, L);
        num_bases += L;
    }
    string size = NStr::UInt8ToString_DataSize(num_bases);
    HBN_LOG("Complement bed segments %zu (%s)", complement_intv_list.size(), size.c_str());
    HBN_LOG("min-intv: %d, max-intv: %d", min_intv, max_intv);
}

void repeat_intv_stats(ChrIntv* a, size_t c, const HbnUnpackedDatabase* reference)
{
    size_t masked_bases = 0;
    int max_intv_size = 0;
    ChrIntv* max_intv = nullptr;
    for (size_t i = 0; i < c; ++i) {
        int L = a[i].to - a[i].from;
        masked_bases += L;
        if (L > max_intv_size) {
            max_intv_size = L;
            max_intv = a + i;
        }
    }
    string size = NStr::UInt8ToString_DataSize(masked_bases);
    fprintf(stderr, "Bed region stats:\n");
    fprintf(stderr, "    Bed intvs: %zu\n", c);
    fprintf(stderr, "    Bed bases: %s\n", size.c_str());
    if (reference) {
        const char* chr_name = reference->SeqName(max_intv->chr_id);
        fprintf(stderr, "    Max-intv: [%s, %d, %d, %d]\n", chr_name, max_intv->from, max_intv->to, max_intv->to - max_intv->from);
    } else {
        fprintf(stderr, "    Max-intv: [%d, %d, %d, %d]\n", max_intv->chr_id, max_intv->from, max_intv->to, max_intv->to - max_intv->from);
    }
}
