#include "lookup_table.hpp"
#include "parasort.h"
#include "pdqsort.h"

#include <algorithm>
#include <limits>
#include <mutex>
#include <pthread.h>

using namespace std;

void 
extract_hash_values_for_one_read(const int seq_id,
    const u8* seq,
    const u64 seq_size,
    const u64 seq_start_offset,
    const int kmer_size,
    const int kmer_window,
    const u64 min_hash,
    const u64 max_hash,
    std::vector<KmerHashAndOffset>& khao_list)
{
    if (seq_size < kmer_size) return;
    const int kIntersect = kmer_size > kmer_window;
    const int kStride = kmer_size - kmer_window;
    const u64 kIntersectMask = kIntersect ? ((U64_ONE << (kStride<<1)) - 1) : 0;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1); 
    KmerHashAndOffset khao;
    khao.hash = U64_MAX;
    khao.ki.seq_id = seq_id;
    khao.ki.seq_offset = -1;

    if (!kIntersect) {
        for (u64 j = 0; j <= seq_size - kmer_size; j += kmer_window) {
            u64 hash = 0;
            for (int k = 0; k < kmer_size; ++k) {
                const u64 pos = j + k;
				const u8 c = seq[pos];
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            khao.hash = hash;
            khao.ki.seq_offset = seq_start_offset + j;
            if (hash >= min_hash && hash < max_hash) khao_list.push_back(khao);
        }
    } else {
        u64 hash = 0;
        for (int j = 0; j < kmer_size; ++j) {
                const u64 pos = j;
				const u8 c = seq[pos];
                hash = (hash << 2) | c;
        }
        hbn_assert(hash < kMaxHashValue);
        khao.hash = hash;
        khao.ki.seq_offset = seq_start_offset;
        if (hash >= min_hash && hash < max_hash) khao_list.push_back(khao);
        for (u64 j = kmer_window; j <= seq_size - kmer_size; j += kmer_window) {
            hash &= kIntersectMask;
            for (int k = kStride; k < kmer_size; ++k) {
                const u64 pos = j + k;
				const u8 c = seq[pos];
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            khao.hash = hash;
            khao.ki.seq_offset = seq_start_offset + j;
            if (hash >= min_hash && hash < max_hash) khao_list.push_back(khao);
        }
    }
}

int
resolve_kmer_occ_cutoff(KmerHashAndOffset* khao_array, 
    u64 khao_array_size, 
    int num_threads,
    int kmer_size,
    double rep_frac, 
    int _max_kmer_occ)
{
    u64 dnk = 0;
    u64 i = 0;
    while (i < khao_array_size) {
        u64 j = i + 1;
        while (j < khao_array_size && khao_array[i].hash == khao_array[j].hash) ++j;
        ++dnk;
        i = j;
    }
    u64 cnt_list_size = dnk;
    int* cnt_list = new int[cnt_list_size];
    const u64 kKmerOccT = numeric_limits<int>::max() / 2;
    u64 cnt_list_idx = 0;
    i = 0;
    while (i < khao_array_size) {
        u64 j = i + 1;
        while (j < khao_array_size && khao_array[i].hash == khao_array[j].hash) ++j;
        u64 cnt = min<u64>(j - i, kKmerOccT);
        cnt_list[cnt_list_idx++] = cnt;
        i = j;
    }
    hbn_assert(cnt_list_idx == cnt_list_size);

    int* a = cnt_list;
    u64 c = cnt_list_size;
    //HBN_LOG("Sort cnt list");
    sort(a, a + c, greater<int>());
    //HBN_LOG("Done");
    const u64 kRepC = cnt_list_size * rep_frac;
    int cutoff = a[kRepC];
    HBN_LOG("Total kmers: %llu", khao_array_size);
    HBN_LOG("distinct kmers: %llu, rep_frac = %g, cutoff = %d", c, rep_frac, cutoff);
    fprintf(stderr, "max kmer occ: %llu, rep_c = %llu\n", a[0], kRepC);
	cutoff = hbn_max(cutoff, _max_kmer_occ);
    HBN_LOG("max kmer occ: %d", cutoff);

    delete[] cnt_list;
    return cutoff;
}

///////////////////////////////////

struct KmerStats
{
u64 total_kmers;
u64 distinct_kmers;
u64 removed_kmers;
u64 removed_distinct_kmers;

    KmerStats() {
        total_kmers = 0;
        distinct_kmers = 0;
        removed_kmers = 0;
        removed_distinct_kmers = 0;
    }

    void add(u64 tk, u64 dk, u64 rk, u64 rdk) {
        total_kmers += tk;
        distinct_kmers += dk;
        removed_kmers += rk;
        removed_distinct_kmers += rdk;
    }

    void stats() {
        double p = total_kmers ? 100.0 * removed_kmers / total_kmers : 0.0;
        fprintf(stderr, "Total kmers: %zu, %zu (%g%%) are removed\n", total_kmers, removed_kmers, p);
        p = distinct_kmers ? 100.0 * removed_distinct_kmers / distinct_kmers : 0.0;
        fprintf(stderr, "Distinct kmers: %zu, %zu (%g%%) are removed\n", distinct_kmers, removed_distinct_kmers, p);
    }
};

__HbnLookupTable::__HbnLookupTable(HbnUnpackedDatabase& db,
    const int kmer_size,
    const int kmer_window,
    const int max_kmer_occ,
    const double rep_frac,
    const int num_threads)
{
    fprintf(stderr, "kmer-size: %d, kmer-window: %d, kmer-occ-cutoff: %d\n", kmer_size, kmer_window, max_kmer_occ);
    vector<KmerHashAndOffset> khao_list;
    const int num_seq = db.NumSeqs();
    for (int i = 0; i < num_seq; ++i) {
        const u8* seq = db.GetSequence(i);
        const int seq_size = db.SeqSize(i);
        for (int k = 0; k < seq_size; ++k) {
            int c = seq[k];
            hbn_assert(c >= 0 && c < 4, "c = %d", c);
        }
        extract_hash_values_for_one_read(i, seq, seq_size, 0,
            kmer_size, kmer_window, 0, U64_MAX, khao_list);
    }
    HBN_LOG("Load %zu kmers", khao_list.size());

    parasort(khao_list.size(), khao_list.data(), num_threads);
    int kmer_occ_cutoff = resolve_kmer_occ_cutoff(khao_list.data(), khao_list.size(), num_threads, kmer_size, rep_frac, max_kmer_occ);

    M_hash_table = new HbnHashTable(kmer_size);
    M_offset_list = (KmerInfo*)calloc(khao_list.size(), sizeof(KmerInfo));
    M_offset_list_size = khao_list.size();
    u64 kmer_idx = 0;

    KmerStats stats;
    KmerHashAndOffset* a = khao_list.data();
    size_t ac = khao_list.size();
    size_t ii = 0;
    HBN_LOG("Add kmers into hash table");
    while (ii < ac) {
        size_t jj = ii + 1;
        while (jj < ac && a[ii].hash == a[jj].hash) ++jj;
        u64 cnt = jj - ii;
        if (cnt > kmer_occ_cutoff) { stats.add(cnt, 1, cnt, 1); ii = jj; continue; }
        stats.add(cnt, 1, 0, 0);
        u64 offset = kmer_idx;
        u64 u = lktbl_pack_hash_value(offset, cnt);
        M_hash_table->add_one_hash_and_value(a[ii].hash, u);
        for (size_t kk = ii; kk < jj; ++kk, ++kmer_idx) M_offset_list[kmer_idx] = a[kk].ki;
        ii = jj;
    }
    stats.stats();
    HBN_LOG("Done.");
    fprintf(stderr, "\n");
}
