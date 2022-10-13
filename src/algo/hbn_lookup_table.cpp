#include "hbn_lookup_table.h"

#include "../corelib/cstr_util.h"
#include "../corelib/sparse_map.h"
#include "hash_list_bucket_sort.h"

#include <cstdlib>
#include <pthread.h>
#include <vector>

using namespace std;

#define LKTBL_OFFSET_BIT    48
#define LKTBL_CNT_BIT   16
#define LKTBL_CNT_MASK ((U64_ONE<<LKTBL_CNT_BIT)-1)

#define lktbl_pack_hash_value(offset_, cnt_) (((offset_) << LKTBL_CNT_BIT) | (cnt_))
#define lktbl_extract_offset(u_) ((u_) >> LKTBL_CNT_BIT)
#define lktbl_extract_cnt(u_) ((u_) & LKTBL_CNT_MASK)

#define USE_LARGE_TABLE 1

#if USE_LARGE_TABLE
struct BobHenkin_u64tou64
{
    u64 operator()(uint64_t key) const {
            key = (~key) + (key << 21); // key = (key << 21) - key - 1; 
            key = key ^ (key >> 24); 
            key = (key + (key << 3)) + (key << 8); // key * 265 
            key = key ^ (key >> 14); 
            key = (key + (key << 2)) + (key << 4); // key * 21 
            key = key ^ (key >> 28); 
            key = key + (key << 31); 
            return key;
    }
};

typedef tsl::sparse_map<u64, u64, BobHenkin_u64tou64> hash2offset_map_t;
#else
static u64*
s_AllocHash2kmiMap(const int kmer_size)
{
    hbn_assert(kmer_size <= 16);
    const u64 item = U64_ONE << (kmer_size * 2);
    u64* hash2kmi_map = (u64*)calloc(item, sizeof(u64));
    return hash2kmi_map;
}
#endif

typedef struct{
    u64 hash;
    u64 offset;
} KmerHashAndOffset;

static void 
s_build_seq_offset_list(u64** p_seq_offset_list, u64* p_seq_offset_list_size, const int qidf, const int qidt, SeqReader* db)
{
    u64* seq_offset_list = *p_seq_offset_list;
    if (seq_offset_list) free(seq_offset_list);
    int nq = SeqReader_NumSeqs(db);
    seq_offset_list = (u64*)malloc(sizeof(u64) * (nq+1));
    u64 offset = 0;  
    for (int i = 0; i < nq; ++i) {
        seq_offset_list[i] = offset;
        offset += SeqReader_SeqSize(db, i);
    }
    seq_offset_list[nq] = offset;
    int seq_offset_list_size = nq;

    *p_seq_offset_list = seq_offset_list;
    *p_seq_offset_list_size = seq_offset_list_size;  
}

static void
s_extract_hash_values(const u64 seq_offset, 
    const u64 subseq_offset, 
    const u8* subseq, 
    const int subseq_size,
    const int kmer_size,
    const int kmer_window,
    KmerHashAndOffset* khao_list,
    u64* p_khao_idx)
{
    const int kIntersect = kmer_size > kmer_window;
    const int kStride = kmer_size - kmer_window;
    const u64 kIntersectMask = kIntersect ? ((U64_ONE << (kStride<<1)) - 1) : 0;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1);    
    KmerHashAndOffset khao;
    u64 khao_idx = *p_khao_idx;

    if (!kIntersect) {
        for (u64 j = 0; j <= subseq_size - kmer_size; j += kmer_window) {
            u64 hash = 0;
            for (int k = 0; k < kmer_size; ++k) {
                const u64 pos = j + k;
				const u8 c = subseq[pos];
                hbn_assert(c < 4);
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            khao.hash = hash;
            khao.offset = seq_offset + subseq_offset + j;
            khao_list[khao_idx++] = khao;
        }
    } else {
        u64 hash = 0;
        for (int j = 0; j < kmer_size; ++j) {
            const u64 pos = j;
			const u8 c = subseq[pos];
            hbn_assert(c < 4, "j = %d, c = %d, pos = %zu", j, c, pos);
            hash = (hash << 2) | c;
        }
        hbn_assert(hash < kMaxHashValue);
        khao.hash = hash;
        khao.offset = seq_offset + subseq_offset;
        khao_list[khao_idx++] = khao;
        for (u64 j = kmer_window; j <= subseq_size - kmer_size; j += kmer_window) {
            hash &= kIntersectMask;
            for (int k = kStride; k < kmer_size; ++k) {
                const u64 pos = j + k;
				const u8 c = subseq[pos];
                hbn_assert(c < 4);
                hash = (hash << 2) | c;
            }
            hbn_assert(hash < kMaxHashValue);
            khao.hash = hash;
            khao.offset = seq_offset + subseq_offset + j;
            khao_list[khao_idx++] = khao;
        }
    }

    *p_khao_idx = khao_idx;
}

static int 
s_calc_khao_cnt_for_one_seq(const int* enzyme_list, 
    const int enzyme_list_size,
    const int seq_id,
    const int seq_size, 
    const int enzyme_size,
    const int kmer_size,
    const int kmer_window)
{
    if (seq_size < kmer_size) return 0;
    return (seq_size - kmer_size) / kmer_window + 1;
}

static u64 
hash_extractor(void* list, const u64 i)
{
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)(list);
    return khao_array[i].hash;
}

static u64 
offset_extractor(void* list, const u64 i)
{
    KmerHashAndOffset* khao_array = (KmerHashAndOffset*)(list);
    return khao_array[i].offset;
}

static void
set_khao_array_item_value(void* src, const u64 src_idx, void* dst, const u64 dst_idx)
{
    KmerHashAndOffset* src_khao_array = (KmerHashAndOffset*)(src);
    KmerHashAndOffset* dst_khao_array = (KmerHashAndOffset*)(dst);
    dst_khao_array[dst_idx] = src_khao_array[src_idx];
}

static void
s_fill_khao_list(SeqReader* db,
    RestrictEnzymeLociList* relist,
    const int kmer_size, 
    const int kmer_window,
    const u64* seq_offset_list,
    KmerHashAndOffset* khao_list,
    const u64 khao_list_size)
{
    const int num_seqs = SeqReader_NumSeqs(db);
    const int enzyme_size = relist->enzyme.enzyme_size;
    u64 khao_idx = 0;
    for (int i = 0; i < num_seqs; ++i) {
        int seq_size = SeqReader_SeqSize(db, i);
        if (seq_size < kmer_size) continue;
        const u8* seq = SeqReader_Seq(db, i, FWD);
        s_extract_hash_values(seq_offset_list[i], 0, seq, seq_size, kmer_size, kmer_window, khao_list, &khao_idx);
    }
    hbn_assert(khao_idx == khao_list_size, "idx = %zu, num = %zu", khao_idx, khao_list_size);
}

static int
s_resolve_max_kmer_occ(const KmerHashAndOffset* khao_array, const size_t khao_array_size, double rep_frac, int max_kmer_occ)
{
    kv_dinit(vec_u64, cnt_list);
    u64 distinct_kmers = 0;
    size_t i = 0;
    while (i < khao_array_size) {
        size_t j = i + 1;
        while (j < khao_array_size && khao_array[i].hash == khao_array[j].hash) ++j;
        ++distinct_kmers;
        size_t cnt = j - i;
        kv_push(u64, cnt_list, cnt);
        i = j;
    }

    u64* a = kv_data(cnt_list);
    u64 c = kv_size(cnt_list);
    ks_introsort_u64(c, a);
    u64 x = 0, y = c;
    while (x < y) {
        u64 t = a[x];
        a[x] = a[y-1];
        a[y-1] = t;
        ++x;
        --y;
    }
    double p = 100.0 * a[0] / khao_array_size;
    HBN_LOG("Max word count: %zu (%g%%)", a[0], p);
    const u64 kRepC = distinct_kmers * rep_frac;
    int cutoff = a[kRepC];
    cutoff = hbn_max(cutoff, max_kmer_occ);
    HBN_LOG("Filter words occurring more than %d times", cutoff);
    kv_destroy(cnt_list);
    return cutoff;
}

static void
s_construct_lookup_table(SeqReader* db,
    RestrictEnzymeLociList* relist,
    const int kmer_size, 
    const int kmer_window,
    double rep_frac,
    const int _max_kmer_occ, 
    const int num_threads,
    const u64* seq_offset_list,
    u64** p_kmer_offset_list,
    u64* p_kmer_offset_list_size,
#if USE_LARGE_TABLE
    hash2offset_map_t** p_hash2offset_map
#else
    u64** p_hash2offset_map
#endif
)
{
    const int num_seqs = SeqReader_NumSeqs(db);
    u64 total_kmers = 0;
    const int enzyme_size = relist->enzyme.enzyme_size;
    for (int i = 0; i < num_seqs; ++i) {
        int seq_size = SeqReader_SeqSize(db, i);
        if (seq_size < kmer_size) continue;
        const int* enzyme_list = relist->reloci_array + relist->seq_reloci_info_array[i].enzyme_loci_offset;
        const int enzyme_list_size = relist->seq_reloci_info_array[i].enzyme_loci_cnt - 1;
        hbn_assert(enzyme_list_size > 1);
        int num_kmers = s_calc_khao_cnt_for_one_seq(enzyme_list, enzyme_list_size, i, seq_size, enzyme_size, kmer_size, kmer_window);
        total_kmers += num_kmers;
    }

    KmerHashAndOffset* khao_list = (KmerHashAndOffset*)calloc(total_kmers, sizeof(KmerHashAndOffset));
    s_fill_khao_list(db, relist, kmer_size, kmer_window, seq_offset_list, khao_list, total_kmers);

    HBN_LOG("Word size: %d", kmer_size);
    HBN_LOG("Word sampling window: %d", kmer_window);
    HBN_LOG("Load %zu words", total_kmers);
    HBN_LOG("Sort word list by hash values with %d CPU threads", num_threads);
    radix_sort(khao_list,
        sizeof(KmerHashAndOffset), 
        total_kmers, 
        num_threads,
        offset_extractor,
        hash_extractor,
        set_khao_array_item_value);
    HBN_LOG("Done.");
    for (u64 i = 0; i < total_kmers - 1; ++i) {
        hbn_assert(khao_list[i].hash <= khao_list[i+1].hash, "i = %zu, %zu --- %zu", i, khao_list[i].hash, khao_list[i+1].hash);
    }

    u64 distinct_kmers = 0;
    u64 removed_distinct_kmers = 0;
    u64 removed_kemrs = 0;
    u64 i = 0;
    int max_kmer_occ = s_resolve_max_kmer_occ(khao_list, total_kmers, rep_frac, _max_kmer_occ);
    while (i < total_kmers) {
        u64 j = i + 1;
        while (j < total_kmers && khao_list[i].hash == khao_list[j].hash) ++j;
        u64 x = j - i;
        if (x > max_kmer_occ) {
            int ruk = 0, rk = 0;
            for (u64 k = i; k < j; ++k) {
                khao_list[k].hash = U64_MAX;
                ruk = 1;
                rk++;
            }
            removed_distinct_kmers += ruk;
            removed_kemrs += rk;
        }
        distinct_kmers++;
        i = j;
    }

    char buf1[64], buf2[64], buf3[64];
    u64_to_string_comma(total_kmers, buf1);
    u64_to_string_comma(removed_kemrs, buf2);
    double perc = 100.0 * removed_kemrs / total_kmers;
    double_to_string(perc, buf3);
    HBN_LOG("Total words: %s", buf1);
    HBN_LOG("Filtered words: %s (%s%%)", buf2, buf3);
    u64_to_string_comma(distinct_kmers, buf1);
    u64_to_string_comma(removed_distinct_kmers, buf2);
    perc = 100.0 * removed_distinct_kmers / distinct_kmers;
    double_to_string(perc, buf3);
    HBN_LOG("Distinct words: %s", buf1);
    HBN_LOG("Filtered distinct words: %s (%s%%)", buf2, buf3);

    const u64 kmer_offset_list_size = total_kmers - removed_kemrs;
    u64* kmer_offset_list = (u64*)calloc(kmer_offset_list_size, sizeof(u64));
    u64* hash_list = (u64*)calloc(kmer_offset_list_size, sizeof(u64));
    u64 kmer_idx = 0;
    for (i = 0; i < total_kmers; ++i) {
        if (khao_list[i].hash == U64_MAX) continue;
        kmer_offset_list[kmer_idx] = khao_list[i].offset;
        hash_list[kmer_idx] = khao_list[i].hash;
        ++kmer_idx;
    }
    hbn_assert(kmer_idx == kmer_offset_list_size);
    free(khao_list);
    khao_list = NULL;
    total_kmers = 0;

    HBN_LOG("Build hash table");
    i = 0;
    kmer_idx = 0;
#if USE_LARGE_TABLE
    hash2offset_map_t* hash2offset_map = new hash2offset_map_t();
    hash2offset_map->reserve(kmer_offset_list_size);
#else 
    u64* hash2offset_map = s_AllocHash2kmiMap(kmer_size);
#endif
    while (i < kmer_offset_list_size) {
        u64 j = i + 1;
        while (j < kmer_offset_list_size && hash_list[j] == hash_list[i]) ++j;
        u64 x = j - i;
        u64 u = lktbl_pack_hash_value(kmer_idx, x);
#if USE_LARGE_TABLE
        hash2offset_map->insert(pair<u64, u64>(hash_list[i], u));
#else
        hash2offset_map[hash_list[i]] = u;
#endif
        kmer_idx += x;
        i = j;
    }
    HBN_LOG("Done\n");
    free(hash_list);
    *p_kmer_offset_list = kmer_offset_list;
    *p_kmer_offset_list_size = kmer_offset_list_size;
    *p_hash2offset_map = hash2offset_map;
}

////////////////////

void HbnLookupTable_KmerOffset2SeqInfo(HbnLookupTable* lktbl, u64 offset, int* id, size_t* start, size_t* end)
{
    int nq = lktbl->m_seq_offset_list_size;
    const u64* m_seq_offset_list = lktbl->m_seq_offset_list;
    int left = 0, mid = 0, right = nq;
    while (left < right) {
        mid = (left + right) >> 1;
        if (offset >= m_seq_offset_list[mid]) {
            if (mid == nq - 1) break;
            if (offset < m_seq_offset_list[mid+1]) break;
            left = mid + 1;
        } else {
            right = mid;
        }
    }
    hbn_assert(mid >= 0);
    hbn_assert(mid < nq);
    *id = mid;
    *start = m_seq_offset_list[mid];
    *end = m_seq_offset_list[mid+1]; 
}

u64* 
HbnLookupTable_ExtractKmerOffsetList(HbnLookupTable* lktbl, u64 hash, u64* cnt)
{
#if USE_LARGE_TABLE
    hash2offset_map_t* m_hash2offset_map = static_cast<hash2offset_map_t*>(lktbl->m_hash2offset_map);
    u64* m_kmer_offset_list = lktbl->m_kmer_offset_list;
    *cnt = 0;
    auto pos = m_hash2offset_map->find(hash);
    if (pos == m_hash2offset_map->end()) return NULL;

    u64 u = pos->second;
    u64 offset = lktbl_extract_offset(u);
    *cnt = lktbl_extract_cnt(u);
    return m_kmer_offset_list + offset;
#else
    *cnt = 0;
    const u64* m_hash2offset_map = (u64*)lktbl->m_hash2offset_map;
    u64* m_kmer_offset_list = lktbl->m_kmer_offset_list;
    const u64 u = m_hash2offset_map[hash];
    if (!u) return NULL;

    u64 offset = lktbl_extract_offset(u);
    *cnt = lktbl_extract_cnt(u);
    return lktbl->m_kmer_offset_list + offset;
#endif
}

const u64*
HbnLookupTable_SeqOffsetList(HbnLookupTable* lktbl)
{
    return lktbl->m_seq_offset_list;
}

HbnLookupTable*
HbnLookupTableFree(HbnLookupTable* lktbl)
{
    if (!lktbl) return NULL;
    u64* m_kmer_offset_list = lktbl->m_kmer_offset_list;
    u64* m_seq_offset_list = lktbl->m_seq_offset_list;
#if USE_LARGE_TABLE
    hash2offset_map_t* m_hash2offset_map = static_cast<hash2offset_map_t*>(lktbl->m_hash2offset_map);
    if (m_hash2offset_map) delete m_hash2offset_map;
#else
    u64* m_hash2offset_map = (u64*)lktbl->m_hash2offset_map;
    if (m_hash2offset_map) free(m_hash2offset_map);
#endif 

    if (m_kmer_offset_list) free(m_kmer_offset_list);
    if (m_seq_offset_list) free(m_seq_offset_list);
    free(lktbl);
    return NULL;
}

HbnLookupTable*
HbnLookupTableNew(SeqReader* db, 
    RestrictEnzymeLociList* relist,
    int kmer_size, 
    int kmer_window,
    double rep_frac,
    int max_kmer_occ, 
    int num_threads)
{
    HbnLookupTable* lktbl = (HbnLookupTable*)calloc(1, sizeof(HbnLookupTable));
    s_build_seq_offset_list(&lktbl->m_seq_offset_list, &lktbl->m_seq_offset_list_size, 0, SeqReader_NumSeqs(db), db);
#if USE_LARGE_TABLE
    hash2offset_map_t* hash2offset_map = NULL;
#else
    u64* hash2offset_map = NULL;
#endif 
    s_construct_lookup_table(db, relist, kmer_size, kmer_window, rep_frac, max_kmer_occ, num_threads, 
        lktbl->m_seq_offset_list, &lktbl->m_kmer_offset_list, 
        &lktbl->m_kmer_offset_list_size, &hash2offset_map);
    lktbl->m_hash2offset_map = hash2offset_map;
    return lktbl;
}
