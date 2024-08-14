#ifndef __LOOKUP_TABLE_HPP
#define __LOOKUP_TABLE_HPP

#include "hbn_aux.h"
#include "unpacked_seqdb.hpp"
#include "sparse_map.h"

#include <vector>

#define DUMP_REPEAT_BED_KMER_STATS 0

#define LKTBL_OFFSET_BIT    48
#define LKTBL_CNT_BIT   16

#define HASH_BUCKETS 16
#define HASH_BUCKET_BITS 4 // 2^4 = 16

inline u64 lktbl_pack_hash_value(u64 offset, u64 cnt)
{
    return (offset << LKTBL_CNT_BIT) | cnt;
}

inline u64 lktbl_extract_offset(u64 u) 
{
    return u >> LKTBL_CNT_BIT;
}

inline u64 lktbl_extract_cnt(u64 u) 
{
    u64 mask = 1;
    mask = (mask << LKTBL_CNT_BIT) - 1;
    return u & mask;
}

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

using large_hash_type = tsl::sparse_map<u64, u64, BobHenkin_u64tou64>;

struct ChrIntv
{
    int chr_id;
    int from;
    int to;
};

struct KmerInfo 
{
    int seq_id;
    int seq_offset;
};

struct KmerHashAndOffset
{
    u64 hash;
    KmerInfo ki;
#if DUMP_REPEAT_BED_KMER_STATS
    u64 bed_id;
#endif

    bool operator < (const KmerHashAndOffset& rhs) const {
        return hash < rhs.hash;
    }

    bool operator == (const KmerHashAndOffset& rhs) const {
        return hash == rhs.hash;
    }
};

struct KmerMatch
{
    u8 qdir;
    u8 is_repeat_kmer;
    int qoff;
    int sid;
    int soff;

    bool operator < ( const KmerMatch& rhs) const {
        return (this->sid < rhs.sid) || (this->sid == rhs.sid && this->soff < rhs.soff);
    }
};

void 
extract_hash_values_for_one_read(const int seq_id,
    const u8* seq,
    const u64 seq_size,
    const u64 seq_start_offset,
    const int kmer_size,
    const int kmer_window,
    const u64 min_hash,
    const u64 max_hash,
    std::vector<KmerHashAndOffset>& khao_list);

int
resolve_kmer_occ_cutoff(KmerHashAndOffset* khao_array, 
    u64 khao_array_size, 
    double rep_frac, 
    int _max_kmer_occ);

/////////////////////////

class HbnHashTable
{
public:
    HbnHashTable(const int kmer_size) {
        hbn_assert(kmer_size > 0 && kmer_size <= kMaxKmerSize, "kmer = %d", kmer_size);
        M_max_hash_value = U64_ONE << (2 * kmer_size);
        if (kmer_size <= kLargeTableKmerCutoff) {
            HBN_LOG("Use small hash table for k = %d", kmer_size);
            M_small_hash = (u64*)calloc(M_max_hash_value, sizeof(u64));
            M_large_hash = nullptr;   
        } else {
            HBN_LOG("Use large hash table for k = %d", kmer_size);
            M_small_hash = nullptr;
            M_large_hash = new large_hash_type;
        }
    }

    ~HbnHashTable() {
        if (M_small_hash) free(M_small_hash);
        if (M_large_hash) delete M_large_hash;
    }

    u64 extract_value(const u64 hash) const {
        hbn_assert(hash < M_max_hash_value);
        if (M_small_hash) {
            return M_small_hash[hash];
        } else {
            auto pos = M_large_hash->find(hash);
            return (pos == M_large_hash->end()) ? 0 : pos->second;
        }
    }

    void add_one_hash_and_value(const u64 hash, const u64 value) {
        hbn_assert(hash < M_max_hash_value);
        if (M_small_hash) {
            M_small_hash[hash] = value;
        } else {
            M_large_hash->insert(std::pair<u64, u64>(hash, value));
        }
    }

    void reserve(const u64 num_hash) {
        if (M_large_hash) M_large_hash->reserve(num_hash);
    }

private:
    static const int kLargeTableKmerCutoff = 16;
    static const int kMaxKmerSize = 31;

    u64                 M_max_hash_value;
    u64*                M_small_hash;
    large_hash_type*    M_large_hash;
};

class __HbnLookupTable
{
public:

    __HbnLookupTable(HbnUnpackedDatabase& db,
        const int kmer_size,
        const int kmer_window,
        const int max_kmer_occ,
        const double rep_frac,
        const int num_threads);  

    __HbnLookupTable(HbnUnpackedDatabase& db,
        ChrIntv* repeat_intva,
        size_t repeat_intvc,
        ChrIntv* non_repeat_intva,
        size_t non_repeat_intvc,
        const int kmer_size,
        const int kmer_window,
        const int max_kmer_occ,
        const double repeat_frac,
        const double non_repeat_frac,
        const int num_threads);

    __HbnLookupTable(): M_hash_table(nullptr), M_offset_list(nullptr), M_offset_list_size(0) {}

    ~__HbnLookupTable() {
        delete M_hash_table;
        delete[] M_offset_list;
    }

    void init_from_repeat_bed_list(HbnUnpackedDatabase& db,
        ChrIntv* intva,
        size_t intvc,
        const int kmer_size,
        const int kmer_window,
        const int max_kmer_occ,
        const double rep_frac,
        const int num_threads);

    void init_from_non_repeat_bed_list(HbnUnpackedDatabase& db,
        ChrIntv* intva,
        size_t intvc,
        const int kmer_size,
        const int kmer_window,
        const int max_kmer_occ,
        const double rep_frac,
        const int num_threads);

    inline void 
    extract_offset_list(const u64 hash, KmerInfo** ol, u64* cnt) const {
        *ol = nullptr;
        *cnt = 0;
        u64 u = M_hash_table->extract_value(hash);
        if (!u) return;
        u64 offset = lktbl_extract_offset(u);
        *cnt = lktbl_extract_cnt(u);
        *ol = M_offset_list + offset;
    }

private:
    HbnHashTable*   M_hash_table;
    KmerInfo*       M_offset_list;
    u64             M_offset_list_size;
};

typedef __HbnLookupTable HbnLookupTable;

#endif // __LOOKUP_TABLE_HPP
