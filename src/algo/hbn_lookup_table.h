#ifndef __HBN_LOOKUP_TABLE_H
#define __HBN_LOOKUP_TABLE_H

#include "seq_loader.h"
#include "../corelib/restrict_enzyme_loci_list.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int m_qidf;
    int m_qidt;
    u64* m_kmer_offset_list;
    u64 m_kmer_offset_list_size;
    u64* m_seq_offset_list;
    u64 m_seq_offset_list_size;
    void* m_hash2offset_map;
} HbnLookupTable;

void HbnLookupTable_KmerOffset2SeqInfo(HbnLookupTable* lktbl, u64 offset, int* id, size_t* start, size_t* end);

u64* 
HbnLookupTable_ExtractKmerOffsetList(HbnLookupTable* lktbl, u64 hash, u64* cnt);

int
HbnLookupTable_qidf(HbnLookupTable* lktbl);

int
HbnLookupTable_qidt(HbnLookupTable* lktbl);

const u64*
HbnLookupTable_SeqOffsetList(HbnLookupTable* lktbl);

HbnLookupTable*
HbnLookupTableFree(HbnLookupTable* lktbl);

HbnLookupTable*
HbnLookupTableNew(SeqReader* db, 
    RestrictEnzymeLociList* relist,
    int kmer_size, 
    int kmer_window,
    double rep_frac,
    int max_kmer_occ, 
    int num_threads);

#ifdef __cplusplus
}
#endif

#endif // __HBN_LOOKUP_TABLE_H