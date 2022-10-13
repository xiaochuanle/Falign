#ifndef __HBN_WORD_FINDER_H
#define __HBN_WORD_FINDER_H

#include "hbn_lookup_table.h"
#include "seq_loader.h"
#include "../corelib/gapped_candidate.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
	u64 hash;
    int context;
    int offset;
} DDFKmer;

typedef kvec_t(DDFKmer) vec_ddfk;

void
sort_ddfk_hash_lt(size_t n, DDFKmer* a);

typedef struct {
    int context;
    int qoff;
    idx soff;
} DDFKmerMatch;

typedef kvec_t(DDFKmerMatch) vec_ddfkm;

void
sort_ddfkm_context_lt(size_t n, DDFKmerMatch* a);

void
sort_ddfkm_soff_lt(size_t n, DDFKmerMatch* a);

typedef struct {
    int sid;
    int km_offset;
    int km_cnt;
    int f;
} SubjectKmInfo;

typedef kvec_t(SubjectKmInfo) vec_skmi;

typedef struct {
    SeqReader* db;
    HbnLookupTable* lktbl;
    vec_ddfk ddfk_list;
    vec_ddfkm ddfkm_list;
    vec_skmi skmi_list;
    int kmer_size;
    int kmer_window;
    int min_kmer_cnt;
    vec_u64 hash_list;
} WordFinderThreadData;

WordFinderThreadData*
WordFinderThreadDataNew(SeqReader* db,
    HbnLookupTable* lktbl,
    const int kmer_size,
    const int kmer_window,
    const int min_ddf_score,
    const double max_ddf);

WordFinderThreadData*
WordFinderThreadDataFree(WordFinderThreadData* data);

void 
WordFinderThreadDataClear(WordFinderThreadData* data);

void
hbn_ddfkm(WordFinderThreadData* data, int read_id);

#ifdef __cplusplus
}
#endif

#endif // __HBN_WORD_FINDER_H