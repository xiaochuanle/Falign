#ifndef __MAP_ONE_VOLUME_H
#define __MAP_ONE_VOLUME_H

#include "../../algo/seq_loader.h"
#include "../../algo/hbn_lookup_table.h"
#include "../../algo/hbn_word_finder.h"
#include "../../corelib/restrict_enzyme_loci_list.h"
#include "chain_align_list.hpp"
#include "hbn_options.h"

#include <pthread.h>

typedef struct {
    int skmi_idx;
    int sid;
    int score;
    int hit_offset;
    int hit_cnt;
    int qbeg, qend;
    int sbeg, send;
} CandidateSubject;

typedef struct {
    int thread_id;
    SeqReader* subjects;
    SeqReader* queries;
    HbnLookupTable* lktbl;
    const HbnProgramOptions* opts;
    WordFinderThreadData* word_data;
    HbnTracebackData* tbck_data;
    RestrictEnzymeLociList* reloci_list;
    QueryVdfEndPointList qvep_list;
    PoreCAlignChainData* pca_chain_data;
    int* qidx;
    pthread_mutex_t* qidx_lock;
    FILE* out;
    pthread_mutex_t* out_lock;
} MapThreadData;

MapThreadData*
MapThreadDataNew(int thread_id,
    SeqReader* subjects,
    SeqReader* queries,
    HbnLookupTable* lktbl,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    int* qidx,
    pthread_mutex_t* qidx_lock,
    FILE* out,
    pthread_mutex_t* out_lock);

MapThreadData*
MapThreadDataFree(MapThreadData* data);

void
map_one_volume(SeqReader* subjects, 
    SeqReader* queries,
    HbnLookupTable* lktbl,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    FILE* out);

#endif // __MAP_ONE_VOLUME_H