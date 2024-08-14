#ifndef __MAP_ONE_VOLUME_H
#define __MAP_ONE_VOLUME_H

#include "../../corelib/restrict_enzyme_loci_list.hpp"
#include "../../corelib/unpacked_seqdb.hpp"
#include "../../sw/hbn_traceback.hpp"
#include "../../sw/hbn_traceback_aux.h"
#include "chain_align_list.hpp"
#include "hbn_options.hpp"
#include "hbn_outputs.hpp"
#include "hbn_word_finder.hpp"

#include <random>

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
    HbnUnpackedDatabase* subjects;
    HbnUnpackedDatabase* queries;
    HbnLookupTable* lktbl;
    const HbnProgramOptions* opts;
    HbnWordFinder* word_finder;
    std::mt19937* gen;
    std::uniform_int_distribution<>* dist;
    HbnTracebackData* tbck_data;
    RestrictEnzymeLociList* reloci_list;
    QueryVdfEndPointList qvep_list;
    PoreCAlignChainData* pca_chain_data;
    size_t query_id_offset;
    int* qidx;
    pthread_mutex_t* qidx_lock;
    HbnOutputs* out;
} MapThreadData;

MapThreadData*
MapThreadDataNew(int thread_id,
    HbnUnpackedDatabase* subjects,
    HbnUnpackedDatabase* queries,
    HbnWordFinder* word_finder,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    size_t query_id_offset,
    int* qidx,
    pthread_mutex_t* qidx_lock,
    HbnOutputs* out);

MapThreadData*
MapThreadDataFree(MapThreadData* data);

void
map_one_volume(HbnUnpackedDatabase* subjects, 
    HbnUnpackedDatabase* queries,
    HbnWordFinder* word_finder,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    size_t query_id_offset,
    HbnOutputs* out);

#endif // __MAP_ONE_VOLUME_H