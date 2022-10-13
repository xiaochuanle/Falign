#include "map_one_volume.hpp"

#include "align_one_read.hpp"
#include "trim_overlap_subseq.hpp"
#include "../../algo/hbn_traceback_aux.h"

#include <algorithm>
#include <vector>

using namespace std;

static const int kQueryBatchSize = 10;

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
    pthread_mutex_t* out_lock)
{
    MapThreadData* data = (MapThreadData*)calloc(1, sizeof(MapThreadData));
    data->thread_id = thread_id;
    data->subjects = subjects;
    data->queries = queries;
    data->lktbl = lktbl;
    data->opts = opts;
    data->word_data = WordFinderThreadDataNew(queries,
            lktbl,
            opts->kmer_size,
            1,
            opts->chain_score,
            opts->ddf);
    data->tbck_data = HbnTracebackDataNew();
    data->reloci_list = reloci_list;
    QueryVdfEndPointList_Init(&data->qvep_list);
    data->pca_chain_data = new PoreCAlignChainData();
    data->qidx = qidx;
    data->qidx_lock = qidx_lock;
    data->out = out;
    data->out_lock = out_lock;

    return data;
}

MapThreadData*
MapThreadDataFree(MapThreadData* data)
{
    data->word_data = WordFinderThreadDataFree(data->word_data);
    data->tbck_data = HbnTracebackDataFree(data->tbck_data);
    QueryVdfEndPointList_Destroy(&data->qvep_list);
    delete data->pca_chain_data;
    free(data);
    return NULL;
}

static BOOL
MapThreadData_GetNextQueryBatch(MapThreadData* data, int* _qidx_from, int* _qidx_to)
{
    int from = 0, to = 0;
    const int num_queries = SeqReader_NumSeqs(data->queries);
    pthread_mutex_lock(data->qidx_lock);
    from = *data->qidx;
    *data->qidx += kQueryBatchSize;
    pthread_mutex_unlock(data->qidx_lock);
    if (from >= num_queries) return FALSE;
    to = hbn_min(from + kQueryBatchSize, num_queries);
    *_qidx_from = from;
    *_qidx_to = to;
    return TRUE;
}

static void*
s_map_thread(void* params)
{
    MapThreadData* data = (MapThreadData*)(params);
    int qidx_from, qidx_to;
    ks_dinit(out_buf);
    while (MapThreadData_GetNextQueryBatch(data, &qidx_from, &qidx_to)) {
        for (int i = qidx_from; i < qidx_to; ++i) {
            align_one_read(data, i, &out_buf);
        }
        pthread_mutex_lock(data->out_lock);
        hbn_fwrite(ks_s(out_buf), 1, ks_size(out_buf), data->out);
        pthread_mutex_unlock(data->out_lock);
        ks_clear(out_buf);
        data->pca_chain_data->release_soa();
        //if ((qidx_to % 1000) == 0) HBN_LOG("%d queries mapped", qidx_to);
    }
    ks_destroy(out_buf);
    return NULL;
}

void
map_one_volume(SeqReader* subjects, 
    SeqReader* queries,
    HbnLookupTable* lktbl,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    FILE* out)
{
    const int num_threads = opts->num_threads;
    MapThreadData* data_array[num_threads];
    pthread_t jobids[num_threads];
    int qidx = 0;
    pthread_mutex_t qidx_lock = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t out_lock = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = MapThreadDataNew(i, subjects, queries, lktbl, opts, reloci_list, &qidx, &qidx_lock, out, &out_lock);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobids + i, NULL, s_map_thread, data_array[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobids[i], NULL);
    }
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = MapThreadDataFree(data_array[i]);
    }
}
