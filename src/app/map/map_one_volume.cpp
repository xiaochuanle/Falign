#include "map_one_volume.hpp"

#include "align_one_read.hpp"
#include "trim_overlap_subseq.hpp"

#include <algorithm>
#include <vector>

using namespace std;

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
    HbnOutputs* out)
{
    MapThreadData* data = new MapThreadData();
    data->thread_id = thread_id;
    data->subjects = subjects;
    data->queries = queries;
    data->word_finder = word_finder;
    data->opts = opts;
    data->tbck_data = new HbnTracebackData();
    data->reloci_list = reloci_list;
    data->pca_chain_data = new PoreCAlignChainData();
    data->query_id_offset = query_id_offset;
    data->qidx = qidx;
    data->qidx_lock = qidx_lock;
    data->out = out;

    random_device rd;
    data->gen = new mt19937(rd());
    data->dist = new uniform_int_distribution<>(0, 1);

    return data;
}

MapThreadData*
MapThreadDataFree(MapThreadData* data)
{
    delete data->gen;
    delete data->dist;
    delete data->tbck_data;
    delete data->pca_chain_data;
    delete data;
    return NULL;
}

static inline int
MapThreadData_GetNextQuery(MapThreadData* data)
{
    int qidx = -1;
    const int num_queries = data->queries->NumSeqs();
    pthread_mutex_lock(data->qidx_lock);
    qidx = *data->qidx;
    *data->qidx += 1;
    pthread_mutex_unlock(data->qidx_lock);
    if (qidx >= num_queries) qidx = -1;
    return qidx;
}

static void
s_setup_query_sequences(const u8* fwd_query, 
    const char* raw_qv,
    const int query_size,
    vector<u8>& _rev_query,
    const u8*& rev_query,
    vector<char>& _fwd_qv,
    const char*& fwd_qv,
    vector<char>& _rev_qv,
    const char*& rev_qv,
    vector<char>& _fwd_raw_query,
    const char*& fwd_raw_query,
    vector<char>& _rev_raw_query,
    const char*& rev_raw_query)
{
    _rev_query.assign(fwd_query, fwd_query + query_size);
    reverse(_rev_query.begin(), _rev_query.end());
    for (auto& c : _rev_query) c = 3 - c;
    rev_query = _rev_query.data();

    fwd_qv = nullptr;
    if (raw_qv) {
        _fwd_qv.assign(raw_qv, raw_qv + query_size);
        for (auto& c : _fwd_qv) c -= 33;
        fwd_qv = _fwd_qv.data();
    }

    rev_qv = nullptr;
    if (raw_qv) {
        _rev_qv.assign(_fwd_qv.rbegin(), _fwd_qv.rend());
        rev_qv = _rev_qv.data();
    }

    _fwd_raw_query.clear();
    _rev_raw_query.clear();
    for (int i = 0; i < query_size; ++i) {
        int fc = fwd_query[i];
        _fwd_raw_query.push_back(DECODE_RESIDUE(fc));
        int rc = rev_query[i];
        _rev_raw_query.push_back(DECODE_RESIDUE(rc));
    }
    fwd_raw_query = _fwd_raw_query.data();
    rev_raw_query = _rev_raw_query.data();
}

static void*
s_map_thread(void* params)
{
    MapThreadData* data = (MapThreadData*)(params);
    vector<PoreCAlign> all_pca_list;
    TrimPcaList trim_pca_list;

    const char* query_name;
    int query_id;
    const u8* fwd_query;
    const char* raw_qv;
    vector<u8> _rev_query;
    vector<char> _fwd_qv;
    vector<char> _rev_qv;
    vector<char> _fwd_raw_query;
    vector<char> _rev_raw_query;
    const u8* rev_query;
    const char* fwd_qv;
    const char* rev_qv;
    const char* fwd_raw_query;
    const char* rev_raw_query;

    int soa_cnt = 0;
    while ((query_id = MapThreadData_GetNextQuery(data)) != -1) {
        query_name = data->queries->SeqName(query_id);
        fwd_query = data->queries->GetSequence(query_id);
        raw_qv = data->queries->GetBaseQualities(query_id);
        const int query_size = data->queries->SeqSize(query_id);
        s_setup_query_sequences(fwd_query, raw_qv, query_size,
            _rev_query, rev_query,
            _fwd_qv, fwd_qv,
            _rev_qv, rev_qv,
            _fwd_raw_query, fwd_raw_query,
            _rev_raw_query, rev_raw_query);
        align_one_read(data, query_name, query_id, fwd_query, rev_query, query_size, all_pca_list, trim_pca_list);
        query_id += data->query_id_offset;
        data->out->dump(data->reloci_list,
            data->subjects,
            query_name,
            query_id,
            fwd_raw_query,
            rev_raw_query,
            fwd_qv,
            rev_qv,
            query_size,
            all_pca_list.data(),
            all_pca_list.size(),
            trim_pca_list);
        if (++soa_cnt == 10)  {
            data->pca_chain_data->release_soa();
            soa_cnt = 0;
        }
    }
    return NULL;
}

void
map_one_volume(HbnUnpackedDatabase* subjects, 
    HbnUnpackedDatabase* queries,
    HbnWordFinder* word_finder,
    const HbnProgramOptions* opts,
    RestrictEnzymeLociList* reloci_list,
    const size_t query_id_offset,
    HbnOutputs* out)
{
    const int num_threads = opts->num_threads;
    MapThreadData* data_array[num_threads];
    pthread_t jobids[num_threads];
    int qidx = 0;
    pthread_mutex_t qidx_lock = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = MapThreadDataNew(i, subjects, queries, word_finder, opts, reloci_list, query_id_offset, &qidx, &qidx_lock, out);
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
