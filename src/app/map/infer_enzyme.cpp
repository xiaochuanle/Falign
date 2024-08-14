#include "infer_enzyme.hpp"

#include "../../corelib/pdqsort.h"
#include "../../corelib/restrict_enzyme_loci_list.hpp"
#include "../../sw/hbn_traceback.hpp"
#include "extend_hit_list.hpp"
#include "hbn_word_finder.hpp"
#include "make_candidate_kmer_chain.hpp"

#include <limits>
#include <map>
#include <mutex>
#include <random>
#include <vector>

using namespace std;

#ifndef U32_MAX
constexpr const u32 U32_MAX = numeric_limits<u32>::max();
#endif

static const char* kEnzymeNameList[] = {
    "DpnII",
    "HindIII",
    "NcoI",
    "NlaIII"
};

static const char* kEnzymeList[] = {
    "^GATC",
    "A^AGCTT",
    "C^CATGG",
    "CATG^"
};
static constexpr const int kEnzymeListSize = 4;

class EnzymeList
{
public:
    EnzymeList(const HbnProgramOptions* options);
    ~EnzymeList();

    RestrictEnzyme** get_enzyme_list() {
        return M_enzyme_list.data();
    }

    void add_one_enzyme_count(u32 hash) {
        std::lock_guard<std::mutex> lg(M_mutex);
        auto pos = M_enzyme_stats.find(hash);
        hbn_assert(pos != M_enzyme_stats.end());
        ++pos->second.second;
    }

    void dump_stats();

    void find_enzyme(const u8* query,
        const int qoff, 
        const int qend,
        const int qsize,
        const u8* subject,
        const int soff,
        const int send,
        const int ssize);

    std::string infer_enzyme() {
        dump_stats();
        
        pair<int, size_t> cnts[kEnzymeListSize];
        for (int i = 0; i < kEnzymeListSize; ++i) {
            u32 hash = M_enzyme_list[i]->enzyme_hash;
            auto pos = M_enzyme_stats.find(hash);
            hbn_assert(pos != M_enzyme_stats.end());
            cnts[i] = pair<int, size_t>(i, pos->second.second);
        }
        pdqsort(cnts, cnts + kEnzymeListSize, [](const pair<int, size_t>& x, const pair<int, size_t>& y) { return x.second > y.second; });
        int c0 = cnts[0].second;
        int c1 = cnts[1].second;
        hbn_assert(c0 >= c1);
        string enzyme;
        if (c0 > c1 * 2) enzyme = kEnzymeList[cnts[0].first];
        return enzyme; 
    }

private:
    const HbnProgramOptions*  M_options;
    std::mutex  M_mutex;
    std::vector<RestrictEnzyme*> M_enzyme_list;
    std::map<u32, std::pair<RestrictEnzyme*, size_t>>   M_enzyme_stats;
};

EnzymeList::EnzymeList(const HbnProgramOptions* options)
{
    M_options = options;
    for (int i = 0; i < kEnzymeListSize; ++i) {
        RestrictEnzyme* enzyme = new RestrictEnzyme();
        RestrictEnzyme_Init(kEnzymeList[i], enzyme);
        M_enzyme_list.push_back(enzyme);
        M_enzyme_stats[enzyme->enzyme_hash] = pair<RestrictEnzyme*, size_t>(enzyme, 0);
    }
}

EnzymeList::~EnzymeList()
{
    for (auto enzyme : M_enzyme_list) delete enzyme;
}

void EnzymeList::dump_stats()
{
    char buf[256];
    for (int i = 0; i < kEnzymeListSize; ++i) {
        u32 hash = M_enzyme_list[i]->enzyme_hash;
        auto pos = M_enzyme_stats.find(hash);
        hbn_assert(pos != M_enzyme_stats.end());
        int cnt = pos->second.second;
        snprintf(buf, 256, "%s(%s)", kEnzymeNameList[i], kEnzymeList[i]);
        print_fixed_width_string(stderr, buf, 20);
        fprintf(stderr, "%d\n", cnt);
    }
}

static bool 
s_enzyme_exists(RestrictEnzyme* enzyme,
    const int kEndDist,
    const int kEnzymeFlankingBases,
    const u8* query,
    int qoff,
    int qsize,
    const u8* subject,
    int soff,
    int ssize)
{
    bool r = (qoff < kEndDist) || (qsize - qoff < kEndDist)
        || (soff < kEndDist) || (ssize - soff < kEndDist);
    if (r) return false;
    
    int qfrom = qoff - kEnzymeFlankingBases;
    int qto = qoff + kEnzymeFlankingBases;
    int p = qfrom;
    bool q_exists = false;
    u32 hash = 0;
    for (; p < qfrom + enzyme->enzyme_size - 1; ++p) hash = (hash << 2) | query[p];
    for (; p < qto; ++p) {
        hash = (hash << 2) | query[p];
        hash = hash & enzyme->enzyme_mask;
        if (hash == enzyme->enzyme_hash) {
            q_exists = true;
            break;
        }
    }
    if (!q_exists) return false;

    int sfrom = soff - kEnzymeFlankingBases;
    int sto = soff + kEnzymeFlankingBases;
    p = sfrom;
    bool s_exists = false;
    hash = 0;
    for (; p < sfrom + enzyme->enzyme_size - 1; ++p) hash = (hash << 2) | subject[p];
    for (; p < sto; ++p) {
        hash = (hash << 2) | subject[p];
        hash = hash & enzyme->enzyme_mask;
        if (hash == enzyme->enzyme_hash) {
            s_exists = true;
            break;
        }        
    }
    if (!s_exists) return false;

    return true;
}

void EnzymeList::find_enzyme(const u8* query,
    const int qoff, 
    const int qend,
    const int qsize,
    const u8* subject,
    const int soff,
    const int send,
    const int ssize)
{
    RestrictEnzyme** enzymes = get_enzyme_list();
    u32 L_exist[kEnzymeListSize]; fill(L_exist, L_exist + kEnzymeListSize, U32_MAX);
    u32 R_exist[kEnzymeListSize]; fill(R_exist, R_exist + kEnzymeListSize, U32_MAX);
    for (int i = 0; i < kEnzymeListSize; ++i) {
	    if (s_enzyme_exists(enzymes[i], M_options->ei_end_dist, M_options->ei_flanking_bases, query, qoff, qsize, subject, soff, ssize)) L_exist[i] = enzymes[i]->enzyme_hash;
    }
    for (int i = 0; i < kEnzymeListSize; ++i) {
	    if (s_enzyme_exists(enzymes[i], M_options->ei_end_dist, M_options->ei_flanking_bases, query, qend, qsize, subject, send, ssize)) R_exist[i] = enzymes[i]->enzyme_hash;
    }

    for (int i = 0; i < kEnzymeListSize; ++i) {
	    if (L_exist[i] == U32_MAX || R_exist[i] == U32_MAX) continue;
	    add_one_enzyme_count(L_exist[i]);
    }
}

////////////////////////////

struct EnzymeInferenceThreadWorkData
{
    int thread_id;
    HbnUnpackedDatabase* subjects;
    HbnUnpackedDatabase* queries;
    const HbnProgramOptions* opts;
    HbnWordFinder* word_finder;
    HbnTracebackData* tbck_data;
    int* qidx;
    pthread_mutex_t* qidx_lock;
    EnzymeList* enzyme_list;
};

EnzymeInferenceThreadWorkData* 
EnzymeInferenceThreadWorkDataNew(int thread_id,
    HbnUnpackedDatabase* subjects,
    HbnUnpackedDatabase* queries,
    HbnWordFinder* word_finder,
    const HbnProgramOptions* opts,
    int* qidx,
    pthread_mutex_t* qidx_lock,
    EnzymeList* enzyme_list)
{
    EnzymeInferenceThreadWorkData* data = new EnzymeInferenceThreadWorkData();
    data->thread_id = thread_id;
    data->subjects = subjects;
    data->queries = queries;
    data->opts = opts;
    data->word_finder = word_finder;
    data->tbck_data = new HbnTracebackData();
    data->qidx = qidx;
    data->qidx_lock = qidx_lock;
    data->enzyme_list = enzyme_list;

    return data;
}

EnzymeInferenceThreadWorkData*
EnzymeInferenceThreadWorkDataFree(EnzymeInferenceThreadWorkData* data)
{
    delete data->tbck_data;
    delete data;
    return NULL;
}

static inline int
EnzymeInferenceThreadWorkData_GetNextQuery(EnzymeInferenceThreadWorkData* data)
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

static void
s_enzyme_inference_detect_candidates(EnzymeInferenceThreadWorkData* data, 
    const int qidx,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    mt19937* gen,
    uniform_int_distribution<>* dist,
    vector<HbnInitHit>& hit_list)
{
    hit_list.clear();

    HbnWordFinder* word_finder = data->word_finder;
    vector<KmerMatch> km_list;
    word_finder->extract_kmer_matches(qidx, fwd_query, rev_query, query_size, km_list, gen, dist);
    vector<HbnInitHit> l_hit_list;
    vector<HbnKmerMatch> fwd_seeds, rev_seeds;

    size_t all_nkm = km_list.size();
    KmerMatch* all_kma = km_list.data();
    size_t i = 0;
    while (i < all_nkm) {
        size_t j = i + 1;
        while (j < all_nkm && all_kma[i].sid == all_kma[j].sid) ++j;
        const int sid = all_kma[i].sid;
        KmerMatch* kma = all_kma + i;
        int kmc = j - i;
        i = j;

        const int subject_size = data->subjects->SeqSize(sid);
        fwd_seeds.clear();
        rev_seeds.clear();
        for (int k = 0; k < kmc; ++k) {
            HbnKmerMatch hkm;
            hkm.qoff = kma[k].qoff;
            hkm.soff = kma[k].soff;
            hkm.length = data->opts->kmer_size;
            if (kma[k].qdir == REV) {
                rev_seeds.push_back(hkm);
            } else {
                fwd_seeds.push_back(hkm);
            }
        }

        HbnKmerMatch* hkma = fwd_seeds.data();
        int hkmc = fwd_seeds.size();
        l_hit_list.clear();
        make_candidate_kmer_chain(hkma,
            hkmc,
            qidx,
            FWD,
            query_size,
            sid,
            subject_size,
            data->opts->kmer_dist,
            data->opts->ddf,
            data->opts->chain_score,
            l_hit_list);
        HbnInitHit* hit_array = l_hit_list.data();
        int hit_count = l_hit_list.size();
        for (int x = 0; x < hit_count; ++x) {
            HbnInitHit hit = l_hit_list[x];
            hit_list.push_back(hit);
        }

        hkma = rev_seeds.data();
        hkmc = rev_seeds.size();
        l_hit_list.clear();
        make_candidate_kmer_chain(hkma,
            hkmc,
            qidx,
            REV,
            query_size,
            sid,
            subject_size,
            data->opts->kmer_dist,
            data->opts->ddf,
            data->opts->chain_score,
            l_hit_list);
        hit_array = l_hit_list.data();
        hit_count = l_hit_list.size();
        for (int x = 0; x < hit_count; ++x) {
            HbnInitHit hit = l_hit_list[x];
            hit_list.push_back(hit);
        }
    }
}

static void
enzyme_inference_for_one_query(EnzymeInferenceThreadWorkData* data, 
    const char* query_name,
    const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    mt19937* gen,
    uniform_int_distribution<>* dist)
{
    //HBN_LOG("mappig query %d:%s:%d", query_id, query_name, query_size);

    vector<HbnInitHit> hit_list;
    s_enzyme_inference_detect_candidates(data, query_id, fwd_query, rev_query, query_size, gen, dist, hit_list);
    if (hit_list.empty()) return;

    HbnInitHit* hita = hit_list.data();
    int hitc = hit_list.size();
    pdqsort(hita, hita + hitc, [](const HbnInitHit& x, const HbnInitHit& y) { return x.score > y.score; });

    for (int i = 0; i < data->opts->ei_num_hits && i < hitc; ++i) {
        HbnInitHit* hit = hita + i;
        const u8* query = (hit->qdir == FWD) ? fwd_query : rev_query;
        const u8* subject = data->subjects->GetSequence(hit->sid);
        const int subject_size = data->subjects->SeqSize(hit->sid);
        if (!align_subseq_enzyme_inference(query_id, hit->qdir, fwd_query, rev_query, query_size, 
            hit->sid, subject, subject_size,
            hit->qbeg, hit->qend, hit->sbeg, hit->send, data->opts->perc_identity, data->tbck_data)) continue;

        int qb = data->tbck_data->qoff;
        int qe = data->tbck_data->qend;
        int sb = data->tbck_data->soff;
        int se = data->tbck_data->send;
        //double pi = data->tbck_data->ident_perc;
        if (qe - qb < data->opts->ei_frag_size || se - sb < data->opts->ei_frag_size) return;
        //fprintf(stderr, "[%d, %d, %d] x [%d, %d, %d], %g\n", qb, qe, query_size,
        //    sb, se, subject_size, pi);
    
        data->enzyme_list->find_enzyme(query, qb, qe, query_size, subject, sb, se, subject_size);
    }
}

static void*
s_enzyme_inference_thread(void* params)
{
    EnzymeInferenceThreadWorkData* data = (EnzymeInferenceThreadWorkData*)(params);

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

    random_device rd;
    mt19937 gen(rd());
    uniform_int_distribution<> dist(0, 1);

    while ((query_id = EnzymeInferenceThreadWorkData_GetNextQuery(data)) != -1) {
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
        enzyme_inference_for_one_query(data, query_name, query_id, fwd_query, rev_query, query_size, &gen, &dist);
    }
    return NULL;
}

///////////////////////////

std::string
infer_enzyme_mt(const HbnProgramOptions* options,
    HbnUnpackedDatabase* subjects,
    HbnUnpackedDatabase* queries,
    HbnWordFinder* word_finder)
{
    const int num_threads = options->num_threads;
    EnzymeList enzyme_list(options);
    EnzymeInferenceThreadWorkData* data_array[num_threads];
    pthread_t jobids[num_threads];
    int qidx = 0;
    pthread_mutex_t qidx_lock = PTHREAD_MUTEX_INITIALIZER;
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = EnzymeInferenceThreadWorkDataNew(i, subjects, queries, word_finder, options, &qidx, &qidx_lock, &enzyme_list);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobids + i, NULL, s_enzyme_inference_thread, data_array[i]);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobids[i], NULL);
    }
    for (int i = 0; i < num_threads; ++i) {
        data_array[i] = EnzymeInferenceThreadWorkDataFree(data_array[i]);
    }

    return enzyme_list.infer_enzyme();
}