#include "../../corelib/lookup_table.hpp"
#include "../../corelib/pdqsort.h"
#include "../../corelib/parasort.h"
#include "window_masker.hpp"

#include <mutex>

using namespace std;

class KmerExtractionWorkThreadData
{
public:
    KmerExtractionWorkThreadData(HbnUnpackedDatabase* reference,
        ChrIntv* intva,
        size_t intvc,
        int kmer_size,
        int kmer_window,
        vector<KmerHashAndOffset>* khao_list)
    {
        M_reference = reference;
        M_kmer_size = kmer_size;
        M_kmer_window = kmer_window;

        M_intva = intva;
        M_intvc = intvc;
        M_intvi = 0;

        M_khao_list = khao_list;
    }

    void add_kmers(KmerHashAndOffset* a, size_t c)
    {
        lock_guard<mutex> _(M_kmer_mutex);
        M_khao_list->insert(M_khao_list->end(), a, a + c);
    }

    HbnUnpackedDatabase* reference()
    {
        return M_reference;
    }

    int kmer_size()
    {
        return M_kmer_size;
    }

    int kmer_window()
    {
        return M_kmer_window;
    }

    bool get_next_intv(ChrIntv& intv, size_t& bed_id) {
        lock_guard<mutex> _(M_intvi_mutex);
        bool r = false;
        if (M_intvi < M_intvc) {
            bed_id = M_intvi;
            intv = M_intva[M_intvi];
            ++M_intvi;
            r = true;
        }
        return r;
    }

private:
    HbnUnpackedDatabase*    M_reference;
    int                     M_kmer_size;
    int                     M_kmer_window;

    ChrIntv*                M_intva;
    size_t                  M_intvc;
    size_t                  M_intvi;
    mutex                   M_intvi_mutex;

    mutex                   M_kmer_mutex;
    vector<KmerHashAndOffset>*  M_khao_list;
};

static void*
s_kmer_extract_thread(void* params)
{
    KmerExtractionWorkThreadData* data = (KmerExtractionWorkThreadData*)(params);
    HbnUnpackedDatabase* reference = data->reference();
    int kmer_size = data->kmer_size();
    int kmer_window = data->kmer_window();
    vector<KmerHashAndOffset> khao_list;
    ChrIntv intv;
    size_t bed_id;
    while (data->get_next_intv(intv, bed_id)) {
        khao_list.clear();
        const u8* chr = reference->GetSequence(intv.chr_id);
        extract_hash_values_for_one_read(intv.chr_id, chr + intv.from, intv.to - intv.from, intv.from, kmer_size, kmer_window, 0, U64_MAX, khao_list);
#if DUMP_REPEAT_BED_KMER_STATS
        for (auto& khao : khao_list) khao.bed_id = bed_id;
#endif
        data->add_kmers(khao_list.data(), khao_list.size());
    }

    return nullptr;
}

int
resolve_kmer_occ_cutoff(KmerHashAndOffset* khao_array, 
    u64 khao_array_size, 
    int num_threads,
    int kmer_size,
    double rep_frac, 
    int _max_kmer_occ);

struct XXX_KmerStats
{
u64 total_kmers;
u64 distinct_kmers;
u64 removed_kmers;
u64 removed_distinct_kmers;

    XXX_KmerStats() {
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

void extract_kmers_from_bed_list_mt(HbnUnpackedDatabase* reference,
    ChrIntv* intva,
    size_t intvc,
    int num_threads,
    int kmer_size,
    int kmer_window,
    std::vector<KmerHashAndOffset>& khao_list)
{
    pthread_t jobs[num_threads];
    KmerExtractionWorkThreadData data(reference, intva, intvc, kmer_size, kmer_window, &khao_list);
    for (int i = 0; i < num_threads; ++i) {
        pthread_create(jobs + i, nullptr, s_kmer_extract_thread, &data);
    }
    for (int i = 0; i < num_threads; ++i) {
        pthread_join(jobs[i], nullptr);
    }
}

/////////////////

static void
s_extract_repeat_region_kmers(HbnUnpackedDatabase& db,
    ChrIntv* intva,
    size_t intvc,
    const int kmer_size,
    const int kmer_window,
    const int max_kmer_occ,
    const double rep_frac,
    const int num_threads,
    vector<KmerHashAndOffset>& khao_list)
{
    extract_kmers_from_bed_list_mt(&db, intva, intvc, num_threads, kmer_size, kmer_window, khao_list);

    {
        KmerHashAndOffset* a = khao_list.data();
        size_t ac = khao_list.size();
        parasort(ac, a, num_threads);
        size_t i = 0;
        XXX_KmerStats stats;
        int kmer_occ_cutoff = resolve_kmer_occ_cutoff(a, ac, num_threads, kmer_size, rep_frac, 1);
        while (i < ac) {
            size_t j = i + 1;
            while (j < ac && a[i].hash == a[j].hash) ++j;
            if (j - i > kmer_occ_cutoff) { 
                for (size_t k = i; k < j; ++k) a[k].hash = U64_MAX;
                stats.add(j - i, 1, j - i, 1); 
                i = j; 
                continue; 
            }
            stats.add(j - i, 1, 0, 0);
            i = j;
        }
        stats.stats();

        size_t n = 0;
        for (i = 0; i < ac; ++i) if (a[i].hash != U64_MAX) a[n++] = a[i];
        khao_list.resize(n);
    }

#if DUMP_REPEAT_BED_KMER_STATS
    {
        size_t N = intvc;
        u8* bed_ids = (u8*)calloc(N, sizeof(u8));
        for (auto& khao : khao_list) bed_ids[khao.bed_id] = 1;
        hbn_dfopen(out, "repeat-beds.txt", "w");
        size_t intvs = 0, bases = 0;
        for (size_t i = 0; i < N; ++i) {
            if (bed_ids[i]) continue;
            int chr_id = intva[i].chr_id;
            const char* chr_name = db.SeqName(chr_id);
            int from = intva[i].from;
            int to = intva[i].to;
            int L = to - from;
            const u8* chr = db.GetSequence(chr_id);
            fprintf(out, "%s:%d-%d-%d\n", chr_name, from, to, L);
            for (int k = from; k < to; ++k) {
                int c = chr[k];
                c = DECODE_RESIDUE(c);
                fprintf(out, "%c", c);
            }
            fprintf(out, "\n");
            ++intvs;
            bases += L;
        }
        hbn_fclose(out);
        free(bed_ids);
        fprintf(stderr, "Repeat beds: %zu, bases: %zu\n", intvs,bases);
    }
#endif
}

static void
s_extract_non_repeat_region_kmers(HbnUnpackedDatabase& db,
    ChrIntv* intva,
    size_t intvc,
    const int kmer_size,
    const int kmer_window,
    const int max_kmer_occ,
    const double rep_frac,
    const int num_threads,
    vector<KmerHashAndOffset>& khao_list)
{
    extract_kmers_from_bed_list_mt(&db, intva, intvc, num_threads, kmer_size, kmer_window, khao_list);

    {
        KmerHashAndOffset* a = khao_list.data();
        size_t ac = khao_list.size();
        parasort(ac, a, num_threads);
        size_t i = 0;
        XXX_KmerStats stats;
        int kmer_occ_cutoff = resolve_kmer_occ_cutoff(a, ac, num_threads, kmer_size, rep_frac, max_kmer_occ);
        while (i < ac) {
            size_t j = i + 1;
            while (j < ac && a[i].hash == a[j].hash) ++j;
            if (j - i > kmer_occ_cutoff) { 
                for (size_t k = i; k < j; ++k) a[k].hash = U64_MAX;
                stats.add(j - i, 1, j - i, 1); 
                i = j; 
                continue; 
            }
            stats.add(j - i, 1, 0, 0);
            i = j;
        }
        stats.stats();

        size_t n = 0;
        for (i = 0; i < ac; ++i) if (a[i].hash != U64_MAX) a[n++] = a[i];
        khao_list.resize(n);
    }
}

__HbnLookupTable::__HbnLookupTable(HbnUnpackedDatabase& db,
    ChrIntv* repeat_intva,
    size_t repeat_intvc,
    ChrIntv* non_repeat_intva,
    size_t non_repeat_intvc,
    const int kmer_size,
    const int kmer_window,
    const int max_kmer_occ,
    const double repeat_frac,
    const double non_repeat_frac,
    const int num_threads)
{
    vector<KmerHashAndOffset> repeat_khao_list;
    if (repeat_intvc > 0)
    s_extract_repeat_region_kmers(db,
        repeat_intva,
        repeat_intvc,
        kmer_size,
        kmer_window,
        max_kmer_occ,
        repeat_frac,
        num_threads,
        repeat_khao_list);
    
    vector<KmerHashAndOffset> non_repeat_khao_list;
    if (non_repeat_intvc > 0)
    s_extract_non_repeat_region_kmers(db,
        non_repeat_intva,
        non_repeat_intvc,
        kmer_size,
        kmer_window,
        max_kmer_occ,
        non_repeat_frac,
        num_threads,
        non_repeat_khao_list);    

    non_repeat_khao_list.insert(non_repeat_khao_list.end(), repeat_khao_list.begin(), repeat_khao_list.end());
    KmerHashAndOffset* a = non_repeat_khao_list.data();
    size_t ac = non_repeat_khao_list.size();
    parasort(ac, a, num_threads);

    M_hash_table = new HbnHashTable(kmer_size);
    M_offset_list = (KmerInfo*)calloc(ac, sizeof(KmerInfo));
    M_offset_list_size = ac;
    u64 kmer_idx = 0;

    size_t distinct_kmers = 0;
    size_t ii = 0;
    while (ii < ac) {
        size_t jj = ii + 1;
        while (jj < ac && a[ii].hash == a[jj].hash) ++jj;
        ++distinct_kmers;
        ii = jj;
    } 
    HBN_LOG("Distinct kmers: %zu", distinct_kmers);
    HBN_LOG("Total kmers: %zu", ac);
    M_hash_table->reserve(distinct_kmers);

    ii = 0;
    HBN_LOG("Add kmers into hash table");
    while (ii < ac) {
        size_t jj = ii + 1;
        while (jj < ac && a[ii].hash == a[jj].hash) ++jj;
        u64 cnt = jj - ii;
        u64 offset = kmer_idx;
        u64 u = lktbl_pack_hash_value(offset, cnt);
        M_hash_table->add_one_hash_and_value(a[ii].hash, u);
        for (size_t kk = ii; kk < jj; ++kk, ++kmer_idx) M_offset_list[kmer_idx] = a[kk].ki;
        ii = jj;
    }
    HBN_LOG("Done.");
    fprintf(stderr, "\n");   
}

void 
kmer_cov_stats(const char* referenece_path, 
    const char* repeat_bed_path,
    int repeat_kmer_size,
    int repeat_kmer_window,
    int repeat_max_kmer_occ,
    double repeat_rep_frac,
    int non_repeat_kmer_size,
    int non_repeat_kmer_window,
    int non_repeat_max_kmer_occ,
    double non_repeat_rep_frac,
    int num_threads)
{
    HbnUnpackedDatabase db(referenece_path);
    db.load_next_batch();

    SeqName2IdMap refname2id;
    const int num_chr = db.NumSeqs();
    for (int i = 0; i < num_chr; ++i) refname2id.add_one_name(db.SeqName(i));
    refname2id.build_name2id_map();

    vector<ChrIntv> repeat_intv_list;
    load_bed_list(repeat_bed_path, refname2id, repeat_intv_list);

    vector<ChrIntv> non_repeat_intv_list;
    compute_complement_intv_list(repeat_intv_list.data(), repeat_intv_list.size(), db, non_repeat_intv_list);

    vector<KmerHashAndOffset> repeat_khao_list;
    s_extract_repeat_region_kmers(db,
        repeat_intv_list.data(),
        repeat_intv_list.size(),
        repeat_kmer_size,
        repeat_kmer_window,
        repeat_max_kmer_occ,
        repeat_rep_frac,
        num_threads,
        repeat_khao_list);
    
    vector<KmerHashAndOffset> non_repeat_khao_list;
    s_extract_non_repeat_region_kmers(db,
        non_repeat_intv_list.data(),
        non_repeat_intv_list.size(),
        non_repeat_kmer_size,
        non_repeat_kmer_window,
        non_repeat_max_kmer_occ,
        non_repeat_rep_frac,
        num_threads,
        non_repeat_khao_list);   

    int chr_id = 22;
    int from = 52820000;
    int to =   65930000;
    int chr_size = db.SeqSize(chr_id);
    const char* chr_name = db.SeqName(chr_id);
    u8* cov_list = (u8*)calloc(chr_size, 1);
    size_t xxx = 0;
    for (auto& khao : repeat_khao_list) {
        if (khao.ki.seq_id != chr_id) continue;
        if (khao.ki.seq_offset < from || khao.ki.seq_offset >= to) continue;
        for (int i = 0; i < repeat_kmer_size; ++i) cov_list[khao.ki.seq_offset+i] = 1;
        ++xxx;
    }
    HBN_LOG("repeat khao: %zu", xxx);
    xxx = 0;
    for (auto& khao : non_repeat_khao_list) {
        if (khao.ki.seq_id != chr_id) continue;
        if (khao.ki.seq_offset < from || khao.ki.seq_offset >= to) continue;
        for (int i = 0; i < repeat_kmer_size; ++i) cov_list[khao.ki.seq_offset+i] = 1;
        ++xxx;
    }
    HBN_LOG("non-repeat khao: %zu", xxx);
    int cov = 0;
    for (int i = from; i < to; ++i) cov += cov_list[i];
    double frac = 1.0 * cov / (to - from);
    HBN_LOG("%s, cov-frac = %g", chr_name, frac);

    int max_cov = 0;
    int i = from;
    while (i < to) {
	    while (i < to && cov_list[i]) ++i;
	    if (i >= to) break;
	    int j = i + 1;
	    while (j < to && cov_list[j] == 0) ++j;
	    if (j - i > max_cov) {
		    max_cov = j - i;
	    }
	    i = j;
    }
    fprintf(stderr, "max-empty: %d\n", max_cov);
    free(cov_list);
}
