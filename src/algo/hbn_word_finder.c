#include "hbn_word_finder.h"

#include "../corelib/cstr_util.h"
#include "../corelib/ksort.h"

#define ddfk_hash_lt(__a, __b) ( ((__a).hash < (__b).hash) || ((__a).hash == (__b).hash && (__a).context < (__b).context) )
KSORT_INIT(ddfk_hash_lt, DDFKmer, ddfk_hash_lt);
void
sort_ddfk_hash_lt(size_t n, DDFKmer* a)
{
    ks_introsort_ddfk_hash_lt(n, a);
}

#define ddfkm_context_lt(a, b) ((a).context < (b).context)
KSORT_INIT(ddfkm_context_lt, DDFKmerMatch, ddfkm_context_lt);
void
sort_ddfkm_context_lt(size_t n, DDFKmerMatch* a)
{
    ks_introsort_ddfkm_context_lt(n, a);
}

#define ddfkm_soff_lt(a, b) ((a).soff < (b).soff)
KSORT_INIT(ddfkm_soff_lt, DDFKmerMatch, ddfkm_soff_lt);
void
sort_ddfkm_soff_lt(size_t n, DDFKmerMatch* a)
{
    ks_introsort_ddfkm_soff_lt(n, a);
}

#define u64_hash_lt(__a, __b) ((__a) < (__b))
KSORT_INIT(u64_hash_lt, u64, u64_hash_lt);

WordFinderThreadData*
WordFinderThreadDataNew(SeqReader* db,
    HbnLookupTable* lktbl,
    const int kmer_size,
    const int kmer_window,
    const int min_ddf_score,
    const double max_ddf)
{
    WordFinderThreadData* data = (WordFinderThreadData*)calloc(1, sizeof(WordFinderThreadData));
    data->db = db;
    data->lktbl = lktbl;
    kv_init(data->ddfk_list);
    kv_init(data->ddfkm_list);
    kv_init(data->skmi_list);
    data->kmer_size = kmer_size;
    data->kmer_window = kmer_window;
    kv_init(data->hash_list);
    return data;
}

WordFinderThreadData*
WordFinderThreadDataFree(WordFinderThreadData* data)
{
    kv_destroy(data->ddfk_list);
    kv_destroy(data->ddfkm_list);
    kv_destroy(data->skmi_list);
    kv_destroy(data->hash_list);
    free(data);
    return NULL;
}

void 
WordFinderThreadDataClear(WordFinderThreadData* data)
{
    kv_clear(data->ddfk_list);
    kv_clear(data->ddfkm_list);
    kv_clear(data->skmi_list);
}

static int
extract_hash_values(const u8* read,
    const int read_size,
    const int kmer_size,
    const int window_size,
    vec_u64* hash_list)
{
    if (read_size < kmer_size) return 0;
    const int intersect = kmer_size > window_size;
	u64 intersect_mask = 0;
	int stride = kmer_size - window_size;
	if (intersect) intersect_mask = (U64_ONE << (stride << 1)) - 1;

	if (!intersect) {
		for (size_t j = 0; j <= read_size - kmer_size; j += window_size) {
			u64 hash = 0;
			for (int k = 0; k < kmer_size; ++k) {
				size_t pos = j + k;
				u8 c = read[pos];
                hbn_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
            kv_push(u64, *hash_list, hash);
		}
	} else {
		u64 hash = 0;
		for (int j = 0; j < kmer_size; ++j) {
			size_t pos = j;
			u8 c = read[pos];
            hbn_assert(c >= 0 && c < 4);
			hash = (hash << 2) | c;
		}
		kv_push(u64, *hash_list, hash);
		for (u64 j = window_size; j <= read_size - kmer_size; j += window_size) {
			hash &= intersect_mask;
			for (int k = stride; k < kmer_size; ++k) {
				size_t pos = j + k;
                hbn_assert(pos < read_size, "p = %d, read_size = %d, j = %d, k = %d, stride = %d", pos, read_size, j, k, stride);
				u8 c = read[pos];
                hbn_assert(c >= 0 && c < 4);
				hash = (hash << 2) | c;
			}
			kv_push(u64, *hash_list, hash);
		}
	}
    return kv_size(*hash_list);
}

static int
collect_ddfkmer_subseq(const int context,
    const u8* read,
    const int read_from,
    const int read_to,
    const int kmer_size,
    const int window_size,
    vec_u64* hash_list,
    vec_ddfk* ddfk_list)
{
    int n = read_to - read_from;
    const int SL = n, SR = n;
    int s = 0;
    DDFKmer ddfk;
    int cnt = 0;
    ddfk.context = context;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1);

    while (s < n) {
        int e = s + SL;
        e = hbn_min(e, n);
        kv_clear(*hash_list);
        int n_kmer = extract_hash_values(read + read_from + s, e - s, kmer_size, window_size, hash_list);
        for (int i = 0; i < n_kmer; ++i) {
            ddfk.offset = read_from + s + i * window_size;
            ddfk.hash = kv_A(*hash_list, i);
            if (ddfk.hash == kMaxHashValue) continue;
            kv_push(DDFKmer, *ddfk_list, ddfk);
            ++cnt;
        }
        s = e + SR;                
    }
    return cnt;
}

static int
collect_ddfkmer_one_query(const int context,
    const u8* read,
    const int read_size,
    const int kmer_size,
    const int kmer_window,
    vec_u64* hash_list,
    vec_ddfk* ddfk_list)
{    
    return collect_ddfkmer_subseq(context, read, 0, read_size, kmer_size, kmer_window, hash_list, ddfk_list);
}

static void
extract_kmer_matches(HbnLookupTable* lktbl, DDFKmer* ddfka, int ddfkc, vec_ddfkm* ddfkm_list)
{
    for (int i = 0; i < ddfkc; ++i) {
        u64 oc = 0;
        u64* ol = HbnLookupTable_ExtractKmerOffsetList(lktbl, ddfka[i].hash, &oc);
        DDFKmerMatch ddfkm;
        for (u64 x = 0; x < oc; ++x) {
            ddfkm.context = ddfka[i].context;
            ddfkm.qoff = ddfka[i].offset;
            ddfkm.soff = ol[x];
            kv_push(DDFKmerMatch, *ddfkm_list, ddfkm);
        }        
    }
}

static BOOL
s_remove_repetative_ddfkm(DDFKmerMatch* ddfkma, int ddfkmc, int min_cnt, SubjectKmInfo* ski)
{
#if 1
	const int kMaxKmerOcc = 20;
    int i = 0;
    while (i < ddfkmc) {
        int j = i + 1;
        while (j < ddfkmc && ddfkma[i].soff == ddfkma[j].soff) ++j;
        if (j - i > kMaxKmerOcc) {
            int cnta[2] = { 0, 0 };
            for (int k = i; k < j; ++k) ++cnta[ddfkma[k].context&1];
            if (cnta[0] > kMaxKmerOcc) {
                for (int k = i; k < j; ++k) {
                    if (!(ddfkma[k].context&1)) ddfkma[k].context=-1;
                }
            } 
            if (cnta[1] > kMaxKmerOcc) {
                for (int k = i; k < j; ++k) {
                    if (ddfkma[k].context&1) ddfkma[k].context = -1;
                }
            }
        }
        i = j;
    }
    i = 0;
    for (int k = 0; k < ddfkmc; ++k) {
        if (ddfkma[k].context >= 0) ddfkma[i++] = ddfkma[k];
    }
    ddfkmc = i;
#endif

    int cnta[2] = { 0, 0 };
    for (int k = 0; k < ddfkmc; ++k) {
        ++cnta[ddfkma[k].context&1];
    }
    ski->f = hbn_max(cnta[0], cnta[1]);
    if (ski->f < min_cnt) return FALSE;
    ski->km_cnt = ddfkmc;
    return TRUE;
}

static void
s_build_skmi_list(HbnLookupTable* lktbl, DDFKmerMatch* ddfkma, int ddfkmc, const int min_cnt, vec_skmi* skmi_list)
{
    sort_ddfkm_soff_lt(ddfkmc, ddfkma);
    int i = 0;
    while (i < ddfkmc) {
        int sid = 0;
        size_t soff = 0, send = 0;
        HbnLookupTable_KmerOffset2SeqInfo(lktbl, ddfkma[i].soff, &sid, &soff, &send);
        int j = i + 1;
        while (j < ddfkmc && ddfkma[j].soff < send) ++j;
        for (int k = i; k < j; ++k) ddfkma[k].soff -= soff;
        if (j - i >= min_cnt) {
            SubjectKmInfo ski = { sid, i, j - i, 0 };
            if (s_remove_repetative_ddfkm(ddfkma + i, j - i, min_cnt, &ski)) {
                kv_push(SubjectKmInfo, *skmi_list, ski);
            }
        }
        i = j;
    }
}

void
hbn_ddfkm(WordFinderThreadData* data, int read_id)
{
    WordFinderThreadDataClear(data);

    int ctx_id = read_id * 2;
    const u8* s = SeqReader_Seq(data->db, read_id, FWD);
    const int sl = SeqReader_SeqSize(data->db, read_id);
    collect_ddfkmer_one_query(ctx_id, s, sl, data->kmer_size, data->kmer_window, &data->hash_list, &data->ddfk_list);

    ++ctx_id;
    s = SeqReader_Seq(data->db, read_id, REV);
    collect_ddfkmer_one_query(ctx_id, s, sl, data->kmer_size, data->kmer_window, &data->hash_list, &data->ddfk_list);

    extract_kmer_matches(data->lktbl, kv_data(data->ddfk_list), kv_size(data->ddfk_list), &data->ddfkm_list);
    s_build_skmi_list(data->lktbl, kv_data(data->ddfkm_list), kv_size(data->ddfkm_list), data->min_kmer_cnt, &data->skmi_list);
}
