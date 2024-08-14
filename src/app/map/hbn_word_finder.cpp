#include "hbn_word_finder.hpp"

#include "../../corelib/pdqsort.h"
#include "window_masker.hpp"

using namespace std;

HbnWordFinder::HbnWordFinder(HbnUnpackedDatabase& db,
    const char* repeat_bed_path,
    const bool skip_repeat_regions,
    int num_threads,
    int kmer_size,
    int kmer_window,
    int max_kmer_occ,
    double repeat_frac,
    double non_repeat_frac)
{
    vector<ChrIntv> repeat_intv_list;
    if (repeat_bed_path) {
        SeqName2IdMap refname2id;
        const int num_chr = db.NumSeqs();
        for (int i = 0; i < num_chr; ++i) refname2id.add_one_name(db.SeqName(i));
        refname2id.build_name2id_map();
        load_bed_list(repeat_bed_path, refname2id, repeat_intv_list);
    } else if (!skip_repeat_regions) {
        fprintf(stderr, "Repeat bed file is not provided. We will create repeat reference regions first, which will take a fiew minutes.\n");
        compute_repeat_regions(&db, 1, repeat_intv_list);
    }
    if (!repeat_intv_list.empty()) {
        HBN_LOG("Repeat region stats:");
        repeat_intv_stats(repeat_intv_list.data(), repeat_intv_list.size(), &db);
    }
    vector<ChrIntv> non_repeat_intv_list;
    compute_complement_intv_list(repeat_intv_list.data(), repeat_intv_list.size(), db, non_repeat_intv_list);
    if (!non_repeat_intv_list.empty()) {
        HBN_LOG("Non-repeat region stats:");
        repeat_intv_stats(non_repeat_intv_list.data(), non_repeat_intv_list.size(), &db);
    }

    M_kmer_size = kmer_size;
    M_lktbl = new HbnLookupTable(db,
        repeat_intv_list.data(),
        repeat_intv_list.size(),
        non_repeat_intv_list.data(),
        non_repeat_intv_list.size(),
        kmer_size,
        kmer_window,
        max_kmer_occ,
        repeat_frac,
        non_repeat_frac,
        num_threads);
}

static void
s_extract_kmer_matches(HbnLookupTable* lktbl,
    int kmer_size,
    const int query_id,
    const int query_strand,
    const u8* query,
    const int query_size,
    mt19937* gen,
    uniform_int_distribution<>* dist,
    vector<KmerMatch>& km_list)
{
    if (query_size < kmer_size) return;
    const u64 kMaxHashValue = U64_ONE << (kmer_size << 1); 
    const u64 kHashMask = kMaxHashValue - 1;

    u64 hash = 0;
    for (int i = 0; i < kmer_size - 1; ++i) {
        u8 c = query[i];
        hash = (hash << 2) | c;
        hash &= kHashMask;
    }

    KmerInfo* ol = nullptr;
    u64 cnt = 0;
    KmerMatch km;
    km.qdir = query_strand;
    for (int i = 0, p = kmer_size - 1; p < query_size; ++i, ++p) {
        u8 c = query[p];
        hash = (hash << 2) | c;
        hash &= kHashMask;
        //if (dist && (*dist)(*gen)) continue;
        lktbl->extract_offset_list(hash, &ol, &cnt);
        km.qoff = i;
        for (u64 k = 0; k < cnt; ++k) {
            km.sid = ol[k].seq_id;
            km.soff = ol[k].seq_offset;
            km_list.push_back(km);
        }
    }
}

void HbnWordFinder::extract_kmer_matches(const int query_id,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size,
    std::vector<KmerMatch>& km_list,
    std::mt19937* gen,
    std::uniform_int_distribution<>* dist)
{
    km_list.clear();

    if (fwd_query) {
        s_extract_kmer_matches(M_lktbl,
            M_kmer_size,
            query_id,
            FWD,
            fwd_query,
            query_size,
            gen,
            dist,
            km_list);
    }

    if (fwd_query) {
        s_extract_kmer_matches(M_lktbl,
            M_kmer_size,
            query_id,
            REV,
            rev_query,
            query_size,
            gen,
            dist,
            km_list);
    }

    pdqsort(km_list.begin(), km_list.end(), [](const KmerMatch& x, const KmerMatch& y) {
        return (x.sid < y.sid) || (x.sid == y.sid && x.soff < y.soff);
    });
    KmerMatch* kma = km_list.data();
    int kmc = km_list.size();
    int i = 0;
    while (i < kmc) {
        int j = i + 1;
        while (j < kmc && kma[i].sid == kma[j].sid && kma[i].soff == kma[j].soff) ++j;
        if (j - i > 20) for (int k = i; k < j; ++k) kma[k].qoff = -1;
        i = j;
    }
    int n = 0;
    for (i = 0; i < kmc; ++i) if (kma[i].qoff >= 0) kma[n++] = kma[i];
    km_list.resize(n);
}
