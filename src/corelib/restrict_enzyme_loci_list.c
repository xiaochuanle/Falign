#include "restrict_enzyme_loci_list.h"

void
RestrictEnzyme_Init(const char* enzyme, RestrictEnzyme* re)
{
    int l = strlen(enzyme);
    //if (l > MAX_ENZYME_SIZE) HBN_ERR("Restrict enzyme is too long: %s", enzyme);

    re->enzyme_size = 0;
    re->break_loci = -1;
    for (int i = 0; i < l; ++i) {
        int c = enzyme[i];
        if (c == '^') {
            re->break_loci = i;
            continue;
        }
        c = nst_nt16_table[c];
        if (c > 3) {
            HBN_ERR("Illegal character '%c' in restrict enzyme '%s'", enzyme[i], enzyme);
        }
        re->enzyme[re->enzyme_size] = enzyme[i];
        re->encoded_enzyme[re->enzyme_size] = c;
        ++re->enzyme_size;
    }
    if (re->break_loci == -1) {
        HBN_ERR("'^' is missing from restrict enzyme '%s'", enzyme);
    }
    re->enzyme_mask = (U64_ONE << (2 * re->enzyme_size)) - 1;
    re->enzyme_hash = 0;
    for (int i = 0; i < re->enzyme_size; ++i) {
        re->enzyme_hash = (re->enzyme_hash << 2) | re->encoded_enzyme[i];
    }

    for (int i = 0; i < re->enzyme_size; ++i) {
        int fc = re->encoded_enzyme[i];
        int rc = re->encoded_enzyme[re->enzyme_size - 1 - i];
        if (fc + rc != 3) {
            HBN_ERR("Illegal restrict enzyme '%s'", enzyme);
        }
    }
}

static void
s_extract_enzyme_loci_list_for_one_subject(const u8* subject,
    const int subject_size,
    RestrictEnzyme* re,
    vec_int* reloci_list)
{
    const u64 kHashMask = re->enzyme_mask;
    const u64 kEnzymeHash = re->enzyme_hash;
    const int enzyme_size = re->enzyme_size;
    kv_clear(*reloci_list);
    kv_push(int, *reloci_list, 0);
    u64 hash = 0;
    for (int i = 1; i < enzyme_size && i < subject_size; ++i) hash = (hash << 2) | subject[i];
    for (int i = enzyme_size; i < subject_size - 1; ++i) {
        hash = (hash << 2) | subject[i];
        hash = hash & kHashMask;
        if (hash == kEnzymeHash) {
            int offset = i + 1 - enzyme_size;
            kv_push(int, *reloci_list, offset);
        }
    }
    kv_push(int, *reloci_list, subject_size);
    kv_push(int, *reloci_list, subject_size + 1);
}

RestrictEnzymeLociList*
RestrictEnzymeLociListNew(SeqReader* subjects, const char* enzyme)
{
    RestrictEnzymeLociList* list = (RestrictEnzymeLociList*)calloc(1, sizeof(RestrictEnzymeLociList));
    RestrictEnzyme_Init(enzyme, &list->enzyme);
    const int n_seq = SeqReader_NumSeqs(subjects);
    list->seq_reloci_info_array = (SeqRestrictEnzymeLociInfo*)calloc(n_seq, sizeof(SeqRestrictEnzymeLociInfo));
    kv_dinit(vec_int, reloci_list);
    size_t n_reloci = 0;
    for (int i = 0; i < n_seq; ++i) {
        const u8* subject = SeqReader_Seq(subjects, i, FWD);
        const int subject_size = SeqReader_SeqSize(subjects, i);
        s_extract_enzyme_loci_list_for_one_subject(subject, subject_size, &list->enzyme, &reloci_list);
        int n = kv_size(reloci_list);
        list->seq_reloci_info_array[i].enzyme_loci_offset = n_reloci;
        list->seq_reloci_info_array[i].enzyme_loci_cnt = n;
        n_reloci += n;
    }
    list->reloci_array = (int*)calloc(n_reloci, sizeof(int));
    size_t i_reloci = 0;
    for (int i = 0; i < n_seq; ++i) {
        const u8* subject = SeqReader_Seq(subjects, i, FWD);
        const int subject_size = SeqReader_SeqSize(subjects, i);
        s_extract_enzyme_loci_list_for_one_subject(subject, subject_size, &list->enzyme, &reloci_list);
        for (size_t k = 0; k < kv_size(reloci_list); ++k) {
            list->reloci_array[i_reloci++] = kv_A(reloci_list, k);
        }
    }    
    kv_destroy(reloci_list);
    hbn_assert(i_reloci == n_reloci);
    list->reloci_cnt = n_reloci;

    return list;
}

RestrictEnzymeLociList*
RestrictEnzymeLociListFree(RestrictEnzymeLociList* list)
{
    free(list->seq_reloci_info_array);
    free(list->reloci_array);
    free(list);
    return NULL;
}

int 
offset_to_enzyme_intv_idx(const int* loci_array, const int loci_cnt, const int offset, int* intv_cnt)
{
    hbn_assert(loci_cnt >= 3);
    int left = 0, right = loci_cnt, mid = 0;
    while (left < right) {
        mid = (left + right) / 2;
        if (offset >= loci_array[mid]) {
            if (mid == loci_cnt - 1) break;
            if (offset < loci_array[mid+1]) break;
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    if (mid == loci_cnt - 1 || mid == loci_cnt - 2) {
        mid = loci_cnt - 3;
    }
    hbn_assert(mid >= 0);
    hbn_assert(offset >= loci_array[mid] && offset <= loci_array[mid+1]);

    if (intv_cnt) *intv_cnt = loci_cnt - 1;
    return mid;
}

void QueryVdfEndPointList_Init(QueryVdfEndPointList* qvep_list)
{
    qvep_list->fwd_query = NULL;
    qvep_list->rev_query = NULL;
    qvep_list->query_size = 0;
    kv_init(qvep_list->fwd_vdf_endpoint_list);
    kv_init(qvep_list->rev_vdf_endpoint_list);
    kv_init(qvep_list->perfect_match_offset_point_list);
    kv_init(qvep_list->perfect_match_end_point_list);
}

void QueryVdfEndPointList_Destroy(QueryVdfEndPointList* qvep_list)
{
    kv_destroy(qvep_list->fwd_vdf_endpoint_list);
    kv_destroy(qvep_list->rev_vdf_endpoint_list);
    kv_destroy(qvep_list->perfect_match_offset_point_list);
    kv_destroy(qvep_list->perfect_match_end_point_list);
}

void 
QueryVdfEndPointList_Setup(QueryVdfEndPointList* qvep_list,
    RestrictEnzyme* enzyme,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size)
{
    qvep_list->enzyme = enzyme;
    qvep_list->fwd_query = fwd_query;
    qvep_list->rev_query = rev_query;
    qvep_list->query_size = query_size;
    kv_clear(qvep_list->fwd_vdf_endpoint_list);
    kv_clear(qvep_list->rev_vdf_endpoint_list);
    s_extract_enzyme_loci_list_for_one_subject(fwd_query, query_size,
        enzyme, &qvep_list->fwd_vdf_endpoint_list);

    s_extract_enzyme_loci_list_for_one_subject(rev_query, query_size,
        enzyme, &qvep_list->rev_vdf_endpoint_list);

#if 0
    HBN_LOG("fwd vdf point %d:", query_size);
	int n_v = kv_size(qvep_list->fwd_vdf_endpoint_list);
	int* a = kv_data(qvep_list->fwd_vdf_endpoint_list);
    for (size_t i = 0; i < n_v; ++i) {
        fprintf(stderr, "%zu\t%d\t", i, kv_A(qvep_list->fwd_vdf_endpoint_list, i));
	if (i > 0 && a[i] < query_size) for (int j = 0; j < enzyme->enzyme_size; ++j) { int c = fwd_query[a[i]+j]; c = "ACGT"[c]; fprintf(stderr, "%c", c); }
	fprintf(stderr, "\n");
    }

    HBN_LOG("rev vdf point %d:", query_size);
    for (size_t i = 0; i < kv_size(qvep_list->rev_vdf_endpoint_list); ++i) {
        fprintf(stderr, "%zu\t%d\n", i, kv_A(qvep_list->rev_vdf_endpoint_list, i));
    }
#endif

    hbn_assert(kv_size(qvep_list->fwd_vdf_endpoint_list) == kv_size(qvep_list->rev_vdf_endpoint_list));
}
