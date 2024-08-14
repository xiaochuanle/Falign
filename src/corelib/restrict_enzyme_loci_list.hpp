#ifndef __RESTRICT_ENZYME_LOCI_LIST_H
#define __RESTRICT_ENZYME_LOCI_LIST_H

#include "hbn_aux.h"
#include "unpacked_seqdb.hpp"

#include <vector>

#define MAX_ENZYME_SIZE 6

#define LEFT_ENZYME_LOCI_MATCH 1
#define RIGHT_ENZYME_LOCI_MATCH 2

typedef struct {
    char enzyme[MAX_ENZYME_SIZE];
    u8 encoded_enzyme[MAX_ENZYME_SIZE];
    int enzyme_size;
    int break_loci;
    u64 enzyme_hash;
    u64 enzyme_mask;
} RestrictEnzyme;

void
RestrictEnzyme_Init(const char* enzyme, RestrictEnzyme* re);

typedef struct {
    size_t enzyme_loci_offset;
    int enzyme_loci_cnt;
} SeqRestrictEnzymeLociInfo;

typedef struct {
    RestrictEnzyme enzyme;
    SeqRestrictEnzymeLociInfo* seq_reloci_info_array;
    int* reloci_array;
    size_t reloci_cnt;
} RestrictEnzymeLociList;

RestrictEnzymeLociList*
RestrictEnzymeLociListFree(RestrictEnzymeLociList* list);

RestrictEnzymeLociList*
RestrictEnzymeLociListNew(HbnUnpackedDatabase* subjects, const char* enzyme);

int 
offset_to_enzyme_intv_idx(const int* loci_array, const int loci_cnt, const int offset, int* intv_cnt);

struct QueryVdfEndPointList {
    RestrictEnzyme* enzyme;
    const u8* fwd_query;
    const u8* rev_query;
    int query_size;
    std::vector<int> fwd_vdf_endpoint_list;
    std::vector<int> rev_vdf_endpoint_list;
};

void 
QueryVdfEndPointList_Setup(QueryVdfEndPointList* qvep_list,
    RestrictEnzyme* enzyme,
    const u8* fwd_query,
    const u8* rev_query,
    const int query_size);

#endif // __RESTRICT_ENZYME_LOCI_LIST_H