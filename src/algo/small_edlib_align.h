#ifndef __SMALL_EDLIB_ALIGN_H
#define __SMALL_EDLIB_ALIGN_H

#include "../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

#define SMALL_EDLIB_MAX_SEQ_SIZE 4096

typedef void* small_edlib_align_struct;

small_edlib_align_struct*
small_edlib_align_struct_new();

small_edlib_align_struct*
small_edlib_align_struct_free(small_edlib_align_struct* data);

int
small_edlib_nw(const u8* query,
		const int query_size,
		const u8* target,
		const int target_size,
		small_edlib_align_struct* _align_data,
        kstring_t* qaln,
        kstring_t* saln);

int
small_edlib_nw_1(const u8* query,
		const int query_size,
		const u8* target,
		const int target_size,
        const int max_dist,
		small_edlib_align_struct* _align_data);

#ifdef __cplusplus
}
#endif

#endif // __SMALL_EDLIB_ALIGN_H