#ifndef __SMALL_EDLIB_ALIGN_HPP
#define __SMALL_EDLIB_ALIGN_HPP

#include "../corelib/hbn_aux.h"

#include <string>
#include <vector>

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
        std::string& qaln,
        std::string& saln);

int
small_edlib_nw_dist(const u8* query,
		const int query_size,
		const u8* target,
		const int target_size,
		small_edlib_align_struct* _align_data);

#endif // __SMALL_EDLIB_ALIGN_HPP