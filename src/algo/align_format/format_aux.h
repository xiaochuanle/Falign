#ifndef __FORMAT_AUX_H
#define __FORMAT_AUX_H

#include "../../corelib/hbn_aux.h"

#ifdef __cplusplus
extern "C" {
#endif

void dump_cigar_string(const char* qaln, const char* saln, const int aln_size, kstring_t* out);

void dump_md_string(const char* qaln, const char* saln, const int aln_size, kstring_t* out);

#ifdef __cplusplus
}
#endif

#endif // __FORMAT_AUX_H