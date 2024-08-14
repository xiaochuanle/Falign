#ifndef __X_TAG_BAM_HPP
#define __X_TAG_BAM_HPP

#include "frag-hap.hpp"

void tag_bam_mt(FragHapInfo* fhia, size_t fhic, const int num_threads, 
    const char* input_bam_path, const char* tagged_bam_path);

void tag_bam_st(FragHapInfo* fhia, size_t fhic, const int num_threads, 
    const char* input_bam_path, const char* tagged_bam_path);

#endif // __X_TAG_BAM_HPP