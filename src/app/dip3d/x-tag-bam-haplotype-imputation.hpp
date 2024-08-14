#ifndef __X_TAG_BAM_HAPLOTYPE_IMPUTATION_HPP
#define __X_TAG_BAM_HAPLOTYPE_IMPUTATION_HPP

#include "frag-hap.hpp"
#include "x-tag-bam-options.hpp"

void haplotype_imputation_mt(BamTagOptions* options, FragHapInfo* fhia, size_t fhic);

#endif // __X_TAG_BAM_HAPLOTYPE_IMPUTATION_HPP