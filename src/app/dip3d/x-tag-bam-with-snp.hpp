#ifndef __X_TAG_BAM_WITH_SNP_HPP
#define __X_TAG_BAM_WITH_SNP_HPP

#include "../../corelib/seq_name2id_map.hpp"
#include "../../corelib/unpacked_seqdb.hpp"
#include "bam-writer.hpp"
#include "frag-hap.hpp"
#include "vcf-reader.hpp"
#include "x-tag-bam-options.hpp"

#include <vector>

void snp_tag_bam_mt(BamTagOptions* options, std::vector<FragHapInfo>& hap_list);

#endif // __X_TAG_BAM_WITH_SNP_HPP