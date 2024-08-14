#include "../../corelib/arg_parse.hpp"
#include "../../corelib/frag_id.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/pdqsort.h"
#include "../../htslib/hts.h"
#include "../../htslib/sam.h"
#include "../../htslib/vcf.h"
#include "../../sw/hbn_traceback_aux.h"
#include "bam-writer.hpp"
#include "frag-hap.hpp"
#include "sam_map_info.hpp"
#include "vcf-reader.hpp"
#include "x-tag-bam-haplotype-imputation.hpp"
#include "x-tag-bam-options.hpp"
#include "x-tag-bam-with-snp.hpp"
#include "x-tag-bam.hpp"

#include <cstdlib>

#include <iostream>
#include <mutex>
#include <set>

using namespace std;

int tag_bam_main(int argc, char* argv[])
{
    BamTagOptions options;
    if (!options.parse_arguments(argc, argv)) {
        options.dump_usage(argc, argv);
        exit (1);
    }
    options.dump();
    create_directory(options.wrk_dir);

    vector<FragHapInfo> snp_tagged_frag_hap_list;
    snp_tag_bam_mt(&options, snp_tagged_frag_hap_list);
    pdqsort(snp_tagged_frag_hap_list.begin(), snp_tagged_frag_hap_list.end(), 
        [](const FragHapInfo& x, const FragHapInfo& y) { return x.read_id < y.read_id; });

    HBN_LOG("SNP-tagged frag and contact stats for %s", options.bam_path);
    frag_and_contact_stats_hap_list(snp_tagged_frag_hap_list.data(), snp_tagged_frag_hap_list.size());
    
    vector<FragHapInfo> imputed_frag_hap_list = snp_tagged_frag_hap_list;
    haplotype_imputation_mt(&options, imputed_frag_hap_list.data(), imputed_frag_hap_list.size());

    HBN_LOG("Imputed frag and contact stats for %s", options.bam_path);
    frag_and_contact_stats_hap_list(imputed_frag_hap_list.data(), imputed_frag_hap_list.size());

    char path[HBN_MAX_PATH_LEN];
    snprintf(path, HBN_MAX_PATH_LEN, "%s/snp-tagged-frag-hap-list", options.wrk_dir);
    save_frag_hap_info_list(snp_tagged_frag_hap_list.data(), snp_tagged_frag_hap_list.size(), path);

    snprintf(path, HBN_MAX_PATH_LEN, "%s/imputed-frag-hap-list", options.wrk_dir);
    save_frag_hap_info_list(imputed_frag_hap_list.data(), imputed_frag_hap_list.size(), path);

    snprintf(path, HBN_MAX_PATH_LEN, "%s/tagged.bam", options.wrk_dir);
    tag_bam_st(imputed_frag_hap_list.data(), imputed_frag_hap_list.size(), options.num_threads, options.bam_path, path);

    return 0;
}
