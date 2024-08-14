#ifndef __HBN_OPTIONS_HPP
#define __HBN_OPTIONS_HPP

#include <stdlib.h>
#include "../../corelib/hbn_aux.h"
#include "../../corelib/read_list.hpp"

enum EOutputFmt {
    eOutputFmt_SAM = 0,
    eOutputFmt_FragSAM,
    eOutputFmt_BAM,
    eOutputFmt_FragBAM,
    eOutputFmt_PAF,
    eOutputFmt_Invalid
} ;

const char* output_format_name(EOutputFmt fmt);

EOutputFmt name_to_output_format(const char* name);

struct HbnProgramOptions {
    /// input sequence options
    size_t          query_batch_size { 4000000000 };
    size_t          query_upto { U64_MAX };

    const char*     _query_batch_size { "-query_batch_size" };
    const char*     _query_upto { "-query_upto" };

    /// ddf scoring options on detecting candidate subject subsequences

    const char*     repeat_bed { nullptr };
    bool            skip_repeat_bed { false };
    int             kmer_size { 15 };
    int             kmer_window { 4 };
    int             max_kmer_occ { 100 };
    double          non_repeat_kmer_frac { 0.0001 };
    double          repeat_kmer_frac { 0.1 };

    double          ddf { 0.2 };
    int             kmer_dist { 800 };
    int             chain_score { 3 };

    const char*     _repeat_bed { "-repeat_bed" };
    const char*     _skip_repeat_bed { "-skip_repeat_bed" };
    const char*     _kmer_size { "-kmer_size" };
    const char*     _kmer_window { "-kmer_window" };
    const char*     _max_kmer_occ { "-max_kmer_occ" };
    const char*     _non_repeat_kmer_frac { "-non_repeat_kmer_frac" };
    const char*     _repeat_kmer_frac { "-repeat_kmer_frac" };
    const char*     _ddf { "-ddf" };
    const char*     _kmer_dist { "-kmer_dist" };
    const char*     _chain_score { "-chain_score" };

    /// enzyme
    const char*     enzyme { nullptr };
    int             ei_num_hits { 3 };
    int             ei_frag_size { 200 };
    int             ei_end_dist { 200 };
    int             ei_flanking_bases { 30 };

    const char*     _enzyme { "-enzyme" };
    const char*     _ei_num_hits { "-ei_num_hits" };
    const char*     _ei_frag_size { "-ei_frag_size" };
    const char*     _ei_end_dist { "-ei_end_dist" };
    const char*     _ei_flanking_bases { "-ei_flanking_bases" };

    /// restrict search or results
    double          perc_identity { 75.0 };
    int             max_hsps_per_subject { 50 };
    int             hitlist_size { 3 };

    const char*     _perc_identity { "-perc_identity" };
    const char*     _max_hsps_per_subject { "-chr_maps" };
    const char*     _hitlist_size { "-target_chrs" };

    /// Miscellaneous options

    int             num_threads { 1 };
    EOutputFmt      outfmt { eOutputFmt_PAF };

    const char*     reference { nullptr };
    const char*     output { "-" };

    const char*     _num_threads { "-num_threads" };
    const char*     _outfmt { "-outfmt" };
    const char*     _output { "-out" };

    /////////////////////////

    void precheck_args(int argc, char* argv[]);

    void simple_usage(int argc, char* argv[]);

    void full_usage(int argc, char* argv[]);

    bool parse(int argc, char* argv[], std::vector<std::string>& query_list);
};

void
s_dump_enzyme_names_and_seqs(FILE* out);

#endif // __HBN_OPTIONS_HPP