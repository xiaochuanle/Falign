#include "hbn_options.hpp"
#include "build_repeat_reference_regions.hpp"
#include "necat_info.hpp"
#include "../../corelib/arg_parse.hpp"

#include <cstring>

using namespace std;

static const char* kOutFmtNameList[] = {
    "sam",
    "frag-sam",
    "bam",
    "frag-bam",
    "paf"
};

const char* output_format_name(EOutputFmt fmt)
{
    if (fmt >= 0 && fmt < eOutputFmt_Invalid) {
        return kOutFmtNameList[fmt];
    }
    return NULL;
}

EOutputFmt name_to_output_format(const char* name)
{
    EOutputFmt fmt = eOutputFmt_Invalid;
    for (int i = 0; i < eOutputFmt_Invalid; ++i) {
        if (strcmp(name, kOutFmtNameList[i]) == 0) {
            fmt = static_cast<EOutputFmt>(i);
            break;
        }
    }
    return fmt;
}

void
s_dump_enzyme_names_and_seqs(FILE* out)
{
    fprintf(out, "Restriction Enzyme_Name and corresponding Enzyme_Seq:\n");
    const int L = 20;
    print_fixed_width_string(out, "Enzyme_Name", L);
    fprintf(out, "    Enzyme_Seq\n");
    print_fixed_width_string(out, "DpnII", L);
    fprintf(out, "    ^GATC\n");
    print_fixed_width_string(out, "HindIII", L);
    fprintf(out, "    A^AGCTT\n");
    print_fixed_width_string(out, "NcoI", L);
    fprintf(out, "    C^CATGG\n");
    print_fixed_width_string(out, "NlaIII", L);
    fprintf(out, "    CATG^\n");
}

static bool
s_parse_output_format_option(int argc, char* argv[], int& i, const char* arg_name, EOutputFmt& fmt)
{
    if (strcmp(argv[i], arg_name)) return false;
    if (i >= argc) HBN_ERR("Argument to option '%s' is missing", arg_name);
    fmt = name_to_output_format(argv[i+1]);
    if (fmt == eOutputFmt_Invalid) HBN_ERR("Illegal output format '%s'", argv[i+1]);
    i += 2;
    return true;
}

bool HbnProgramOptions::parse(int argc, char* argv[], std::vector<std::string>& query_list)
{
    precheck_args(argc, argv);
    int i = 1;
    while (i < argc) {
        bool r = (argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (parse_data_size_arg_value(argc, argv, i, _query_batch_size, query_batch_size)) continue;
        if (parse_data_size_arg_value(argc, argv, i, _query_upto, query_upto)) continue;

        if (parse_string_arg_value(argc, argv, i, _repeat_bed, repeat_bed)) continue;
        if (parse_bool_arg_value(argc, argv, i, _skip_repeat_bed, skip_repeat_bed)) continue;
        if (parse_int_arg_value(argc, argv, i, _kmer_size, kmer_size)) continue;
        if (parse_int_arg_value(argc, argv, i, _kmer_window, kmer_window)) continue;
        if (parse_int_arg_value(argc, argv, i, _max_kmer_occ, max_kmer_occ)) continue;
        if (parse_real_arg_value(argc, argv, i, _non_repeat_kmer_frac, non_repeat_kmer_frac)) continue;
        if (parse_real_arg_value(argc, argv, i, _repeat_kmer_frac, repeat_kmer_frac)) continue;
        if (parse_real_arg_value(argc, argv, i, _ddf, ddf)) continue;
        if (parse_int_arg_value(argc, argv, i, _kmer_dist, kmer_dist)) continue;
        if (parse_int_arg_value(argc, argv, i, _chain_score, chain_score)) continue;

        if (parse_string_arg_value(argc, argv, i, _enzyme, enzyme)) continue;
        if (parse_int_arg_value(argc, argv, i, _ei_num_hits, ei_num_hits)) continue;
        if (parse_int_arg_value(argc, argv, i, _ei_frag_size, ei_frag_size)) continue;
        if (parse_int_arg_value(argc, argv, i, _ei_end_dist, ei_end_dist)) continue;
        if (parse_int_arg_value(argc, argv, i, _ei_flanking_bases, ei_flanking_bases)) continue;

        if (parse_real_arg_value(argc, argv, i, _perc_identity, perc_identity)) continue;
        if (parse_int_arg_value(argc, argv, i, _max_hsps_per_subject, max_hsps_per_subject)) continue;
        if (parse_int_arg_value(argc, argv, i, _hitlist_size, hitlist_size)) continue;

        if (parse_string_arg_value(argc, argv, i, _output, output)) continue;
        if (s_parse_output_format_option(argc, argv, i, _outfmt, outfmt)) continue;
        if (parse_int_arg_value(argc, argv, i, _num_threads, num_threads)) continue;

        fprintf(stderr, "ERROR: Unrecognised option '%s'\n", argv[i]);
        return false;
    }

    if (i >= argc) return false;
    reference = argv[i];
    ++i;

    for(; i < argc; ++i) query_list.push_back(argv[i]);

    return !query_list.empty();
}

void HbnProgramOptions::simple_usage(int argc, char* argv[])
{
    FILE* out = stderr;
    fprintf(out, "\n");
    fprintf(out, "USAGE\n");

    fprintf(out, "\n");
    fprintf(out, "Building reference repeats:\n");
    fprintf(out, "  %s build-repeat reference repeat-bed\n", argv[0]);

    fprintf(out, "\n");
    fprintf(out, "Mapping sequencing reads to reference:\n");
    fprintf(out, "  %s [OPTIONS] <reference> <read1> [... <readN>]\n", argv[0]);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION\n");
    fprintf(out, "  Alignment toolkit for long noisy chromosome conformation capture (3C) reads\n");

    fprintf(out, "\n");
    fprintf(out, "VERSION\n");
    const char* version = HBN_PACKAGE_VERSION;
    fprintf(out, "  %s\n", version);

    fprintf(out, "\n");
    fprintf(out, "Type option '-help' to see details of optional arguments\n");    
}

void HbnProgramOptions::full_usage(int argc, char* argv[])
{
    string size;
    FILE* out = stderr;
    fprintf(out, "\n");
    fprintf(out, "USAGE\n");
    fprintf(out, "  %s [OPTIONS] <reference> <read1> [... <readN>]\n", argv[0]);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION\n");
    fprintf(out, "  Alignment toolkit for long noisy chromosome conformation capture (3C) reads\n");

    fprintf(out, "\n");
    fprintf(out, "VERSION\n");
    const char* version = HBN_PACKAGE_VERSION;
    fprintf(out, "  %s\n", version);

    fprintf(out, "\n");
    fprintf(out, "OPTIONAL ARGUMENTS\n");

    fprintf(out, "  %s\n", "-h");
    fprintf(out, "    Print USAGE and DESCRIPTION; ignore all other parameters\n");
    fprintf(out, "  %s\n", "-help");
    fprintf(out, "    Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\n");

    fprintf(out, "\n");
    fprintf(out, "  *** Input sequence options\n");
    fprintf(out, "  %s <DataSize>\n", _query_batch_size);
    fprintf(out, "    Batch size of queries (in bp) processed at one time\n");
    size = NStr::UInt8ToString_DataSize(query_batch_size);
    fprintf(out, "    Default = '%s'\n", size.c_str());
    fprintf(out, "  %s <DataSize>\n", _query_upto);
    fprintf(out, "    Map at most this size of queries and skip the rest queries\n");

    fprintf(out, "\n");
    fprintf(out, "  *** Candidate detection options\n");
    fprintf(out, "  %s <path>\n", _repeat_bed);
    fprintf(out, "    A bed file containing reference repeat regions\n");
    fprintf(out, "    If not provided, repeat regions will be created before mapping\n");
    fprintf(out, "  %s\n", _skip_repeat_bed);
    fprintf(out, "    If set, reference repeat regions will not be created and used\n");
    fprintf(out, "  %s <Integer>\n", _kmer_size);
    fprintf(out, "    Kmer size (length of best perfect match)\n");
    fprintf(out, "    Default = '%d'\n", kmer_size);
    fprintf(out, "  %s <Integer>\n", _kmer_window);
    fprintf(out, "    Kmer sampling window size in reference sequences\n");
    fprintf(out, "    Default = '%d'\n", kmer_window);
    fprintf(out, "  %s <Integer>\n", _max_kmer_occ);
    fprintf(out, "    Filter out kmers occuring larger than this value\n");
    fprintf(out, "    Default = '%d'\n", max_kmer_occ);
    fprintf(out, "  %s <Real>\n", _non_repeat_kmer_frac);
    fprintf(out, "    Filter out this fraction of most frequently occurring kmers at non-repeat reference regions\n");
    fprintf(out, "    Default = '%g'\n", non_repeat_kmer_frac);
    fprintf(out, "  %s <Real>\n", _repeat_kmer_frac);
    fprintf(out, "    Filter out this fraction of most frequently occurring kmers at repeat reference regions\n");
    fprintf(out, "    Default = '%g'\n", repeat_kmer_frac);
    fprintf(out, "  %s <Integer>\n", _kmer_dist);
    fprintf(out, "    Read and reference distance between adjacent\n");
    fprintf(out, "    matched kmers in a chain must not greater than this value\n");
    fprintf(out, "    Default = '%d'\n", kmer_dist);
    fprintf(out, "  %s <Real>\n", _ddf);
    fprintf(out, "    DDF factor\n");
    fprintf(out, "    Default = '%g'\n", ddf);
    fprintf(out, "  %s <Integer>\n", _chain_score);
    fprintf(out, "    Minimum score of alignment candidates\n");
    fprintf(out, "    Default = '%d'\n", chain_score);

    fprintf(out, "\n");
    fprintf(out, "  *** Restrict enzyme inference options\n");
    fprintf(out, "  %s <EnzymeSeq>\n", _enzyme);
    fprintf(out, "    The restrict enzyme used for generating Pore-C reads\n");
    fprintf(out, "    If provided, the enzyme inference process will be skipped\n");
    fprintf(out, "    If not provided, a enzyme inference process is performed for each query file\n");
    s_dump_enzyme_names_and_seqs(stderr);
    fprintf(out, "  %s <Integer>\n", _ei_num_hits);
    fprintf(out, "    For each query, consider its INT alignment candidates with highest DDF scores for enzyme inference\n");
    fprintf(out, "    Default = '%d'\n", ei_num_hits);
    fprintf(out, "  %s <Integer>\n", _ei_frag_size);
    fprintf(out, "    Mininum fragment length (bp) considered for enzyme inference\n");
    fprintf(out, "    Default = '%d'\n", ei_frag_size);
    fprintf(out, "  %s <Integer>\n", _ei_end_dist);
    fprintf(out, "    Minimum distance between mapping positions and query ends\n");
    fprintf(out, "    Default = '%d'\n", ei_end_dist);
    fprintf(out, "  %s <Integer>\n", _ei_flanking_bases);
    fprintf(out, "    Looking for restriction enzymes in INT flanking bases of mapping positions\n");
    fprintf(out, "    Default = '%d'\n", ei_flanking_bases);

    fprintf(out, "\n");
    fprintf(out, "  *** Restrict search or results\n");
    fprintf(out, "  %s <Real, 0..100>\n", _perc_identity);
    fprintf(out, "    Minimum percentage of identity of alignments\n");
    fprintf(out, "    Default = '%g'\n", perc_identity);
    fprintf(out, "  %s <Integer>\n", _max_hsps_per_subject);
    fprintf(out, "    Maximum number of alignments per reference sequence\n");
    fprintf(out, "    Deafult = '%d'\n", max_hsps_per_subject);
    fprintf(out, "  %s <Integer>\n", _hitlist_size);
    fprintf(out, "    Maximum number of aligned reference sequences to keep\n");
    fprintf(out, "    Default = '%d'\n", hitlist_size);

    fprintf(out, "\n");
    fprintf(out, "  *** Miscellaneous options\n");
    fprintf(out, "  %s <Integer>\n", _num_threads);
    fprintf(out, "    Number of threads (CPUs) to use in the search\n");
    fprintf(out, "    Default = '%d'\n", num_threads);
    fprintf(out, "  %s <String, Permissible values:", _outfmt);
    for (int i = 0; i < eOutputFmt_Invalid; ++i) fprintf(out, " '%s'", kOutFmtNameList[i]);
    fprintf(out, ">\n");
    fprintf(out, "    Output format\n");
    fprintf(out, "    Default = '%s'\n", kOutFmtNameList[outfmt]);
    fprintf(out, "  %s <file>\n", _output);
    fprintf(out, "    Output file name.\n");
    fprintf(out, "    Default = '%s'\n", output);
}

void HbnProgramOptions::precheck_args(int argc, char* argv[])
{
    for (int i = 1; i < argc; ++i) {
        bool r = (argv[i][0] != '-') || (argv[i][0] == '-' && strlen(argv[i]) == 1);
        if (r) break;

        if (strcmp(argv[i] + 1, "h") == 0) {
            simple_usage(argc, argv);
            exit(0);
        }
        if (strcmp(argv[i] + 1, "help") == 0) {
            full_usage(argc, argv);
            exit(0);
        }
    }
}