#include "hbn_options.h"

#include "../../corelib/cstr_util.h"
#include "../../corelib/hbn_package_version.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

using namespace std;

static const char* kOutFmtNameList[] = {
    "sam",
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

const char* kDefault_QueryBatchSize = "1g";
const char* kOpt_QueryBatchSize = "query_batch_size";

const char* kOpt_MapUpto = "query_up_to";

const int kMinKmerSize = 8;
const int kMaxKmerSize = 32;
const int kDefault_KmerSize = 17;
const char* kOpt_KmerSize = "kmer_size";

const int kMinKmerWindow = 1;
const int kDefault_KmerWindow = 5;
const char* kOpt_KmerWindow = "kmer_window";

const double kMinRepFrac = 0.0;
const double kMaxRepFrac = 1.0;
const double kDefault_RepFrac = 0.001;
const char* kOpt_RepFrac = "rep_frac";

const int kMinMaxKmerOcc = 1;
const int kDefault_MaxKmerOcc = 200;
const char* kOpt_MaxKmerOcc = "max_kmer_occ";

const double kMinDDF = 0.0;
const double kDefault_DDF = 0.20;
const char* kOpt_DDF = "ddf";

const int kMinKmerDist = 1;
const int kDefaultKmerDist = 800;
const char* kOpt_KmerDist = "kmer_dist";

const int kMinChainScore = 0;
const int kDefaultChainScore = 2;
const char* kOpt_ChainScore = "chain_score";

const double kMinPercIdentity = 0.0;
const double kMaxPercIdentity = 100.0;
const double kDefault_PercIdentity = 75.0;
const char* kOpt_Perc_Identity = "perc_identity";

const int kMinMaxAlignsPercSubject = 1;
const int kDefault_MaxAlignsPerSubject = 50;
const char* kOpt_MaxAlignsPerSubject = "aligns_per_subject";

const int kMinMaxTargetSubjects = 1;
const int kDefault_MaxTargetSubjects = 3;
const char* kOpt_MaxTargetSubjects = "target_subjects";

const int kMinNumThreads = 1;
const int kDefault_NumThreads = 1;
const char* kOpt_NumThreads = "num_threads";

const EOutputFmt kDefault_OutFmt = eOutputFmt_SAM;
const char* kOpt_OutFmt = "outfmt";

const int kDefault_DumpByFile = 0;
const char* kOpt_DumpByFile = "dump_by_file";

const char* kDefault_Output = "-";
const char* kOpt_Output = "out";

const char* kDefault_OutputDir = NULL;
const char* kOpt_OutputDir = "out_dir";

const char* kArgHelp = "h";
const char* kArgFullHelp = "help";

void
init_options(HbnProgramOptions* opts)
{
    opts->query_batch_size = datasize_to_u64(kDefault_QueryBatchSize);
    opts->query_upto = U64_MAX;

    opts->kmer_size = kDefault_KmerSize;
    opts->kmer_window = kDefault_KmerWindow;
    opts->rep_frac = kDefault_RepFrac;
    opts->max_kmer_occ = kDefault_MaxKmerOcc;
    opts->ddf = kDefault_DDF;
    opts->kmer_dist =kDefaultKmerDist;
    opts->chain_score = kDefaultChainScore;

    opts->perc_identity = kDefault_PercIdentity;
    opts->max_hsps_per_subject = kDefault_MaxAlignsPerSubject;
    opts->hitlist_size = kDefault_MaxTargetSubjects;

    opts->num_threads = kDefault_NumThreads;
    opts->outfmt = kDefault_OutFmt;
    opts->enzyme = NULL;
    opts->output = kDefault_Output;
}

void dump_usage_simple(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "\n");
    fprintf(out, "USAGE\n");
    fprintf(out, "  %s [OPTIONS] <enzyme_seq> <reference> <read1> [... <readN>]\n", pn);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION\n");
    fprintf(out, "  Alignment toolkit for long noisy chromosome conformation capture (3C) reads\n");

    fprintf(out, "\n");
    fprintf(out, "VERSION\n");
    const char* version = HBN_PACKAGE_VERSION;
    fprintf(out, "  %s\n", version);

    fprintf(out, "\n");
    fprintf(out, "*** Examples of familiar <enzyme_seq>:\n");
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

    fprintf(out, "\n");
    fprintf(out, "Type option '-help' to see details of optional arguments\n");    
}

void
dump_usage_full(const char* pn)
{
    FILE* out = stderr;
    fprintf(out, "\n");
    fprintf(out, "USAGE\n");
    fprintf(out, "  %s [OPTIONS] <enzyme_seq> <reference> <read1> [... <readN>]\n", pn);

    fprintf(out, "\n");
    fprintf(out, "DESCRIPTION\n");
    fprintf(out, "  PORE-C sequence alignment toolkit\n");

    fprintf(out, "\n");
    fprintf(out, "VERSION\n");
    const char* version = HBN_PACKAGE_VERSION;
    fprintf(out, "  %s\n", version);

    fprintf(out, "\n");
    fprintf(out, "*** Examples of familiar <enzyme_seq>:\n");
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

    fprintf(out, "\n");
    fprintf(out, "OPTIONAL ARGUMENTS\n");

    fprintf(out, "  -%s\n", kArgHelp);
    fprintf(out, "    Print USAGE and DESCRIPTION; ignore all other parameters\n");
    fprintf(out, "  -%s\n", kArgFullHelp);
    fprintf(out, "    Print USAGE, DESCRIPTION and ARGUMENTS; ignore all other parameters\n");

    fprintf(out, "\n");
    fprintf(out, "  *** Input sequence options\n");
    fprintf(out, "  -%s <DataSize>\n", kOpt_QueryBatchSize);
    fprintf(out, "    Batch size of queries (in bp) processed at one time\n");
    fprintf(out, "    Default = '%s'\n", kDefault_QueryBatchSize);
    fprintf(out, "  -%s <DataSize>\n", kOpt_MapUpto);
    fprintf(out, "    Map at most this size of queries and skip the rest queries\n");

    fprintf(out, "\n");
    fprintf(out, "  *** Candidate detection options\n");
    fprintf(out, "  -%s <Integer, %d..%d>\n", kOpt_KmerSize, kMinKmerSize, kMaxKmerSize);
    fprintf(out, "    Kmer size (length of best perfect match)\n");
    fprintf(out, "    Default = '%d'\n", kDefault_KmerSize);
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_KmerWindow, kMinKmerWindow);
    fprintf(out, "    Kmer sampling window size in reference sequences\n");
    fprintf(out, "    Default = '%d'\n", kDefault_KmerWindow);
    fprintf(out, "  -%s <Real, %g..%g>\n", kOpt_RepFrac, kMinRepFrac, kMaxRepFrac);
    fprintf(out, "    Filter out this fraction of most frequently occurring kmers\n");
    fprintf(out, "    Default = '%g'\n", kDefault_RepFrac);
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_KmerDist, kMinKmerDist);
    fprintf(out, "    Read distance and reference distance between adjacent\n");
    fprintf(out, "    matched kmers in a chain must not greater than this value\n");
    fprintf(out, "    Default = '%d'\n", kDefaultKmerDist);
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_MaxKmerOcc, kMinMaxKmerOcc);
    fprintf(out, "    Filter out kmers occuring larger than this value\n");
    fprintf(out, "    Default = '%d'\n", kDefault_MaxKmerOcc);
    fprintf(out, "  -%s <Real, >=%g>\n", kOpt_DDF, kMinDDF);
    fprintf(out, "    DDF factor\n");
    fprintf(out, "    Default = '%g'\n", kDefault_DDF);
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_ChainScore, kMinChainScore);
    fprintf(out, "    Minimum score of alignment candidates\n");
    fprintf(out, "    Default = '%d'\n", kDefaultChainScore);

    fprintf(out, "\n");
    fprintf(out, "  *** Restrict search or results\n");
    fprintf(out, "  -%s <Real, %g..%g>\n", kOpt_Perc_Identity, kMinPercIdentity, kMaxPercIdentity);
    fprintf(out, "    Minimum percentage of identity of alignments\n");
    fprintf(out, "    Default = '%g'\n", kDefault_PercIdentity);
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_MaxAlignsPerSubject, kMinMaxAlignsPercSubject);
    fprintf(out, "    Maximum number of alignments per reference sequence\n");
    fprintf(out, "    Deafult = '%d'\n", kDefault_MaxAlignsPerSubject);
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_MaxTargetSubjects, kMinMaxTargetSubjects);
    fprintf(out, "    Maximum number of aligned reference sequences to keep\n");
    fprintf(out, "    Default = '%d'\n", kDefault_MaxTargetSubjects);

    fprintf(out, "\n");
    fprintf(out, "  *** Miscellaneous options\n");
    fprintf(out, "  -%s <Integer, >=%d>\n", kOpt_NumThreads, kMinNumThreads);
    fprintf(out, "    Number of threads (CPUs) to use in the search\n");
    fprintf(out, "    Default = '%d'\n", kDefault_NumThreads);
    fprintf(out, "  -%s <String, Permissible values:", kOpt_OutFmt);
    for (int i = 0; i < eOutputFmt_Invalid; ++i) fprintf(out, " '%s'", kOutFmtNameList[i]);
    fprintf(out, ">\n");
    fprintf(out, "    Output format\n");
    fprintf(out, "    Default = '%s'\n", kOutFmtNameList[kDefault_OutFmt]);
    fprintf(out, "  -%s\n", kOpt_DumpByFile);
    fprintf(out, "    Output alignment results of each read file to its own results file\n");
    fprintf(out, "    Read file 'read.fq.gz' will produce result file 'read.fq.gz.sam' or 'read.fq.gz.paf'\n");
    fprintf(out, "    Default = %s\n", kDefault_DumpByFile ? "enabled" : "disabled");
    fprintf(out, "  -%s <Directory>\n", kOpt_OutputDir);
    fprintf(out, "    Alignment result files will be output into this directory\n");
    fprintf(out, "    Works only if option '-%s' is enabled\n", kOpt_DumpByFile);
    fprintf(out, "  -%s <file>\n", kOpt_Output);
    fprintf(out, "    Output file name.\n");
    fprintf(out, "    Default = '%s'\n", kDefault_Output);
}

void
s_PreCheckCmdLineArgs(int argc, char* argv[])
{
    for (int i = 1; i < argc; ++i) {
        if (argv[i][0] != '-') continue;
        if (strcmp(argv[i] + 1, kArgHelp) == 0) {
            dump_usage_simple(argv[0]);
            exit(0);
        }
        if (strcmp(argv[i] + 1, kArgFullHelp) == 0) {
            dump_usage_full(argv[0]);
            exit(0);
        }
    }
}

BOOL 
s_argument_is_supplied(int argc, char* argv[], int arg_i, int n_arg)
{
    if (arg_i + n_arg >= argc) {
        HBN_LOG("Mandatory argument to option '%s' is missing\n", argv[arg_i]);
        return FALSE;
    }
    return TRUE;
}

int parse_arguments(int argc, char* argv[], HbnProgramOptions* opts)
{
    s_PreCheckCmdLineArgs(argc, argv);
    init_options(opts);
    int i = 1;
    while (i < argc) {
        if (argv[i][0] != '-') break;

        if (strcmp(argv[i] + 1, kOpt_QueryBatchSize) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->query_batch_size = NStr::StringToUInt8_DataSize(NStr::CTempString(argv[i+1]));
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_MapUpto) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->query_upto = NStr::StringToUInt8_DataSize(NStr::CTempString(argv[i+1]));
            HBN_LOG("Map at most %s query bases", argv[i+1]);
            i += 2;
            continue;            
        }

        if (strcmp(argv[i] + 1, kOpt_KmerSize) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->kmer_size = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->kmer_size < kMinKmerSize || opts->kmer_size > kMaxKmerSize) {
                HBN_LOG("Illegal value '%d' to option '%s': must be in range [%d, %d]", opts->kmer_size, argv[i], kMinKmerSize, kMaxKmerSize);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_KmerWindow) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->kmer_window = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->kmer_window < kMinKmerWindow) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->kmer_window, argv[i], kMinKmerWindow);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_RepFrac) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->rep_frac = NStr::StringToDouble(NStr::CTempString(argv[i+1]));
            if (opts->rep_frac < kMinRepFrac || opts->rep_frac > kMaxRepFrac) {
                HBN_LOG("Illegal value '%g' to option '%s': must be in range [%g, %d]", opts->rep_frac, argv[i], kMinRepFrac, kMaxRepFrac);
                abort();
            }
            i += 2;
            continue;
        }    

        if (strcmp(argv[i] + 1, kOpt_KmerDist) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->kmer_dist = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->kmer_dist < kMinKmerDist) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->kmer_dist, argv[i], kMinKmerDist);
                abort();
            }
            i += 2;
            continue;
        }       

        if (strcmp(argv[i] + 1, kOpt_MaxKmerOcc) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->max_kmer_occ = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->max_kmer_occ < kMinMaxKmerOcc) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->max_kmer_occ, argv[i], kMinMaxKmerOcc);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_DDF) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->ddf = NStr::StringToDouble(NStr::CTempString(argv[i+1]));
            if (opts->ddf < kMinDDF) {
                HBN_LOG("Illegal value '%g' to option '%s': must be >=%g", opts->ddf, argv[i], kMinDDF);
                abort();
            }
            i += 2;
            continue;
        }          

        if (strcmp(argv[i] + 1, kOpt_ChainScore) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->chain_score = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->chain_score < kMinChainScore) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->chain_score, argv[i], kMinChainScore);
                abort();
            }
            i += 2;
            continue;
        } 

        if (strcmp(argv[i] + 1, kOpt_Perc_Identity) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->perc_identity = NStr::StringToDouble(NStr::CTempString(argv[i+1]));
            if (opts->perc_identity < kMinPercIdentity || opts->perc_identity > kMaxPercIdentity) {
                HBN_LOG("Illegal value '%g' to option '%s': must be in range [%g, %g]", opts->perc_identity, argv[i], kMinPercIdentity, kMaxPercIdentity);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_MaxAlignsPerSubject) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->max_hsps_per_subject = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->max_hsps_per_subject < kMinMaxAlignsPercSubject) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->max_hsps_per_subject, argv[i], kMinMaxAlignsPercSubject);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_MaxTargetSubjects) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->hitlist_size = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->hitlist_size < kMinMaxTargetSubjects) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->hitlist_size, argv[i], kMinMaxTargetSubjects);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_NumThreads) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->num_threads = NStr::StringToInt(NStr::CTempString(argv[i+1]));
            if (opts->num_threads < kMinNumThreads) {
                HBN_LOG("Illegal value '%d' to option '%s': must be >=%d", opts->num_threads, argv[i], kMinNumThreads);
                abort();
            }
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_OutFmt) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->outfmt = name_to_output_format(argv[i+1]);
            if (opts->outfmt == eOutputFmt_Invalid) HBN_ERR("Illegal output format name '%s'", argv[i+1]);
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_DumpByFile) == 0) {
            opts->dump_by_file = 1;
            i += 1;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_OutputDir) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->output_dir = argv[i+1];
            i += 2;
            continue;
        }

        if (strcmp(argv[i] + 1, kOpt_Output) == 0) {
            if (!s_argument_is_supplied(argc, argv, i, 1)) return 0;
            opts->output = argv[i+1];
            i += 2;
            continue;
        }

        HBN_LOG("Unrecognised option '%s'", argv[i]);
        return 0;
    }
    if (i + 3 > argc) return 0;
    opts->enzyme = argv[i];
    return i + 1;
}
