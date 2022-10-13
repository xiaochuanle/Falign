#include "hbn_options.h"
#include "../../algo/align_format/sam.h"
#include "../../algo/seq_loader.h"
#include "../../algo/hbn_lookup_table.h"
#include "../../corelib/cstr_util.h"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/restrict_enzyme_loci_list.h"
#include "map_one_volume.hpp"
#include <sys/utsname.h>

static int mapped_queries = 0;
static size_t mapped_bases = 0;

static void
s_validate_restrict_enzyme(const char* enzyme)
{
    RestrictEnzyme re;
    RestrictEnzyme_Init(enzyme, &re);
}

static void
extract_pg_name(const char* argv_0, char pg_name[])
{
    int n = strlen(argv_0);
    int s = n;
    while (s) {
        --s;
        if (argv_0[s] == '/') {
            ++s;
            break;
        }
    }

    int i = 0;
    for (; s < n; ++i, ++s) pg_name[i] = argv_0[s];
    pg_name[i] = '\0';
}

static void 
dump_sam_hdr(
    int argc,
    char* argv[],
    SeqReader* subjects,
    FILE* out)
{
    sam_dump_hdr_hd(out);

    int n = SeqReader_NumSeqs(subjects);
    for (int i = 0; i < n; ++i) {
        const char* name = SeqReader_SeqName(subjects, i);
        int size = SeqReader_SeqSize(subjects, i);
        sam_dump_hdr_sq(name, size, out);
    }

    char pg_name[256];
    extract_pg_name(argv[0], pg_name);
    sam_dump_hdr_pg(pg_name, pg_name, HBN_PACKAGE_VERSION, argc, argv, out);
}

static void
s_map_one_db_volume(SeqReader* subjects, 
    const HbnProgramOptions* opts, 
    int query_file_cnt,
    char* query_file_list[],
    int argc,
    char* argv[])
{
    hbn_dfopen(out, opts->output, "w");
    if (opts->outfmt == eOutputFmt_SAM) dump_sam_hdr(argc, argv, subjects, out);
    RestrictEnzymeLociList* reloci_list = RestrictEnzymeLociListNew(subjects, opts->enzyme);
    HbnLookupTable* lktbl = HbnLookupTableNew(subjects, 
        reloci_list,
        opts->kmer_size,
        opts->kmer_window,
        opts->rep_frac,
        opts->max_kmer_occ,
        opts->num_threads);
    SeqReader* queries = SeqReaderNew(query_file_cnt, query_file_list);

    while (SeqReaderLoad(queries, opts->query_batch_size, TRUE)) {
        hbn_assert(kv_size(queries->fwd_unpacked_seq) == kv_size(queries->rev_unpacked_seq));
	    map_one_volume(subjects, queries, lktbl, opts, reloci_list, out);
        mapped_queries += SeqReader_NumSeqs(queries);
        mapped_bases += kv_size(queries->fwd_unpacked_seq);
        char size_buf[64];
        u64_to_string_datasize(mapped_bases, size_buf);
        HBN_LOG("%8d (%s) queries mapped", mapped_queries, size_buf);
        if (mapped_bases >= opts->query_upto) break;
    }

    lktbl = HbnLookupTableFree(lktbl);
    SeqReaderFree(queries);
    RestrictEnzymeLociListFree(reloci_list);
    hbn_fclose(out);
}

static void
s_make_results_file_name(const HbnProgramOptions* opts,
    const char* query_file_name,
    char path[])
{
    path[0] ='\0';
    if (opts->output_dir) sprintf(path, "%s/", opts->output_dir);

    int n = strlen(query_file_name);
    int i = n;
    while (i) {
        --i;
        if (query_file_name[i] == '/') {
            ++i;
            break;
        }
    }
    strcat(path, query_file_name + i);
    if (opts->outfmt == eOutputFmt_SAM) {
        strcat(path, ".sam");
    } else if (opts->outfmt == eOutputFmt_PAF) {
        strcat(path, ".paf");
    }
}

static void
s_map_one_db_volume_by_query_file(SeqReader* subjects, 
    const HbnProgramOptions* opts, 
    int query_file_cnt,
    char* query_file_list[],
    int argc,
    char* argv[])
{
    RestrictEnzymeLociList* reloci_list = RestrictEnzymeLociListNew(subjects, opts->enzyme);
    HbnLookupTable* lktbl = HbnLookupTableNew(subjects, 
        reloci_list,
        opts->kmer_size,
        opts->kmer_window,
        opts->rep_frac,
        opts->max_kmer_occ,
        opts->num_threads);

    SeqReader* queries = SeqReaderNew(query_file_cnt, query_file_list);

    while (1) {
        const char* name = SeqReader_CurrentFileName(queries);
        char path[HBN_MAX_PATH_LEN];
        s_make_results_file_name(opts, name, path);
        hbn_dfopen(out, path, "w");
        if (opts->outfmt == eOutputFmt_SAM) dump_sam_hdr(argc, argv, subjects, out);
        while (SeqReaderLoad(queries, opts->query_batch_size, TRUE)) {
            hbn_assert(kv_size(queries->fwd_unpacked_seq) == kv_size(queries->rev_unpacked_seq));
            map_one_volume(subjects, queries, lktbl, opts, reloci_list, out);
            mapped_queries += SeqReader_NumSeqs(queries);
            mapped_bases += kv_size(queries->fwd_unpacked_seq);
            char size_buf[64];
            u64_to_string_datasize(mapped_bases, size_buf);
            HBN_LOG("%8d (%s) queries mapped", mapped_queries, size_buf);
        }        
        hbn_fclose(out);
        if (!SeqReader_OpenNextFile(queries)) break;
    }
    lktbl = HbnLookupTableFree(lktbl);
    SeqReaderFree(queries);
    RestrictEnzymeLociListFree(reloci_list);
}

static void
s_dump_config()
{
    struct utsname _os_info_buf;
    struct utsname* os_info = NULL;
    if (uname(&_os_info_buf) == 0) os_info = &_os_info_buf;
    int n_threads = hbn_get_cpu_count();
    u64 mem_bytes = system_ram_bytes();
    char mem_size[64];
    u64_to_string_datasize(mem_bytes, mem_size);

    fprintf(stderr, "\n");
    fprintf(stderr, "Program:       falign\n");
    fprintf(stderr, "Version:       %s\n", HBN_PACKAGE_VERSION);
    fprintf(stderr, "Description:   Alignment toolkit for long noisy chromosome conformation capture (3C) reads\n");
    if (os_info) {
    fprintf(stderr, "System:        %s; %s; %s; %s\n", os_info->sysname, os_info->release, os_info->version, os_info->machine);
    }
    fprintf(stderr, "System CPUs:   %d\n", n_threads);
    fprintf(stderr, "System memory: %zu (%s)\n", mem_bytes, mem_size);
    fprintf(stderr, "\n");
}

int 
main(int argc, char* argv[])
{
    HbnProgramOptions* opts = (HbnProgramOptions*)calloc(1, sizeof(HbnProgramOptions));
    int ref_arg_idx = parse_arguments(argc, argv, opts);
    if (!ref_arg_idx) {
        dump_usage_simple(argv[0]);
        return EXIT_FAILURE;
    }
    s_validate_restrict_enzyme(opts->enzyme);
    s_dump_config();

    if (!opts->dump_by_file) {
        SeqReader* subjects = SeqReaderNew(1, argv + ref_arg_idx);
        int qry_arg_idx = ref_arg_idx + 1;
        while (SeqReaderLoad(subjects, U64_MAX, FALSE)) {
            s_map_one_db_volume(subjects, opts, argc - qry_arg_idx, argv + qry_arg_idx, argc, argv);
        }
        SeqReaderFree(subjects);
    } else {
        if (opts->output_dir) create_directory(opts->output_dir);
        SeqReader* subjects = SeqReaderNew(1, argv + ref_arg_idx);
        int qry_arg_idx = ref_arg_idx + 1;
        while (SeqReaderLoad(subjects, U64_MAX, FALSE)) {
            s_map_one_db_volume_by_query_file(subjects, opts, argc - qry_arg_idx, argv + qry_arg_idx, argc, argv);
        }
        SeqReaderFree(subjects);
    }

    free(opts);
    return 0;
}