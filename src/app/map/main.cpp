#include "hbn_options.hpp"
#include "../../corelib/lookup_table.hpp"
#include "../../corelib/cstr_util.h"
#include "../../corelib/restrict_enzyme_loci_list.hpp"
#include "build_repeat_reference_regions.hpp"
#include "infer_enzyme.hpp"
#include "map_one_volume.hpp"
#include "necat_info.hpp"
#include "hbn_word_finder.hpp"

#include <sys/utsname.h>

using namespace std;

static size_t mapped_queries = 0;
static size_t mapped_bases = 0;

static void
s_validate_restrict_enzyme(const char* enzyme)
{
    RestrictEnzyme re;
    RestrictEnzyme_Init(enzyme, &re);
}

static void
s_map_one_db_volume(HbnUnpackedDatabase* subjects, 
    const HbnProgramOptions* opts, 
    int argc,
    char* argv[],
    vector<string>& query_list)
{
    HbnOutputs* out = new HbnOutputs(argc, argv, subjects, opts->output, opts->outfmt);
    RestrictEnzymeLociList* reloci_list = RestrictEnzymeLociListNew(subjects, opts->enzyme);
    HbnWordFinder* word_finder = new HbnWordFinder(*subjects,
        opts->repeat_bed,
        opts->skip_repeat_bed,
        opts->num_threads,
        opts->kmer_size,
        opts->kmer_window,
        opts->max_kmer_occ,
        opts->repeat_kmer_frac,
        opts->non_repeat_kmer_frac);
    HbnUnpackedDatabase* queries = new HbnUnpackedDatabase(query_list, opts->query_batch_size);

    while (queries->load_next_batch()) {
	    map_one_volume(subjects, queries, word_finder, opts, reloci_list, mapped_queries, out);
        mapped_queries += queries->NumSeqs();
        mapped_bases += queries->NumBases();
        char size_buf[64];
        u64_to_string_datasize(mapped_bases, size_buf);
        HBN_LOG("%12zu (%s) queries mapped", mapped_queries, size_buf);
        if (mapped_bases >= opts->query_upto) break;
    }

    delete word_finder;
    delete queries;
    RestrictEnzymeLociListFree(reloci_list);
    delete out;
}

static void
s_map_one_db_volume_with_enzyme_inference(HbnUnpackedDatabase* subjects, 
    const HbnProgramOptions* opts, 
    int argc,
    char* argv[],
    vector<string>& query_list)
{
    HbnOutputs* out = new HbnOutputs(argc, argv, subjects, opts->output, opts->outfmt);
    HbnWordFinder* word_finder = new HbnWordFinder(*subjects,
        opts->repeat_bed,
        opts->skip_repeat_bed,
        opts->num_threads,
        opts->kmer_size,
        opts->kmer_window,
        opts->max_kmer_occ,
        opts->repeat_kmer_frac,
        opts->non_repeat_kmer_frac);
    string last_enzyme;
    RestrictEnzymeLociList* reloci_list = nullptr;

    for (auto& query : query_list) {
        HbnUnpackedDatabase* queries = new HbnUnpackedDatabase(query.c_str(), opts->query_batch_size);
        queries->load_next_batch();
        HBN_LOG("Infer restriction enzyme for %s", query.c_str());
        string enzyme = infer_enzyme_mt(opts, subjects, queries, word_finder);
        if (enzyme.empty()) {
            HBN_LOG("FAIL at infere enzyme for %s", query.c_str());
            fprintf(stderr, "We skip this query file\n");
            delete queries;
            continue;
        } else {
            HBN_LOG("Successfully infer enzyme: %s", enzyme.c_str());
        }
        if (enzyme != last_enzyme) {
            last_enzyme = enzyme;
            delete reloci_list;
            reloci_list = RestrictEnzymeLociListNew(subjects, enzyme.c_str());
        }

        while (1) {
	        map_one_volume(subjects, queries, word_finder, opts, reloci_list, mapped_queries, out);
            mapped_queries += queries->NumSeqs();
            mapped_bases += queries->NumBases();
            char size_buf[64];
            u64_to_string_datasize(mapped_bases, size_buf);
            HBN_LOG("%12zu (%s) queries mapped", mapped_queries, size_buf);
            if (mapped_bases >= opts->query_upto) break;
            if (!queries->load_next_batch()) break;
        }

        delete queries;
    }

    delete reloci_list;
    delete word_finder;
    delete out;
}

int main(int argc, char* argv[])
{
    build_reference_repeat_regions(argc, argv);
    
    HbnProgramOptions _opts;
    HbnProgramOptions* opts = &_opts;
    vector<string> query_list;
    if (!opts->parse(argc, argv, query_list)) {
        opts->simple_usage(argc, argv);
        return EXIT_FAILURE;
    }
    HbnRunningInfo hbnrun;

    HbnUnpackedDatabase* subjects = new HbnUnpackedDatabase(opts->reference, U64_MAX);
    subjects->load_next_batch();
    if (opts->enzyme) {
        s_validate_restrict_enzyme(opts->enzyme);
        s_map_one_db_volume(subjects, opts, argc, argv, query_list);
    } else {
        s_map_one_db_volume_with_enzyme_inference(subjects, opts, argc, argv, query_list);
    }

    delete subjects;
    return 0;
}