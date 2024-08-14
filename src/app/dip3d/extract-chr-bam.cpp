#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../corelib/pdqsort.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "../../htslib/sam.h"
#include "bam-writer.hpp"

#include <algorithm>
#include <vector>

using namespace std;

int extract_chr_bam_main(int argc, char* argv[])
{
    if (argc != 5) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s input-bam chr chr-bam\n", argv[0], argv[1]);
        exit (1);
    }
    const char* input_bam_path = argv[2];
    const char* chr = argv[3];
    const char* chr_bam_path = argv[4];
    const int num_threads = 8;
    const size_t chunk_size = 1000000000;

    samFile* in = sam_open(input_bam_path, "rb");
    hts_set_threads(in, num_threads);
    sam_hdr_t* hdrin = sam_hdr_read(in);
    BamWriter* out = new BamWriter(chr_bam_path, hdrin, chunk_size, num_threads);
    const int chr_id = sam_hdr_name2tid(hdrin, chr);
    if (chr_id == -1) {
        HBN_ERR("Chromosome %s does not exist in %s", chr, input_bam_path);
    } else if (chr_id == -2) {
        HBN_ERR("FAIL at parsing BAM header from %s", input_bam_path);
    }
    size_t num_frags = 0, num_bases = 0;
    bam1_t* bam = nullptr;
    while (1) {
        if (!bam) bam = bam_init1();
        int r = sam_read1(in, hdrin, bam);
        if (r == -1) {
            bam_destroy1(bam);
            break;
        }
        if (r < 0) HBN_ERR("FAIL at reading BAM record");
	if (bam->core.tid != chr_id) {
		continue;
	}
        ++num_frags;
        num_bases += bam->core.l_qseq;
        out->save_one_bam(bam);
	bam = nullptr;
    }
    delete out;
    sam_hdr_destroy(hdrin);
    sam_close(in);

    string size = NStr::UInt8ToString_DataSize(num_bases);
    HBN_LOG("Extract %zu BAM records (%s) for %s", num_frags, size.c_str(), chr);

    return 0;
}
