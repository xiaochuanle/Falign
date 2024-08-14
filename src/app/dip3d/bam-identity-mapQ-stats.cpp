#include "bam-writer.hpp"

using namespace std;

int bam_identity_mapq_stats_main(int argc, char* argv[])
{
    if (argc != 3) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s bam\n", argv[0], argv[1]);
        exit(1);
    }
    const char* bam_path = argv[2];

    samFile* in = sam_open(bam_path, "rb");
    if (!in) HBN_ERR("Could not open file %s for reading", bam_path);
    hts_set_threads(in, 8);
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", bam_path);

    char pitagname[2] = { 'p', 'i' };
    size_t sum_mapQ = 0;
    double sum_identity = 0;
    size_t num_frag = 0;
    size_t num_bases = 0;

    bam1_t* bam = nullptr;
    while (1) {
        if (!bam) bam = bam_init1();
        int r = sam_read1(in, hdr, bam);
        if (r == -1) break;
        if (r < 0) HBN_ERR("FAIL at reading BAM record from %s", bam_path);

        int frag_size = bam->core.l_qseq;
        int mapQ = bam->core.qual;
        uint8_t* pitag = bam_aux_get(bam, pitagname);
        if (!pitag) HBN_ERR("Could not retrieve alignment identity from %c%c BAM tag", pitagname[0], pitagname[1]);
        double pi = bam_aux2f(pitag);         
    
        ++num_frag;
        sum_mapQ += mapQ;
        sum_identity += pi;
        num_bases += frag_size;
    }
    if (bam) bam_destroy1(bam); 
    bam = nullptr;
    sam_close(in);
    sam_hdr_destroy(hdr);

    int avg_frag_size = num_bases / num_frag;
    int avg_mapQ = sum_mapQ / num_frag;
    double avg_identity = sum_identity / num_frag;
    fprintf(stderr, "Frags: %zu, avg-frag-size: %d, avg-mapQ: %d, avg-identity: %g\n", num_frag, avg_frag_size, avg_mapQ, avg_identity);

    return 0;
}