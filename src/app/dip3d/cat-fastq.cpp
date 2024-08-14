#include "../../corelib/fasta.hpp"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

using namespace std;

int cat_fastq_main(int argc, char* argv[])
{
    int num_files = 0;
    size_t num_reads = 0;
    size_t num_bases = 0;

    for (int i = 2; i < argc; ++i) {
        fprintf(stderr, "%s\n", argv[i]);
        ++num_files;
        HbnFastaReader in(argv[i]);
        while (in.ReadOneSeq() != -1) {
            ++num_reads;
            num_bases += in.sequence().size();

            fprintf(stdout, "@%s\n", in.name().c_str());
            fprintf(stdout, "%s\n", in.sequence().c_str());
            fprintf(stdout, "+\n");
            fprintf(stdout, "%s\n", in.qual().c_str());
        }
    }

    string size = NStr::UInt8ToString_DataSize(num_bases);
    fprintf(stderr, "Files: %d, reads: %zu, bases: %s\n", num_files, num_reads, size.c_str());

    return 0;
}