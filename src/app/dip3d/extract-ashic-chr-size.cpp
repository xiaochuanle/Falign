#include "../../corelib/fasta.hpp"

using namespace std;

int extract_ashic_chr_size_main(int argc, char* argv[])
{
    if (argc != 4) {
        fprintf(stderr, "USAGE:");
        fprintf(stderr, "%s %s reference chr-size-path\n", argv[0], argv[1]);
        return 1;
    }
    const char* reference = argv[2];
    const char* output = argv[3];

    hbn_dfopen(out, output, "w");
    HbnFastaReader in(reference);
    while (in.ReadOneSeq() != -1) {
        fprintf(out, "%s\t%zu\n", in.name().c_str(), in.sequence().size());
    }
    hbn_fclose(out);

    return 0;
}
