#include "../../corelib/fasta.hpp"

#include <string>
#include <unordered_set>

using namespace std;

int extract_chr_main(int argc, char* argv[])
{
    if (argc < 5) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s input-fasta output-fasta chr-1 [..chr-n]\n", argv[0], argv[1]);
        exit (1);
    }
    const char* input = argv[2];
    const char* output = argv[3];
    unordered_set<string> chrs;
    for (int i = 4; i < argc; ++i) chrs.insert(argv[i]);

    HbnFastaReader in(input);
    hbn_dfopen(out, output, "w");
    while (in.ReadOneSeq() != -1) {
        if (!chrs.count(in.name())) continue;
        fprintf(stderr, "%s\t%zu\n", in.name().c_str(), in.sequence().size());
        fprintf(out, ">%s\n", in.name().c_str());
        fprintf(out, "%s\n", in.sequence().c_str());
    }
    hbn_fclose(out);

    return 0;
}