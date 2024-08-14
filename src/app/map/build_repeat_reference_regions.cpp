#include "build_repeat_reference_regions.hpp"
#include "window_masker.hpp"

using namespace std;

const char* const kReferenceRepeatBuildCmd = "build-repeat";

void build_reference_repeat_regions(int argc, char* argv[])
{
    if (argc < 2) return;
    if (strcmp(argv[1], kReferenceRepeatBuildCmd)) return;
    if (argc != 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "  %s %s reference repeat-bed\n", argv[0], argv[1]);
        exit (EXIT_FAILURE);
    }
    const char* reference_path = argv[2];
    const char* repeat_bed_path = argv[3];
    
    HbnUnpackedDatabase reference(reference_path);
    reference.load_next_batch();
    vector<ChrIntv> repeat_intv_list;
    compute_repeat_regions(&reference, 1, repeat_intv_list);
    repeat_intv_stats(repeat_intv_list.data(), repeat_intv_list.size(), &reference);

    hbn_dfopen(out, repeat_bed_path, "w");
    for (auto& intv : repeat_intv_list) {
        const char* chr_name = reference.SeqName(intv.chr_id);
        fprintf(out, "%s\t%d\t%d", chr_name, intv.from, intv.to);
        if (intv.to - intv.from < 500) {
            fprintf(out, "\t");
            const u8* chr = reference.GetSequence(intv.chr_id);
            for (int i = intv.from; i < intv.to; ++i) {
                int c = chr[i];
                c = DECODE_RESIDUE(c);
                fprintf(out, "%c", c);
            }
        }
        fprintf(out, "\n");
    }
    hbn_fclose(out);

    exit(EXIT_SUCCESS);
}