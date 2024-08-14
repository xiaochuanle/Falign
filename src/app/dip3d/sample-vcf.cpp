#include "../../corelib/line_reader.hpp"

#include <random>
#include <string>

using namespace std;

int sample_vcf_main(int argc, char* argv[])
{
    if (argc != 5) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s input-vcf sample-frac sampled-vcf\n", argv[0], argv[1]);
        exit(1);
    }

    const char* input_vcf_path = argv[2];
    const double frac = atof(argv[3]);
    const char* output_vcf_path = argv[4];

    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<> dist;

    HbnLineReader in(input_vcf_path);
    hbn_dfopen(out, output_vcf_path, "w");
    size_t total_vcf = 0, sampled_vcf = 0;
    while (in.ReadOneLine()) {
        NStr::CTempString line = *in;
        if (line.empty()) continue;
        if (line[0] == '#') {
            hbn_fwrite(line.data(), 1, line.size(), out);
            fprintf(out, "\n");
            continue;
        }
        ++total_vcf;
        if (dist(gen) > frac) continue;
        ++sampled_vcf;
        hbn_fwrite(line.data(), 1, line.size(), out);
        fprintf(out, "\n");
    }
    hbn_fclose(out);

    double p = 1.0 * sampled_vcf / total_vcf;
    fprintf(stderr, "Total VCF records: %zu\n", total_vcf);
    fprintf(stderr, "Sampled VCF records: %zu (%g)\n", sampled_vcf, p);

    return 0;
}