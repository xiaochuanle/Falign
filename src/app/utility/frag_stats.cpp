#include "paf_reader.hpp"
#include "../../corelib/fasta.h"

#include <map>
#include <string>
#include <vector>

using namespace std;

int frag_stats_main(int argc, char* argv[])
{
    const char* paf_path = argv[2];
    const char* fq_path = argv[3];

    map<string, int> frag_stats;
    vector<PAF> paf_list;
    vector<string> ref_name_list;
    string read_name;
    
    PAFReader paf(paf_path);
    while(paf.read_next_paf_list(paf_list, read_name, ref_name_list)) {
        frag_stats[read_name] = paf_list.size();
    }

    HbnFastaReader* fq = HbnFastaReaderNew(fq_path);
    while (!HbnLineReaderAtEof(fq->line_reader)) {
        HbnFastaReaderReadOneSeq(fq);
        int cnt = 0;
        const char* s = ks_s(fq->comment);
        const int sl = ks_size(fq->comment);
        for (int i = 0; i < sl; ++i) if (s[i] == ';') ++cnt;
        ++cnt;

        read_name.assign(ks_s(fq->name), ks_size(fq->name));
        auto pos = frag_stats.find(read_name);
        if (pos == frag_stats.end()) continue;

        if (pos->second - cnt > 0) fprintf(stderr, "%s\t%d\t%d\n", read_name.c_str(), pos->second, cnt);
    }
    HbnFastaReaderFree(fq);

    return 0;
}
