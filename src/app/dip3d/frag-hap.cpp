#include "frag-hap.hpp"

#include "../../corelib/frag_id.hpp"
#include "../../corelib/line_reader.hpp"
#include "../../corelib/pdqsort.h"
#include "bam-writer.hpp"

using namespace std;

void save_frag_hap_info_list(FragHapInfo* a, size_t c, const char* output)
{
    hbn_dfopen(out, output, "w");
    fprintf(out, "## Format:\n");
    fprintf(out, "## read-id\tfrag-id\tchr-id\tchr-pos\tfrag-length\tmapQ\tidentity\thaplotype\n");
    for (size_t i = 0; i < c; ++i) dump_one_frag_hap_info(a[i], out, '\n');
    hbn_fclose(out);
}

void load_frag_hap_info_list(const char* path, std::vector<FragHapInfo>& fhi_list)
{
    HbnLineReader in(path);
    string sline;
    FragHapInfo fhi;
    while (in.ReadOneLine()) {
        NStr::CTempString line = *in;
        if (line.empty()) continue;
        if (line[0] == '#') continue;
        sline.assign(line.data(), line.size());
        parse_one_frag_hap_info(sline.c_str(), fhi);
        fhi_list.push_back(fhi);
    }
}

void load_frag_hap_info_list_from_bam(const char* path, std::vector<FragHapInfo>& fhi_list)
{
    BamChunkReader in(path, 8);
    FragHapInfo fhi;
    char tagname[2];
    while (1) {
        bam1_t* bam = in.get_next_bam();
        if (!bam) break;
        in.advanve();
        const char* qname = bam_get_qname(bam);
        extract_frag_id_from_name(qname, strlen(qname), &fhi.read_id, &fhi.frag_id, &fhi.subject_id, &fhi.subject_offset);
        fhi.frag_size = bam->core.l_qseq;
        fhi.mapQ = bam->core.qual;

        tagname[0] = 'p';
        tagname[1] = 'i';
        uint8_t* pitag = bam_aux_get(bam, tagname);
        if (!pitag) HBN_ERR("pi tag is not present in %s", path);
        fhi.identity = bam_aux2f(pitag);

        tagname[0] = 'H';
        tagname[1] = 'P';
        uint8_t* hptag = bam_aux_get(bam, tagname);
        fhi.hp = 0;
        if (hptag) fhi.hp = bam_aux2i(hptag);
        
        fhi_list.push_back(fhi);

        bam_destroy1(bam);
    }
}

void frag_and_contact_stats_hap_list(FragHapInfo* fhia, size_t fhic)
{
    pdqsort(fhia, fhia + fhic, [](const FragHapInfo& x, const FragHapInfo& y) { 
        return (x.read_id < y.read_id) || (x.read_id == y.read_id && x.frag_id < y.frag_id); });
    
    FragAndContactStats stats;
    size_t i = 0;
    while (i < fhic) {
        size_t j = i + 1;
        while (j < fhic && fhia[i].read_id == fhia[j].read_id) ++j;
        stats.add_one_list(fhia + i, j - i);
        i = j;
    }

    stats.dump_stats();
}