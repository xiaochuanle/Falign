#include "x-tag-bam.hpp"

#include "../../corelib/frag_id.hpp"
#include "bam-writer.hpp"

#include <map>
#include <mutex>
#include <vector>

using namespace std;

class BamTagWorkThreadData
{
public:
    BamTagWorkThreadData(map<pair<int, int>, int>* frag_haps, vector<bam1_t*>* bam_list) {
        M_frag_haps = frag_haps;
        M_bam_list = bam_list;
        M_bam_idx = 0;
    }

    map<pair<int, int>, int>* frag_haps() {
        return M_frag_haps;
    }

    bam1_t* get_next_bam() {
        lock_guard<mutex> _(M_bam_mutex);
        if (M_bam_idx >= M_bam_list->size()) return nullptr;
        return (*M_bam_list)[M_bam_idx++];
    }

private:
    map<pair<int, int>, int>*   M_frag_haps;
    vector<bam1_t*>*            M_bam_list;
    size_t                      M_bam_idx;
    mutex                       M_bam_mutex;
};

static void*
tag_bam_thread(void* params)
{
    BamTagWorkThreadData* data = (BamTagWorkThreadData*)(params);
    map<pair<int, int>, int>* frag_haps = data->frag_haps();
    char hptagname[2] = { 'H', 'P' };
    bam1_t* bam = nullptr;
    int read_id = -1, frag_id = -1, subject_id = -1, subject_offset = -1;

    while ((bam = data->get_next_bam())) {
        uint8_t* hptag = bam_aux_get(bam, hptagname);
        if (hptag) bam_aux_del(bam, hptag);

        const char* qn = bam_get_qname(bam);
        const int qnl = strlen(qn);
        extract_frag_id_from_name(qn, qnl, &read_id, &frag_id, &subject_id, &subject_offset);
        auto pos = frag_haps->find(pair<int, int>(read_id, frag_id));
        if (pos == frag_haps->end()) continue;

        int hp = pos->second;
        if (bam_aux_update_int(bam, hptagname, hp)) HBN_ERR("FAIL at HP tag to BAM record");    
    }

    return nullptr;
}

void tag_bam_mt(FragHapInfo* fhia, size_t fhic, const int num_threads, 
    const char* input_bam_path, const char* tagged_bam_path)
{
    HBN_LOG("Haplotagging %s", input_bam_path);
    map<pair<int, int>, int> frag_haps;
    for (size_t i = 0; i < fhic; ++i) {
        if (fhia[i].hp < 1) continue;
        frag_haps[pair<int, int>(fhia[i].read_id, fhia[i].frag_id)] = fhia[i].hp;
    }
    BamChunkReader in(input_bam_path, 8);
    vector<bam1_t*> bam_list;
    BamWriter out(tagged_bam_path, in.sam_hdr(), 1000000000, 8);
    pthread_t jobs[num_threads];
    string size;
    while (1) {
        size_t num_frags = 0, num_bases = 0;
        bool eof = false;
        while (1) {
            bam1_t* bam = in.get_next_bam();
            if (!bam) {
                eof = true;
                break;
            }
            in.advanve();
            bam_list.push_back(bam);
            ++num_frags;
            num_bases += bam->core.l_qseq;
            if (num_bases >= 10000000000) break;
        }
        if (bam_list.empty()) break;
        size = NStr::UInt8ToString_DataSize(num_bases);
        HBN_LOG("Load %zu frags (%s)", num_frags, size.c_str());

        BamTagWorkThreadData data(&frag_haps, &bam_list);
        for (int i = 0; i < num_threads; ++i) {
            pthread_create(jobs + i, nullptr, tag_bam_thread, &data);
        }
        for (int i = 0; i < num_threads; ++i) {
            pthread_join(jobs[i], nullptr);
        }

        for (auto bam : bam_list) out.save_one_bam(bam);
        bam_list.clear();
        if (eof) break;
    }
}

void tag_bam_st(FragHapInfo* fhia, size_t fhic, const int num_threads, 
    const char* input_bam_path, const char* tagged_bam_path)
{
    HBN_LOG("Haplotagging %s", input_bam_path);
    map<pair<int, int>, int> frag_haps;
    for (size_t i = 0; i < fhic; ++i) {
        if (fhia[i].hp < 1) continue;
        frag_haps[pair<int, int>(fhia[i].read_id, fhia[i].frag_id)] = fhia[i].hp;
    }

    samFile* in = sam_open(input_bam_path, "rb");
    hts_set_threads(in, 8);
    if (!in) abort();
    sam_hdr_t* hdr = sam_hdr_read(in);
    if (!hdr) HBN_ERR("Could not read BAM header from %s", input_bam_path);
    samFile* out = sam_open(tagged_bam_path, "wb");
    if (!out) HBN_ERR("Could not open %s for writing", tagged_bam_path);
    hts_set_threads(out, 8);
    int r = sam_hdr_write(out, hdr);
    if (r) HBN_ERR("FAIL at writing BAM header to '%s'", tagged_bam_path);

    bam1_t* bam = bam_init1();
    char hptagname[2] = { 'H', 'P' };
    int read_id = -1, frag_id = -1, subject_id = -1, subject_offset = -1;
    while (1) {
        r = sam_read1(in, hdr, bam);
        if (r == -1) break;
        if (r < 0) HBN_ERR("FAIL at reading BAM record from %s", input_bam_path);

        const char* qn = bam_get_qname(bam);
        const int qnl = strlen(qn);
        extract_frag_id_from_name(qn, qnl, &read_id, &frag_id, &subject_id, &subject_offset);
        auto pos = frag_haps.find(pair<int, int>(read_id, frag_id));
        if (pos != frag_haps.end()) {
            int hp = pos->second;
            if (bam_aux_update_int(bam, hptagname, hp)) HBN_ERR("FAIL at add HP tag to BAM record");   
        }

        if (sam_write1(out, hdr, bam) == -1) HBN_ERR("FAIL at writing BAM reocrd to %s", tagged_bam_path);
    }

    bam_destroy1(bam);
    bam_hdr_destroy(hdr);
    sam_close(out);
    sam_close(in);
}