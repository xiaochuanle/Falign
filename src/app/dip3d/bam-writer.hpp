#ifndef __BAM_WRITER_HPP
#define __BAM_WRITER_HPP

#include "../../corelib/arg_parse.hpp"
#include "../../corelib/hbn_aux.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"
#include "../../htslib/sam.h"

#include <algorithm>
#include <limits>
#include <vector>

struct BamInfo
{
    int sid;
    int soff;
    bam1_t* bam;

    bool operator < (const BamInfo& rhs) const {
        return (sid < rhs.sid) || (sid == rhs.sid && soff < rhs.soff);
    }

    BamInfo(int _sid, int _soff, bam1_t* _bam):
        sid(_sid),
        soff(_soff),
        bam(_bam) {}
};

class BamChunkReader
{
public:
    BamChunkReader(const char* bam_path, const int num_threads) {
        M_in = sam_open(bam_path, "rb");
        hts_set_threads(M_in, num_threads);
        M_sam_hdr = sam_hdr_read(M_in);
        if (!M_sam_hdr) HBN_ERR("FAIL at reading BAM header from %s", bam_path);
        M_eof = false;
        M_bam_list = new bam1_t*[kBamChunkSize];

        load_next_bam_batch();
    }

    ~BamChunkReader() {
        sam_close(M_in);
        sam_hdr_destroy(M_sam_hdr);
        delete[] M_bam_list;
    }

    bam1_t* get_next_bam() {
        if (M_bam_list_idx >= M_bam_list_size) load_next_bam_batch();
        if (M_eof && M_bam_list_idx >= M_bam_list_size) return nullptr;
        return M_bam_list[M_bam_list_idx];
    }

    void advanve() {
        ++M_bam_list_idx;
    }

    sam_hdr_t* sam_hdr() {
        return M_sam_hdr;
    }

private:
    void load_next_bam_batch() {
        M_bam_list_size = 0;
        M_bam_list_idx = 0;

        size_t n = 0;
        while (n < kBamChunkSize) {
            bam1_t* bam = bam_init1();
            int r = sam_read1(M_in, M_sam_hdr, bam);
            if (r == -1) {
                bam_destroy1(bam);
                M_eof = true;
                break;
            }
            if (r < 0) HBN_ERR("FAIL at reading BAM record");
            M_bam_list[n] = bam;
            ++n;
        }
        M_bam_list_size = n;
        M_bam_list_idx = 0;
    }

private:
    static constexpr const size_t kBamChunkSize = 1024;
    samFile*    M_in;
    sam_hdr_t*  M_sam_hdr;
    bool        M_eof;
    bam1_t**    M_bam_list;
    size_t      M_bam_list_size;
    size_t      M_bam_list_idx;
};

static bam1_t*
select_next_bam(BamChunkReader** chunk_list, const int num_chunks)
{
    int best_sid = std::numeric_limits<int>::max();
    int best_soff = std::numeric_limits<int>::max();
    int best_i = -1;
    for (int i = 0; i < num_chunks; ++i) {
        bam1_t* bam = chunk_list[i]->get_next_bam();
        if (!bam) continue;
        int sid = bam->core.tid;
        int soff = bam->core.pos;
        bool r = (sid < best_sid) || (sid == best_sid && soff < best_soff);
        if (r) {
            best_sid = sid;
            best_soff = soff;
            best_i = i;
        }
    }
    if (best_i == -1) return nullptr;
    bam1_t* best_bam = chunk_list[best_i]->get_next_bam();
    chunk_list[best_i]->advanve();
    return best_bam;
}

class BamWriter
{
public:
    BamWriter(const char* output_path, sam_hdr_t* sam_hdr, const size_t chunk_size, const int num_threads) {
        M_out = sam_open(output_path, "wb");
        hts_set_threads(M_out, num_threads);
        M_sam_hdr = sam_hdr_dup(sam_hdr);
        int r = sam_hdr_write(M_out, M_sam_hdr);
        if (r) HBN_ERR("FAIL at writing BAM header to '%s'", output_path);
        M_chunk_size = chunk_size;
        M_stored_size = 0;
    }

    ~BamWriter() {
        dump_bams();
        sam_close(M_out);
        sam_hdr_destroy(M_sam_hdr);
    }

    void save_one_bam(bam1_t* bam) {
        M_bam_list.push_back(bam);
        M_stored_size += bam->core.l_qseq;
        if (M_stored_size >= M_chunk_size) {
            dump_bams();
        }
    }

private:
    void dump_bams() {
        std::string size = NStr::UInt8ToString_DataSize(M_stored_size);
        //HBN_LOG("Save %zu BAM records (%s)", M_bam_list.size(), size.c_str());
        for (auto bam : M_bam_list) {
            if (sam_write1(M_out, M_sam_hdr, bam) == -1) HBN_ERR("FAIL at writing BAM reocrd");
        }
        for (auto bam : M_bam_list) {
            bam_destroy1(bam);
        }
        M_bam_list.clear();
        M_stored_size = 0;
        //HBN_LOG("Done");
    }

private:
    samFile*                M_out;
    sam_hdr_t*              M_sam_hdr;
    size_t                  M_chunk_size;
    size_t                  M_stored_size;
    std::vector<bam1_t*>    M_bam_list;
};

static inline void
make_bam_chunk_path(const char* sorted_bam_path, const int chunk_id, char path[])
{
    snprintf(path, HBN_MAX_PATH_LEN, "%s.tmp.%d", sorted_bam_path, chunk_id);
}

#endif // __BAM_WRITER_HPP
