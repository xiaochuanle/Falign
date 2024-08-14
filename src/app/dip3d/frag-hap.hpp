#ifndef __FRAG_HAP_HPP
#define __FRAG_HAP_HPP

#include "../../corelib/hbn_aux.h"

#include <cstdlib>

#include <vector>

struct FragHapInfo
{
    int read_id;
    int frag_id;
    int subject_id;
    int subject_offset;
    int frag_size;
    int mapQ;
    double identity;
    int hp;
    int amb_hp;
};

static inline void
dump_one_frag_hap_info(FragHapInfo& fhi, FILE* out, const char newline = '\n')
{
    fprintf(out, "%d\t%d\t%d\t%d\t%d\t%d\t%g\t%d%c", 
        fhi.read_id, fhi.frag_id, fhi.subject_id, fhi.subject_offset, fhi.frag_size, fhi.mapQ, fhi.identity, fhi.hp, newline);
}

static inline void
parse_one_frag_hap_info(const char* line, FragHapInfo& fhi)
{
    HBN_SCANF(sscanf, line, 8, "%d%d%d%d%d%d%lf%d",
        &fhi.read_id, &fhi.frag_id, &fhi.subject_id, &fhi.subject_offset, &fhi.frag_size, &fhi.mapQ, &fhi.identity, &fhi.hp);
}

void save_frag_hap_info_list(FragHapInfo* a, size_t c, const char* output);

void load_frag_hap_info_list(const char* path, std::vector<FragHapInfo>& fhi_list);

void load_frag_hap_info_list_from_bam(const char* path, std::vector<FragHapInfo>& fhi_list);

/////////////////

struct ContactStats
{
    ContactStats(int max_dist) {
        M_max_dist = max_dist;
        M_contacts = 0;
        M_phased_contacts = 0;
        M_h1_contacts = 0;
        M_h2_contacts = 0;
        M_h_trans_contacts = 0;
    }

    void add_one_list(FragHapInfo* a, int c) {
        for (int i = 0; i < c; ++i) {
            int hi = a[i].hp;
            for (int j = i + 1; j < c; ++j) {
                hbn_assert(a[i].read_id == a[j].read_id);
                hbn_assert(a[i].frag_id != a[j].frag_id);
                if (a[i].subject_id != a[j].subject_id) continue;
                if (abs(a[i].subject_offset - a[j].subject_offset) > M_max_dist) continue;
                int hj = a[j].hp;

                ++M_contacts;
                if (hi < 1 || hj < 1) continue;
                ++M_phased_contacts;
                if (hi == hj && hi == 1) {
                    ++M_h1_contacts;
                } else if (hi == hj && hi == 2) {
                    ++M_h2_contacts;
                } else {
                    hbn_assert(hi != hj);
                    ++M_h_trans_contacts;
                }
            }
        }
    }

    void dump_stats() {
        if (M_contacts == 0) return;
        double phased_frac = 100.0 * M_phased_contacts / M_contacts;
        double h1_frac = 0.0, h2_frac = 0.0, h_trans_frac = 0.0;
        if (M_phased_contacts > 0) {
            h1_frac = 100.0 * M_h1_contacts / M_phased_contacts;
            h2_frac = 100.0 * M_h2_contacts / M_phased_contacts;
            h_trans_frac = 100.0 * M_h_trans_contacts / M_phased_contacts;
        }

        HBN_LOG("Contact stats for contact dist <= %d", M_max_dist);
        fprintf(stderr, "   Total contacts: %zu\n", M_contacts);
        fprintf(stderr, "   Phased contacts: %zu (%g%%)\n", M_phased_contacts, phased_frac);
        fprintf(stderr, "   H1 contacts: %zu (%g%%)\n", M_h1_contacts, h1_frac);
        fprintf(stderr, "   H2 contacts: %zu (%g%%)\n", M_h2_contacts, h2_frac);
        fprintf(stderr, "   H-trans contacts: %zu (%g%%)\n", M_h_trans_contacts, h_trans_frac);
    }

    int     M_max_dist;
    size_t  M_contacts;
    size_t  M_phased_contacts;
    size_t  M_h1_contacts;
    size_t  M_h2_contacts;
    size_t  M_h_trans_contacts;
};

struct FragAndContactStats
{
    FragAndContactStats() {
        M_contact_stats.push_back(ContactStats(5000000));
        M_contact_stats.push_back(ContactStats(50000000));
        M_contact_stats.push_back(ContactStats(INT32_MAX));

        M_reads = 0;
        M_reads_with_phased_frags = 0;
        M_frags = 0;
        M_phased_frags = 0;
        M_h1_frags = 0;
        M_h2_frags = 0;
    }

    void add_one_list(FragHapInfo* a, int c) {
        M_frags += c;
        bool has_phased_frags = false;
        for (int i = 0; i < c; ++i) {
            if (a[i].hp < 1) continue;
            ++M_phased_frags;
            if (a[i].hp == 1) ++M_h1_frags;
            if (a[i].hp == 2) ++M_h2_frags;
            has_phased_frags = true;
        }
        for (auto& x : M_contact_stats) x.add_one_list(a, c);
        ++M_reads;
        M_reads_with_phased_frags += has_phased_frags;
    }

    void dump_stats() {
        fprintf(stderr, "   Total reads: %zu\n", M_reads);
        if (M_reads == 0) return;
        double frac = 100.0 * M_reads_with_phased_frags / M_reads;
        fprintf(stderr, "   Reads with phased frags: %zu (%g%%)\n", M_reads_with_phased_frags, frac);
        fprintf(stderr, "   Total frags: %zu\n", M_frags);
        if (M_frags == 0) return;

        double phased_frac = 100.0 * M_phased_frags / M_frags;
        double h1_frac = 100.0 * M_h1_frags / M_frags;
        double h2_frac = 100.0 * M_h2_frags / M_frags;
        fprintf(stderr, "   Phased frags: %zu (%g%%)\n", M_phased_frags, phased_frac);
        fprintf(stderr, "   H1 frags: %zu (%G%%)\n", M_h1_frags, h1_frac);
        fprintf(stderr, "   H2 frags: %zu (%g%%)\n", M_h2_frags, h2_frac);
        for (auto& c : M_contact_stats) c.dump_stats();
    }

    std::vector<ContactStats> M_contact_stats;
    size_t M_reads;
    size_t M_reads_with_phased_frags;
    size_t M_frags;
    size_t M_phased_frags;
    size_t M_h1_frags;
    size_t M_h2_frags;
};

void frag_and_contact_stats_hap_list(FragHapInfo* fhia, size_t fhic);

#endif // __FRAG_HAP_HPP
