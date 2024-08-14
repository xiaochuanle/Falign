#include "../../corelib/pdqsort.h"
#include "frag-hap.hpp"

#include <vector>

using namespace std;

struct ContactHapInfo
{
    int read_id;
    int frag_id1;
    int frag_id2;
    int hp1;
    int hp2;
};

struct FragEval
{
    size_t total_frag;
    size_t tagged_frag;
    size_t acc_frag;

    size_t total_contact;
    size_t tagged_contact;
    size_t acc_contact;

    FragEval() {
        total_frag = 0;
        tagged_frag = 0;
        acc_frag = 0;

        total_contact = 0;
        tagged_contact = 0;
        acc_contact = 0;
    }

    void eval() {
        HBN_LOG("Frag eval:");
        double acc = 1.0 * acc_frag / tagged_frag;
        HBN_LOG("Accuracy: %g", acc);
        double recall = 1.0 * acc_frag / total_frag;
        HBN_LOG("Recall: %g", recall);

        HBN_LOG("Contact eval:");
        acc = 1.0 * acc_contact / tagged_contact;
        HBN_LOG("Accuracy: %g", acc);
        recall = 1.0 * acc_contact / total_contact;
        HBN_LOG("Recall: %g", recall);
    }
};

static void
s_make_contacts(FragHapInfo* fhia, size_t fhic, vector<ContactHapInfo>& contacts)
{
    ContactHapInfo contact;
    size_t i = 0;
    while (i < fhic) {
        size_t j = i + 1;
        while (j < fhic && fhia[i].read_id == fhia[j].read_id) ++j;

        for (size_t x = i; x < j; ++x) {
            for (size_t y = x + 1; y < j; ++y) {
                contact.read_id = fhia[i].read_id;
                contact.frag_id1 = fhia[x].frag_id;
                contact.frag_id2 = fhia[y].frag_id;
                contact.hp1 = fhia[x].hp;
                contact.hp2 = fhia[y].hp;
                contacts.push_back(contact);
            }
        }

        i = j;
    }
}

static void
s_eval_one_list(FragHapInfo* gt, FragHapInfo* cmp, size_t N, FragEval* eval)
{
    for (size_t i = 0; i < N; ++i) {
        hbn_assert(gt[i].read_id == cmp[i].read_id);
        hbn_assert(gt[i].frag_id == cmp[i].frag_id);
        if (gt[i].hp == 0) continue;
        ++eval->total_frag;
        if (cmp[i].hp == 0) continue;
        ++eval->tagged_frag;
        if (gt[i].hp == cmp[i].hp) ++eval->acc_frag;
    }

    vector<ContactHapInfo> gt_contacts;
    s_make_contacts(gt, N, gt_contacts);
    vector<ContactHapInfo> cmp_contacts;
    s_make_contacts(cmp, N, cmp_contacts);

    size_t M = gt_contacts.size();
    for (size_t i = 0; i < M; ++i) {
        ContactHapInfo& c1 = gt_contacts[i];
        ContactHapInfo& c2 = cmp_contacts[i];
        hbn_assert(c1.read_id == c2.read_id);
        hbn_assert(c1.frag_id1 == c2.frag_id1);
        hbn_assert(c1.frag_id2 == c2.frag_id2);
        if (c1.hp1 == 0 || c1.hp2 == 0) continue;
        ++eval->total_contact;

        if (c2.hp1 == 0 || c2.hp2 == 0) continue;
        ++eval->tagged_contact;

        if (c1.hp1 == c2.hp1 && c1.hp2 == c2.hp2) ++eval->acc_contact;
    }
}

int eval_hap_tag_main(int argc, char* argv[])
{
    if (argc < 4) {
        fprintf(stderr, "USAGE:\n");
        fprintf(stderr, "%s %s gt1 cmp1 [... gt-n cmp-n]\n", argv[0], argv[1]);
        exit (1);
    }
    FragEval eval;
    int i = 2;
    vector<FragHapInfo> gt_list, cmp_list;
    while (i < argc) {
        const char* gt_path = argv[i];
        const char* cmp_path = argv[i+1];
        HBN_LOG("Compare %s and %s", gt_path, cmp_path);
        gt_list.clear();
        cmp_list.clear();
        load_frag_hap_info_list(gt_path, gt_list);
        load_frag_hap_info_list(cmp_path, cmp_list);

        pdqsort(gt_list.begin(), gt_list.end(), [](const FragHapInfo& x, const FragHapInfo& y) { return (x.read_id < y.read_id) || (x.read_id == y.read_id && x.frag_id < y.frag_id); });
        pdqsort(cmp_list.begin(), cmp_list.end(), [](const FragHapInfo& x, const FragHapInfo& y) { return (x.read_id < y.read_id) || (x.read_id == y.read_id && x.frag_id < y.frag_id); });
        s_eval_one_list(gt_list.data(), cmp_list.data(), gt_list.size(), &eval);

        i += 2;
    }

    eval.eval();

    return 0;
}