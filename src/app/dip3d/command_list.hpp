#ifndef __COMMAND_LIST_HPP
#define __COMMAND_LIST_HPP

#include <map>
#include <string>

typedef int main_fun_type(int argc, char* argv[]);

int bam_identity_mapq_stats_main(int argc, char* argv[]);
int cat_fastq_main(int argc, char* argv[]);
int compress_vcf_main(int argc, char* argv[]);
int eval_hap_tag_main(int argc, char* argv[]);
int extract_ashic_chr_size_main(int argc, char* argv[]);
int extract_chr_bam_main(int argc, char* argv[]);
int extract_chr_main(int argc, char* argv[]);
int extract_mvp_het_snp_main(int argc, char* argv[]);
int frag_to_ashic_read_pair_main(int argc, char* argv[]);
int haplo_tag_stats_main(int argc, char* argv[]);
int haplotype_consensus_main(int argc, char* argv[]);
int index_bam_main(int argc, char* argv[]);
int index_fasta_main(int argc, char* argv[]);
int index_vcf_main(int argc, char* argv[]);
int make_pore_c_frag_pair_main(int argc, char* argv[]);
int merge_snp_bam_main(int argc, char* argv[]);
int paf_frag_stats_main(int argc, char* argv[]);
int sample_vcf_main(int argc, char* argv[]);
int select_chr_snp_bam_main(int argc, char* argv[]);
int select_snp_bam_main(int argc, char* argv[]);
int sort_bam_main(int argc, char* argv[]);
int split_bam_main(int argc, char* argv[]);
int split_vcf_main(int argc, char* argv[]);
int tag_bam_main(int argc, char* argv[]);

class Dip3dCommandList
{
public:
    Dip3dCommandList() {
        add_one_cmd("bam-identity-mapQ-stats", bam_identity_mapq_stats_main);
        add_one_cmd("cat-fastq", cat_fastq_main);
        add_one_cmd("compress-vcf", compress_vcf_main);
        add_one_cmd("eval-frag-hap-tag", eval_hap_tag_main);
        add_one_cmd("extract-ashic-chr-size", extract_ashic_chr_size_main);
        add_one_cmd("extract-chr-bam", extract_chr_bam_main);
        add_one_cmd("extract-chr", extract_chr_main);
        add_one_cmd("extract-mvp-het-snp", extract_mvp_het_snp_main);
        add_one_cmd("frag-to-ashic-read-pair", frag_to_ashic_read_pair_main);
        add_one_cmd("haplo-tag-stats", haplo_tag_stats_main);
        add_one_cmd("haplotype-consensus", haplotype_consensus_main);
        add_one_cmd("index-bam", index_bam_main);
        add_one_cmd("index-fasta", index_fasta_main);
        add_one_cmd("index-vcf", index_vcf_main);
        add_one_cmd("make-pore-c-frag-pair", make_pore_c_frag_pair_main);
        add_one_cmd("merge-snp-bam", merge_snp_bam_main);
        add_one_cmd("paf-frag-stats", paf_frag_stats_main);
        add_one_cmd("sample-vcf", sample_vcf_main);
        add_one_cmd("select-chr-snp-bam", select_chr_snp_bam_main);
        add_one_cmd("select-snp-bam", select_snp_bam_main);
        add_one_cmd("sort-bam", sort_bam_main);
        add_one_cmd("split-bam", split_bam_main);
        add_one_cmd("split-vcf", split_vcf_main);
        add_one_cmd("tag-bam", tag_bam_main);
    }

    void add_one_cmd(const char* cmdname, main_fun_type* fun) {
        M_cmds[std::string(cmdname)] = fun;
    }

    void dump_cmds(FILE* out = stderr) {
        fprintf(stderr, "\n");
        fprintf(stderr, "COMMANDS:\n");
        for (auto& cmd : M_cmds) {
            fprintf(stderr, "  %s\n", cmd.first.c_str());
        }
    }

    int run_cmd(int argc, char* argv[]) {
        auto pos = M_cmds.find(std::string(argv[1]));
        if (pos == M_cmds.end()) {
            fprintf(stderr, "ERROR: Unrecognised command '%s'\n", argv[1]);
            return 1;
        }
        return (*pos->second)(argc, argv);
    }

private:
    std::map<std::string, main_fun_type*>    M_cmds;
};

#endif // __COMMAND_LIST_HPP