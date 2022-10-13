#include "sam.h"

#include "format_aux.h"

void sam_dump_hdr_hd(FILE* out)
{
    fprintf(out, "@HD");
    fprintf(out, "\t");
    fprintf(out, "VN:%s", SAM_VERSION);
    fprintf(out, "\t");
    fprintf(out, "SO:unknown");
    fprintf(out, "\t");
    fprintf(out, "GO:query");
    fprintf(out, "\n");
}

void sam_dump_hdr_sq(const char* name, const int size, FILE* out)
{
    fprintf(out, "@SQ");
    fprintf(out, "\t");
    fprintf(out, "SN:%s", name);
    fprintf(out, "\t");
    fprintf(out, "LN:%d", size);
    fprintf(out, "\n");
}

void sam_dump_hdr_pg(const char* pg_id, const char* pg_name, const char* pg_version, int argc, char* argv[], FILE* out)
{
    fprintf(out, "@PG");
    fprintf(out, "\t");
    fprintf(out, "ID:%s", pg_id);
    fprintf(out, "\t");
    fprintf(out, "PN:%s", pg_name);
    fprintf(out, "\t");
    fprintf(out, "VN:%s", pg_version);
    fprintf(out, "\t");
    fprintf(out, "CL:");
    for (int i = 0; i < argc; ++i) {
        fprintf(out, "%s", argv[i]);
        if (i != argc - 1) fprintf(out, " ");
    }
    fprintf(out, "\n");
}

void sam_dump_hdr_rg(const char* rg_id, const char* sample_name, FILE* out)
{
    fprintf(out, "@RG");
    fprintf(out, "\t");
    fprintf(out, "ID:%s", rg_id);
    if (sample_name) {
        fprintf(out, "\t");
        fprintf(out, "SM:%s", sample_name);
    }
    fprintf(out, "\n");
}

////////////////////////////////////////////////////////////////

static void
dump_sam_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    if (qoff) ksprintf(out, "%dS", qoff);
    dump_cigar_string(qaln, saln, aln_size, out);
    if (qend < qsize) ksprintf(out, "%dS", qsize - qend);
}

void
dump_one_sam_result(
    const u8* qseq,
    const int qid,
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const int sid,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_quality,
    const int map_score,
    const double perc_identity,
    const double eff_perc_identity,
    const int is_complete_map,
    const BOOL dump_md,
    const char* rg_sample,
    kstring_t* out)
{
    const char tab = '\t';
    int flag = 0;
    if (qdir ==  REV) flag |= 0x10; // reverse query strand

    ksprintf(out, "%s", qname); /// 1) query name --- string  [!-?A-~]{1,254}
    kputc(tab, out);
    ksprintf(out, "%d", flag); /// 2) flag --- int  [0, 2^16-1]
    kputc(tab, out);
    ksprintf(out, "%s", sname); /// 3) subject name --- string  \*|[:rname:^*=][:rname:]*
    kputc(tab, out);
    ksprintf(out, "%d", soff + 1); /// 4) left most subject position (1-based) --- int  [0, 2^31-1]
    kputc(tab, out);
    if (map_quality >= 0) {  /// 5) mapq --- int [0, 2^8-1]
        hbn_assert(map_quality >= 0 && map_quality < 256);
        ksprintf(out, "%d", map_quality);
    } else {
        ksprintf(out, "%d", 60);
    }
    kputc(tab, out);
    dump_sam_cigar(qoff, qend, qsize, qaln, saln, aln_size, out); /// 6) cigar --- string  \*|([0-9]+[MIDNSHPX=])+
    kputc(tab, out);
    ksprintf(out, "*"); /// 7) reference name of the mate/next read --- string  \*|=|[:rname:^*=][:rname]*
    kputc(tab, out);
    ksprintf(out, "0"); /// 8) position of the mate/next read --- int [0, 2^31-1]
    kputc(tab, out);
    ksprintf(out, "%d", 0); /// 9) observed template LENgth --- int [-2^31 + 1, 2^31 - 1]
    kputc(tab, out);
    /// 10) segment SEQuence --- string \*|[A-Za-z=.]+
    for (int i = 0; i < qsize; ++i) {
        int c = qseq[i];
        c = DECODE_RESIDUE(c);
        kputc(c, out);
    }
    kputc(tab, out);
    ksprintf(out, "*"); /// 11) ASCII of Phred-scaled base QUALity+33 --- string [!-~]+

    int num_ident = 0;
    const char* q = qaln;
    const char* s = saln;
    for (int i = 0; i < aln_size; ++i) {
        if (q[i] == s[i]) ++num_ident;
    }
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", map_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", aln_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", map_score); ///  dp score
    if (dump_md) {
        kputc(tab, out);
        dump_md_string(qaln, saln, aln_size, out);
    }
    if (rg_sample) {
        kputc(tab, out);
        ksprintf(out, "RG:Z:%s", rg_sample);
    }
    /// qid
    kputc(tab, out);
    ksprintf(out, "qi:i:%d", qid);
    /// qdir
    kputc(tab, out);
    ksprintf(out, "qd:i:%d", qdir);
    /// qs
    kputc(tab, out);
    ksprintf(out, "qs:i:%d", qoff);
    /// qe
    kputc(tab, out);
    ksprintf(out, "qe:i:%d", qend);
    // ql
    kputc(tab, out);
    ksprintf(out, "ql:i:%d", qsize);
    /// sid
    kputc(tab, out);
    ksprintf(out, "si:i:%d", sid);
    /// ss
    kputc(tab, out);
    ksprintf(out, "ss:i:%d", soff);
    /// se
    kputc(tab, out);
    ksprintf(out, "se:i:%d", send);
    ///sl
    kputc(tab, out);
    ksprintf(out, "sl:i:%d", ssize);
    /// identity
    kputc(tab, out);
    ksprintf(out, "mc:f:%g", perc_identity);
    /// effective identity
    if (eff_perc_identity >= 0.0) {
        kputc(tab, out);
        ksprintf(out, "ec:f:%g", eff_perc_identity);
    }
    /// complete map
    if (is_complete_map >= 0) {
        kputc(tab, out);
        ksprintf(out, "fm:i:%d", is_complete_map);
    }
    kputc('\n', out);
}