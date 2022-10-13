#include "paf.h"

#include "format_aux.h"

void
print_paf_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    ksprintf(out, "cg:Z:");
    if (qoff) ksprintf(out, "%dS", qoff);
    dump_cigar_string(qaln, saln, aln_size, out);
    if (qend < qsize) ksprintf(out, "%dS", qsize - qend);
}

#if 0
                                         1 | string | Query sequence name                                     |
                                      |  2 |  int   | Query sequence length                                   |
                                      |  3 |  int   | Query start coordinate (0-based)                        |
                                      |  4 |  int   | Query end coordinate (0-based)                          |
                                      |  5 |  char  | `+' if query/target on the same strand; `-' if opposite |
                                      |  6 | string | Target sequence name                                    |
                                      |  7 |  int   | Target sequence length                                  |
                                      |  8 |  int   | Target start coordinate on the original strand          |
                                      |  9 |  int   | Target end coordinate on the original strand            |
                                      | 10 |  int   | Number of matching bases in the mapping                 |
                                      | 11 |  int   | Number bases, including gaps, in the mapping            |
                                      | 12 |  int   | Mapping quality (0-255 with 255 for missing)  
#endif 

void
dump_one_paf_result(
    const char* qname,
    const int qdir,
    const int qoff,
    const int qend,
    const int qsize,
    const char* sname,
    const int soff,
    const int send,
    const int ssize,
    const char* qaln,
    const char* saln,
    const int aln_size,
    const int map_score,
    const BOOL dump_cigar,
    const BOOL dump_md,
    kstring_t* out)
{
    const char tab = '\t';
    ksprintf(out, "%s", qname); /// 1) query name
    kputc(tab, out);
    ksprintf(out, "%d", qsize); /// 2) query length
    kputc(tab, out);
    ksprintf(out, "%d", qoff); /// 3) query start coordinate (0-based)
    kputc(tab, out);
    ksprintf(out, "%d", qend); /// 4) query end coordinate (0-based)
    kputc(tab, out);
    /// 5) '+' if query and subject on the same strand; '-' if opposite
    /// the subject sequence is always on forward strand
    if (qdir == FWD) {
        kputc('+', out);
    } else {
        kputc('-', out);
    }
    kputc(tab, out);
    ksprintf(out, "%s", sname); /// 6) subject name
    kputc(tab, out);
    ksprintf(out, "%d", ssize); /// 7) subject length
    kputc(tab,out);
    ksprintf(out, "%d", soff); /// 8) subject start coordinate (0-based)
    kputc(tab, out);
    ksprintf(out, "%d", send); /// 9) subject end coordinate (0-based)
    kputc(tab, out);
    
    int num_ident = 0;
    const char* q = qaln;
    const char* s = saln;
    for (int i = 0; i < aln_size; ++i) if (q[i] == s[i]) ++num_ident;
    ksprintf(out, "%d", num_ident); /// 10) number of matching bases in the alignment
    kputc(tab, out);
    ksprintf(out, "%d", aln_size); /// 11) number of bases in the alignment (including gaps and mismatch bases)
    kputc(tab, out);
    ksprintf(out, "%d", 60); /// 12) mapq
    kputc(tab, out);
    ksprintf(out, "s1:i:%d", map_score); /// chaining score
    kputc(tab, out);
    ksprintf(out, "NM:i:%d", aln_size - num_ident); /// gaps and mismatches in the alignment
    kputc(tab, out);
    ksprintf(out, "AS:i:%d", map_score); ///  dp score
    if (dump_cigar) {
        kputc(tab, out);
        print_paf_cigar(qoff, qend, qsize, qaln, saln, aln_size, out);
    }
    if (dump_md) {
        kputc(tab, out);
        dump_md_string(qaln, saln, aln_size, out);
    }
    kputc('\n', out);
}