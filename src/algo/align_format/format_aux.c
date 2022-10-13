#include "format_aux.h"

#include "../hbn_traceback_aux.h"

void dump_cigar_string(const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    const char* q = qaln;
    const char* s = saln;
    int i = 0;
    while (i < aln_size) {
        char type = 'N';
        int cnt = 0;
        if (q[i] == GAP_CHAR) {  // delete from subject (gap in query)
            type = 'D';
            while (i < aln_size && q[i] == GAP_CHAR) {
		if (s[i] == GAP_CHAR) dump_align_string(q, s, aln_size, stderr);
                hbn_assert(s[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else if (s[i] == GAP_CHAR) { // insert into subject (gap in subject)
            type = 'I';
            while (i < aln_size && s[i] == GAP_CHAR) {
                hbn_assert(q[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else { // substitution
            type = 'M';
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (i < aln_size && q[i] != GAP_CHAR && s[i] != GAP_CHAR) {
                ++i;
                ++cnt;
            }
        }
        ksprintf(out, "%d%c", cnt, type);
    }
}

void dump_md_string(const char* qaln, const char* saln, const int aln_size, kstring_t* out)
{
    ksprintf(out, "MD:Z:");
    const char* q = qaln;
    const char* s = saln;
    int i = 0, l_md = 0;
    while (i < aln_size) {
        int j = i;
        if (s[i] == GAP_CHAR) { // insert into subject (gap in  query)
            while (j < aln_size && s[j] == GAP_CHAR) {
                hbn_assert(q[j] != GAP_CHAR);
                ++j;
            }
        } else if (q[i] == GAP_CHAR) { // delete from subject (gap in subject)
            while (j < aln_size && q[j] == GAP_CHAR) {
                hbn_assert(s[j] != GAP_CHAR);
                ++j;
            }
            ksprintf(out, "%d^", l_md);
            kputsn(&s[i], j - i, out);
            l_md = 0;
        } else { // substitution
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (j < aln_size && q[j] != GAP_CHAR && s[j] != GAP_CHAR) ++j;
            for (int k = i; k < j; ++k) {
                if (q[k] != s[k]) {
                    ksprintf(out, "%d%c", l_md, s[k]);
                    l_md = 0;
                } else {
                    ++l_md;
                }
            }
        }
        i = j;
    }   
    if (l_md > 0) ksprintf(out, "%d", l_md); 
}
