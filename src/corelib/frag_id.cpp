#include "frag_id.hpp"

#include "hbn_aux.h"

#include <cctype>
#include <cstdlib>

using namespace std;

static const int kReadIdSize = 10;
static const int kFragIdSize = 3;
static const int kSubjectIdSize = 10;
static const int kSubjectOffsetSize = 10;

void
frag_id_to_string(const int read_id, const int frag_id, const int subject_id, const int subject_offset, char buf[])
{
    char buf1[64], buf2[64], buf3[64], buf4[64];
    u64_to_fixed_width_string_r(read_id, buf1, kReadIdSize);
    u64_to_fixed_width_string_r(frag_id, buf2, kFragIdSize);
    u64_to_fixed_width_string_r(subject_id, buf3, kSubjectIdSize);
    u64_to_fixed_width_string_r(subject_offset, buf4, kSubjectOffsetSize);
    sprintf(buf, "%s:%s:%s:%s", buf1, buf2, buf3, buf4);
}

static void
s_illegal_frag_name_format(const int line, const char* s, const int sl)
{
    int read_id = 1, frag_id = 2, subject_id = 3, subject_offset = 4;
    char buf[FRAD_ID_SIZE];
    frag_id_to_string(read_id, frag_id, subject_id, subject_offset, buf);
    fprintf(stderr, "[%d] Illegal frag name format '", line);
    hbn_fwrite(s, 1, sl, stderr);
    fprintf(stderr, "'\n");
    fprintf(stderr, "A legal frag name with read_id = %d, frag_id = %d, subject_id = %d and subject_offset = %d should end with '%s'\n", 
        read_id, frag_id, subject_id, subject_offset, buf);
    abort();
}

void
extract_frag_id_from_name(const char* s, const int sl, int* read_id, int* frag_id, int* subject_id, int* subject_offset)
{
    int p = sl;
    while (p) {
        --p;
        if (s[p] == '_') break;
    }
    if (s[p] != '_') s_illegal_frag_name_format(__LINE__, s, sl);
    ++p;
    if (p >= sl) s_illegal_frag_name_format(__LINE__, s, sl);
    if (!isdigit(s[p])) s_illegal_frag_name_format(__LINE__, s, sl);
    if (read_id) *read_id = atoi(s + p);

    while (p < sl && s[p] != ':') ++p;
    if (s[p] != ':') s_illegal_frag_name_format(__LINE__, s, sl);
    ++p;
    if (p >= sl) s_illegal_frag_name_format(__LINE__, s, sl);
    if (!isdigit(s[p])) s_illegal_frag_name_format(__LINE__, s, sl);
    if (frag_id) *frag_id = atoi(s + p);

    while (p < sl && s[p] != ':') ++p;
    if (s[p] != ':') s_illegal_frag_name_format(__LINE__, s, sl);
    ++p;
    if (p >= sl) s_illegal_frag_name_format(__LINE__, s, sl);
    if (!isdigit(s[p])) s_illegal_frag_name_format(__LINE__, s, sl);
    if (subject_id) *subject_id = atoi(s + p);

    while (p < sl && s[p] != ':') ++p;
    if (s[p] != ':') s_illegal_frag_name_format(__LINE__, s, sl);
    ++p;
    if (p >= sl) s_illegal_frag_name_format(__LINE__, s, sl);
    if (!isdigit(s[p])) s_illegal_frag_name_format(__LINE__, s, sl);
    if (subject_offset) *subject_offset = atoi(s + p);
}