#include "split_string_by_char.hpp"

using namespace std;

static inline int
s_string_find_first_not_of(const char* s, const int sl, const char delim, int pos)
{
    for (int i = pos; i < sl; ++i) {
        if (s[i] != delim) return i;
    }
    return sl;
}

static inline int
s_string_find_first_of(const char* s, const int sl, const char delim, int pos)
{
    for (int i = pos; i < sl; ++i) {
        if (s[i] == delim) return i;
    }
    return sl;
}

void
split_string_by_char(const char* s, const int sl, const char delim, std::vector<std::pair<const char*, int>>& tokens)
{
    tokens.clear();
    int last_pos = s_string_find_first_not_of(s, sl, delim, 0);
    int pos = s_string_find_first_of(s, sl, delim, last_pos);
    while (last_pos < sl) {
        tokens.push_back(pair<const char*, int>(s + last_pos, pos - last_pos));
        last_pos = s_string_find_first_not_of(s, sl, delim, pos);
        pos = s_string_find_first_of(s, sl, delim, last_pos);
    }
}