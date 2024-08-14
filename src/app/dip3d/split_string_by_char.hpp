#ifndef __SPLIT_STRING_BY_CHAR
#define __SPLIT_STRING_BY_CHAR

#include <utility>
#include <vector>

void
split_string_by_char(const char* s, const int sl, const char delim, std::vector<std::pair<const char*, int>>& tokens);

#endif // __SPLIT_STRING_BY_CHAR