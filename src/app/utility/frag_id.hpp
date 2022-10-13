#ifndef __READ_FRAG_ID_HPP
#define __READ_FRAG_ID_HPP

#include <map>
#include <string>

#define FRAD_ID_SIZE 128

void
frag_id_to_string(const int read_id, const int frag_id, const int subject_id, const int subject_offset, char buf[]);

void
extract_frag_id_from_name(const char* s, const int sl, int* read_id, int* frag_id, int* subject_id, int* subject_offset);

#endif // __READ_FRAG_ID_HPP