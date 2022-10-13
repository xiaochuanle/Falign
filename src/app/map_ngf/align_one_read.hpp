#ifndef __ALIGN_ONE_READ_HPP
#define __ALIGN_ONE_READ_HPP

#include "map_one_volume.hpp"

void
align_one_read(MapThreadData* data, const int qidx, kstring_t* out);

void
align_one_read_1(MapThreadData* data, const int qidx, kstring_t* out);

#endif // __ALIGN_ONE_READ_HPP