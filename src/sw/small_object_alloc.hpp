#ifndef __SMALL_OBJECT_ALLOC_HPP
#define __SMALL_OBJECT_ALLOC_HPP

#include "../corelib/hbn_aux.h"

#include <vector>

typedef struct {
    char* data;
    u32 avail_data;
    u32 object_size;
} MemoryChunk;

typedef struct {
    std::vector<MemoryChunk*>   chunk_list;
    u32             first_avail_chunk;
    u32             object_size;
} SmallObjectAlloc;

SmallObjectAlloc*
SmallObjectAllocNew(const u32 object_size);

SmallObjectAlloc*
SmallObjectAllocFree(SmallObjectAlloc* alloc);

void
SmallObjectAllocClear(SmallObjectAlloc* alloc);

void*
SmallObjectAllocAlloc(SmallObjectAlloc* alloc, const u32 num_objects);

void
SmallObjectAllocSetup(SmallObjectAlloc* alloc);

#endif // __SMALL_OBJECT_ALLOC_HPP