#ifndef __SMALL_OBJECT_ALLOC_H
#define __SMALL_OBJECT_ALLOC_H

#include "hbn_aux.h"
#include "kvec.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    char* data;
    u32 avail_data;
    u32 object_size;
} MemoryChunk;

typedef kvec_t(MemoryChunk*) vec_mem_chunk;

typedef struct {
    vec_mem_chunk   chunk_list;
    u32             first_avail_chunk;
    u32             object_size;
    void**          deallocated_object_list;
    int             deallocated_objects;
	int		num_alloc;
} SmallObjectAlloc;

SmallObjectAlloc*
SmallObjectAllocNew(const u32 object_size);

SmallObjectAlloc*
SmallObjectAllocFree(SmallObjectAlloc* alloc);

void
SmallObjectAllocClear(SmallObjectAlloc* alloc);

void*
SmallObjectAllocAlloc(SmallObjectAlloc* alloc, const u32 num_objects);

void*
SmallObjectAllocAllocOne(SmallObjectAlloc* alloc);

void
SmallocObjactAllocDeallocOne(SmallObjectAlloc* alloc, void* object);

void
SmallObjectAllocStat(SmallObjectAlloc* alloc);

#ifdef __cplusplus
}
#endif

#endif // __SMALL_OBJECT_ALLOC_H
