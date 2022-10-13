#ifndef __SEQ_LOADER_H
#define __SEQ_LOADER_H

#include "../corelib/cstr_util.h"
#include "../corelib/hbn_aux.h"
#include "../corelib/fasta.h"

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct {
    size_t name_offset;
    int seq_size;
    size_t seq_offset;
    size_t qv_offset;
} SeqInfo;

typedef kvec_t(SeqInfo) vec_seq_info;

typedef struct {
    vec_u8 fwd_unpacked_seq;
    vec_u8 rev_unpacked_seq;
    size_t seq_size;
    size_t max_seq_size;
    vec_seq_info seq_info_list;
    kstring_t name_list;
    kstring_t qv_list;
    int num_seqs;

    HbnFastaReader* fasta_reader;
    int file_name_idx;
    vec_int file_name_offset_list;
    kstring_t file_name_list;
    int load_in_terms_of_file;
} SeqReader;

const u8*
SeqReader_Seq(SeqReader* reader, int id, int strand);

const char*
SeqReader_QV(SeqReader* reader, int id);

int
SeqReader_SeqSize(SeqReader* reader, int id);

size_t 
SeqReader_SeqOffset(SeqReader* reader, int id);

const char*
SeqReader_SeqName(SeqReader* reader, int id);

int 
SeqReader_NumSeqs(SeqReader* reader);

size_t
SeqReaderLoad(SeqReader* reader, size_t volume_size, BOOL load_rev_seq);

BOOL 
SeqReader_OpenNextFile(SeqReader* reader);

const char*
SeqReader_CurrentFileName(SeqReader* reader);

void
SeqReaderClear(SeqReader* reader);

SeqReader*
SeqReaderNew(int argc, char* argv[]);

SeqReader*
SeqReaderFree(SeqReader* reader);

#ifdef __cplusplus
}
#endif 

#endif // __SEQ_LOADER_H