#include "seq_loader.h"
#include "../corelib/db_format.h"

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

static BOOL 
s_extract_read_list_from_dir(const char* dir_path, kstring_t* name_list, vec_int* name_offset_list)
{
    struct stat statbuf;
    if (stat(dir_path, &statbuf)) {
        HBN_ERR("Fail to stat path '%s': %s", dir_path, strerror(errno));
    }
    if (!S_ISDIR(statbuf.st_mode)) return FALSE;

    DIR* dp = NULL;
    struct dirent* entry = NULL;
    if ((dp = opendir(dir_path)) == NULL) {
        HBN_ERR("Fail to open directory '%s': %s", dir_path, strerror(errno));
    }
    struct stat file_stat;
    char path[HBN_MAX_PATH_LEN];
    while ((entry = readdir(dp)) != NULL) {
        sprintf(path, "%s/%s", dir_path, entry->d_name);
        if (stat(path, &file_stat)) HBN_ERR("Fail to stat path '%s': %s", path, strerror(errno));
        if (S_ISDIR(file_stat.st_mode)) continue;
        int offset = ks_size(*name_list);
        kv_push(int, *name_offset_list, offset);
        kputsn(path, strlen(path), name_list);
        kputc('\0', name_list);
    }
    closedir(dp);
    return TRUE;
}

static void
s_extract_read_list_from_list_file(const char* file_list_path, kstring_t* name_list, vec_int* name_offset_list)
{
    HbnLineReader* line_reader = HbnLineReaderNew(file_list_path);
    while (!HbnLineReaderAtEof(line_reader)) {
        HbnLineReaderReadOneLine(line_reader);
        if (ks_empty(line_reader->line)) continue;
        kputc('\0', &line_reader->line);
        if (access(ks_s(line_reader->line), F_OK)) {
            HBN_LOG("Sequence file '%s' does not exist\n", ks_s(line_reader->line));
            abort();
        }
        int offset = ks_size(*name_list);
        kv_push(int, *name_offset_list, offset);
        kputsn(ks_s(line_reader->line), ks_size(line_reader->line), name_list);
    }
    HbnLineReaderFree(line_reader);
}

static void s_make_read_list(int argc, char* argv[], kstring_t* name_list, vec_int* name_offset_list)
{
    for (int i = 0; i < argc; ++i) {
        if (s_extract_read_list_from_dir(argv[i], name_list, name_offset_list)) continue;

        if (access(argv[i], F_OK)) {
            HBN_LOG("Sequence file '%s' does not exist.\n", argv[i]);
            abort();
        }

        EDbFormat format = hbn_guess_db_format(argv[i]);
        if (format == eDbFormatEmptyFile) {
            HBN_LOG("File '%s' contains no sequence data, skip it.\n", argv[i]);
            continue;
        }

        if (format == eDbFormatFasta || format == eDbFormatFastq) {
            int offset = ks_size(*name_list);
            kv_push(int, *name_offset_list, offset);
            kputsn(argv[i], strlen(argv[i]), name_list);
            kputc('\0', name_list);
        } else if (format == eDbFormatUnknown) {
            s_extract_read_list_from_list_file(argv[i], name_list, name_offset_list);
        }
    }
}

SeqReader*
SeqReaderFree(SeqReader* reader)
{
    if (!reader) return NULL;
    kv_destroy(reader->fwd_unpacked_seq);
    kv_destroy(reader->rev_unpacked_seq);
    ks_destroy(reader->name_list);
    ks_destroy(reader->qv_list);
    kv_destroy(reader->seq_info_list);

    if (reader->fasta_reader) reader->fasta_reader = HbnFastaReaderFree(reader->fasta_reader);
    kv_destroy(reader->file_name_offset_list);
    ks_destroy(reader->file_name_list);
    free(reader);
    return NULL;
}

SeqReader*
SeqReaderNew(int argc, char* argv[])
{
    SeqReader* reader = (SeqReader*)calloc(1, sizeof(SeqReader));
    kv_init(reader->fwd_unpacked_seq);
    kv_init(reader->rev_unpacked_seq);
    reader->num_seqs = 0;
    ks_init(reader->name_list);
    ks_init(reader->qv_list);
    kv_init(reader->seq_info_list);

    reader->fasta_reader = NULL;
    reader->file_name_idx = -1;
    kv_init(reader->file_name_offset_list);
    ks_init(reader->file_name_list);

    s_make_read_list(argc, argv, &reader->file_name_list, &reader->file_name_offset_list);   
    if (kv_empty(reader->file_name_offset_list)) {
        SeqReaderFree(reader);
        return NULL;
    }

    reader->load_in_terms_of_file = 0;
    SeqReader_OpenNextFile(reader);
    return reader;
}

BOOL 
SeqReader_OpenNextFile(SeqReader* reader)
{
    ++reader->file_name_idx;
    if (reader->file_name_idx >= kv_size(reader->file_name_offset_list)) return FALSE;
    const char* name = ks_s(reader->file_name_list) + kv_A(reader->file_name_offset_list, reader->file_name_idx);
    if (reader->fasta_reader) reader->fasta_reader = HbnFastaReaderFree(reader->fasta_reader);
    reader->fasta_reader = HbnFastaReaderNew(name);
    return TRUE;
}

const char*
SeqReader_CurrentFileName(SeqReader* reader)
{
    hbn_assert (reader->file_name_idx < kv_size(reader->file_name_offset_list));
    const char* name = ks_s(reader->file_name_list) + kv_A(reader->file_name_offset_list, reader->file_name_idx);
    return name;
}

void
SeqReaderClear(SeqReader* reader)
{
    kv_clear(reader->fwd_unpacked_seq);
    kv_clear(reader->rev_unpacked_seq);
    reader->num_seqs = 0;
    kv_clear(reader->seq_info_list);
    ks_clear(reader->name_list);
    ks_clear(reader->qv_list);
}

size_t
SeqReaderLoad(SeqReader* reader, size_t volume_size, BOOL load_rev_seq)
{
    SeqReaderClear(reader);
    size_t load_size = 0;
    while (load_size < volume_size) {
        if (HbnLineReaderAtEof(reader->fasta_reader->line_reader)) {
            if (reader->load_in_terms_of_file) break;
            if (!SeqReader_OpenNextFile(reader)) break;
        }
        HbnFastaReaderReadOneSeq(reader->fasta_reader);
        if (ks_empty(reader->fasta_reader->sequence)) continue;

        SeqInfo seq_info;
        seq_info.seq_size = ks_size(reader->fasta_reader->sequence);
        seq_info.seq_offset = kv_size(reader->fwd_unpacked_seq);
        seq_info.name_offset = ks_size(reader->name_list);
        seq_info.qv_offset = U64_MAX;
        kputsn(ks_s(reader->fasta_reader->name), ks_size(reader->fasta_reader->name), &reader->name_list);
        kputc('\0', &reader->name_list);
        if (!ks_empty(reader->fasta_reader->qual)) {
            seq_info.qv_offset = ks_size(reader->qv_list);
            kputsn(ks_s(reader->fasta_reader->qual), ks_size(reader->fasta_reader->qual), &reader->qv_list);
        }
        kv_push(SeqInfo, reader->seq_info_list, seq_info);
        load_size += seq_info.seq_size;

        const char* q = ks_s(reader->fasta_reader->sequence);
        for (int i = 0; i < seq_info.seq_size; ++i) {
            u8 c = q[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            kv_push(u8, reader->fwd_unpacked_seq, c);
        }
        if (!load_rev_seq) continue;
        for (int i = seq_info.seq_size - 1; i >= 0; --i) {
            u8 c = q[i];
            c = nst_nt4_table[c];
            if (c > 3) c = 0;
            c = 3 - c;
            kv_push(u8, reader->rev_unpacked_seq, c);
        }
    }
    if (!load_size) return load_size;

    reader->num_seqs = kv_size(reader->seq_info_list);
    char size_buf[64];
    u64_to_string_datasize(load_size, size_buf);
    HBN_LOG("Load %d sequences (%s)", reader->num_seqs, size_buf);
    return load_size;
}

int 
SeqReader_NumSeqs(SeqReader* reader)
{
    return reader->num_seqs;
}

const char*
SeqReader_SeqName(SeqReader* reader, int id)
{
    hbn_assert(id >= 0 && id < reader->num_seqs);
    return ks_s(reader->name_list) + kv_A(reader->seq_info_list, id).name_offset;
}

int
SeqReader_SeqSize(SeqReader* reader, int id)
{
    hbn_assert(id >= 0 && id < reader->num_seqs);
    return kv_A(reader->seq_info_list, id).seq_size;   
}

size_t 
SeqReader_SeqOffset(SeqReader* reader, int id)
{
    hbn_assert(id >= 0 && id < reader->num_seqs);
    return kv_A(reader->seq_info_list, id).seq_offset; 
}

const char*
SeqReader_QV(SeqReader* reader, int id)
{
    hbn_assert(id >= 0 && id < reader->num_seqs);
    if (kv_A(reader->seq_info_list, id).qv_offset == U64_MAX) return NULL;
    return ks_s(reader->qv_list) + kv_A(reader->seq_info_list, id).qv_offset;
}

const u8*
SeqReader_Seq(SeqReader* reader, int id, int strand)
{
    hbn_assert(id >= 0 && id < reader->num_seqs);
    if (strand == FWD) {
        if (kv_empty(reader->fwd_unpacked_seq)) {
            HBN_ERR("No forward strand sequence loaded");
        }
        return kv_data(reader->fwd_unpacked_seq) + kv_A(reader->seq_info_list, id).seq_offset;
    }
    hbn_assert(strand == REV);
    if (kv_empty(reader->rev_unpacked_seq)) {
        HBN_ERR("No reverse strand sequence loaded");
    }
    return kv_data(reader->rev_unpacked_seq) + kv_A(reader->seq_info_list, id).seq_offset;    
}