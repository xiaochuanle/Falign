#include "necat_format.h"

#include "../../corelib/hbn_package_version.h"

static void
extract_pg_name(const char* argv_0, char pg_name[])
{
    int n = strlen(argv_0);
    int s = n;
    while (s) {
        --s;
        if (argv_0[s] == '/') {
            ++s;
            break;
        }
    }

    int i = 0;
    for (; s < n; ++i, ++s) pg_name[i] = argv_0[s];
    pg_name[i] = '\0';
}

void 
dump_sam_hdr(const char* rg_id,
    const char* rg_sample,
    const CSeqDB* db,
    RawReadReader* rdb,
    int argc,
    char* argv[],
    FILE* out)
{
    sam_dump_hdr_hd(out);

    if (db) {
        for (int i = 0; i < db->dbinfo.num_seqs; ++i) {
            const char* subject_name = seqdb_seq_name(db, i);
            int subject_size = seqdb_seq_size(db, i);
            sam_dump_hdr_sq(subject_name, subject_size, out);
        }
    } else {
        hbn_assert(rdb);
        for (int i = 0; i < rdb->dbinfo.num_seqs; ++i) {
            const char* subject_name = RawReadReader_ReadName(rdb, i);
            int subject_size = RawReadReader_ReadSize(rdb, i);
            sam_dump_hdr_sq(subject_name, subject_size, out);
        }
    }

    if (rg_id) sam_dump_hdr_rg(rg_id, rg_sample, out);

    char pg_name[256];
    extract_pg_name(argv[0], pg_name);
    sam_dump_hdr_pg(pg_name, pg_name, HBN_PACKAGE_VERSION, argc, argv, out);
}