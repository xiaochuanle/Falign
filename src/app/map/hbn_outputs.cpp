#include "hbn_outputs.hpp"

#include "../../corelib/frag_id.hpp"
#include "../../corelib/hbn_package_version.h"
#include "../../corelib/pdqsort.h"
#include "../../sw/hbn_traceback_aux.h"

#include <algorithm>
#include <sstream>

using namespace std;

static void
bam_dump_hdr_hd(sam_hdr_t* hdr)
{
    int r = sam_hdr_add_line(hdr, "HD", "VN", SAM_VERSION, "SO", "unknown", "GO", "query", nullptr);
    if (r) HBN_ERR("Fail at adding BAM HD header line");
}

static void
bam_dump_hdr_sq(const char* name, const int size, sam_hdr_t* hdr)
{
    char sizebuf[256];
    snprintf(sizebuf, 256, "%d", size);
    int r = sam_hdr_add_line(hdr, "SQ", "SN", name, "LN", sizebuf, nullptr);
    if (r) HBN_ERR("Fail at adding BAM SQ header line");
}

static void 
bam_dump_hdr_pg(const char* pg_id, const char* pg_name, const char* pg_version, int argc, char* argv[], sam_hdr_t* hdr)
{
    char* cl = stringify_argv(argc, argv);
    int r = sam_hdr_add_line(hdr, "PG", "ID", pg_id, "PN", pg_name, "VN", pg_version, "CL", cl, nullptr);
    free(cl);
    if (r) HBN_ERR("FAIL at adding BAM PG header");
}

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

static void init_hbn_bam_hdr(int argc, char* argv[], HbnUnpackedDatabase* subjects, sam_hdr_t* hdr)
{
    bam_dump_hdr_hd(hdr);

    int n = subjects->NumSeqs();
    for (int i = 0; i < n; ++i) {
        const char* name = subjects->SeqName(i);
        const int size = subjects->SeqSize(i);
        bam_dump_hdr_sq(name, size, hdr);
    }

    char pg_name[256];
    extract_pg_name(argv[0], pg_name);
    bam_dump_hdr_pg(pg_name, pg_name, HBN_PACKAGE_VERSION, argc, argv, hdr);
}

static void 
dump_cigar_string(const char* qaln, const char* saln, const int aln_size, ostringstream& os)
{
    const char* q = qaln;
    const char* s = saln;
    int i = 0;
    while (i < aln_size) {
        char type = 'N';
        int cnt = 0;
        if (q[i] == GAP_CHAR) {  // delete from subject (gap in query)
            type = 'D';
            while (i < aln_size && q[i] == GAP_CHAR) {
		if (s[i] == GAP_CHAR) dump_align_string(q, s, aln_size, stderr);
                hbn_assert(s[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else if (s[i] == GAP_CHAR) { // insert into subject (gap in subject)
            type = 'I';
            while (i < aln_size && s[i] == GAP_CHAR) {
                hbn_assert(q[i] != GAP_CHAR);
                ++i;
                ++cnt;
            }
        } else { // substitution
            type = 'M';
            hbn_assert(q[i] != GAP_CHAR && s[i] != GAP_CHAR);
            while (i < aln_size && q[i] != GAP_CHAR && s[i] != GAP_CHAR) {
                ++i;
                ++cnt;
            }
        }
        os << cnt << type;
    }
}

static void
dump_sam_cigar(const int qoff, const int qend, const int qsize,
    const char* qaln, const char* saln, const int aln_size, ostringstream& os)
{
    if (qoff) os << qoff << 'S';
    dump_cigar_string(qaln, saln, aln_size, os);
    if (qend < qsize) os << qsize - qend << 'S';
}

static BOOL 
validate_cigar_and_seq_size(const char* cigar, const int query_size)
{
    const int c_n = strlen(cigar);
    int c_i = 0;
    int q_base = 0;
    while (c_i < c_n) {
        size_t j = c_i + 1;
        while (j < c_n && isdigit(cigar[j])) ++j;
        char op = cigar[j];
        int num = atoi(cigar + c_i);
        switch (op)
        {
        case 'S':
            q_base += num;
            break;
        case 'M':
            q_base += num;
            break;
        case 'D':
            break;
        case 'I':
            q_base += num;
            break;
        default:
            break;
        }
        c_i = j + 1;
    }    
    if (q_base != query_size) {
        fprintf(stderr, "%s\n", cigar);
        HBN_LOG("cigar and seq_size is inconsistent: %d v.s. %d", q_base, query_size);
        return false;
    }
    return true;
}

static BOOL
truncate_align_bad_ends_4(const char* qaln,
    const char* saln,
    const int aln_size,
    int* qoff,
    int* qend,
    int* soff,
    int* send,
    const char** qas_,
    const char** qae_,
    const char** sas_,
    const char** sae_)
{
    const char* const qas = qaln;
    const char* const qae = qaln + aln_size;
    const char* const sas = saln;
    const char* const sae = saln + aln_size;
    int qcnt = 0, scnt = 0;
    const char* qa = qas;
    const char* sa = sas;
    int m = 0;
    const int M = 4;

    while (qa < qae) {
        int qc = *qa; ++qa;
        int sc = *sa; ++sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == M) break;
    }
    if (m < M || qa == qae) return 0;
    *qas_ = qa - M;
    *sas_ = sa - M;
    *qoff += qcnt - M;
    *soff += scnt - M;  

    qcnt = 0;
    scnt = 0;
    qa = qae;
    sa = sae;
    m = 0;
    while (qa > qas) {
        --qa;
        --sa;
        int qc = *qa;
        int sc = *sa;
        if (qc != '-') ++qcnt;
        if (sc != '-') ++scnt;
        m = (qc == sc) ? (m+1) : 0;
        if (m == M) break;
    }
    hbn_assert(m == M, "m = %d, qcnt = %d, scnt = %d", m, qcnt, scnt);
    if (qa == qas) return 0;

    *qae_ = qa + M;
    *sae_ = sa + M;
    *qend -= qcnt - M;
    *send -= scnt - M; 
    return TRUE;
}

static void
s_validate_aligned_string(const char* source_file,
						const char* source_func,
						const int source_line,
						int qid,
						const char* query,
						const int qoff,
						const int qend,
						const char* query_mapped_string,
						int tid,
						const u8* target,
						const int toff,
						const int tend,
						const char* target_mapped_string,
						const size_t align_size,
					    const BOOL right_extend)
{
	//return;
	int x = qoff, y = toff;
	for (size_t i = 0; i != align_size; ++i) {
		const char qc = query_mapped_string[i];
		if (qc != GAP_CHAR) {
            int c = right_extend ? query[x] : query[-x];
			const char qc1 = c;
            if (qc != qc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, qc = %c, qc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d, align_size = %lu\n",
					  source_file,
					  source_func,
					  source_line,
					  qid,
					  tid,
					  right_extend,
					  i,
					  x,
					  y,
					  qc,
					  qc1,
					  qoff,
					  qend,
					  toff,
					  tend,
					  align_size);
                abort();
            }		  
			++x;
		}
		const char tc = target_mapped_string[i];
		if (tc != GAP_CHAR) {
            int c = right_extend ? target[y] : target[-y];
			const char tc1 = DECODE_RESIDUE(c);
            if (tc != tc1) {
			    fprintf(stderr, "[%s, %s, %d] qid = %d, tid = %d, right_extend = %d, i = %lu, x = %d, y = %d, tc = %c, tc1 = %c, qoff = %d, qend = %d, toff = %d, tend = %d\n",
						  source_func,
						  source_func,
						  source_line,
						  qid,
						  tid,
						  right_extend,
						  i,
						  x,
						  y,
						  tc,
						  tc1,
						  qoff,
						  qend,
						  toff,
						  tend);
                dump_align_string(query_mapped_string, target_mapped_string, align_size, stderr);
                abort();
            }
			++y;
		}
	}
}

static bool
s_truncate_pca_align_bad_ends(HbnUnpackedDatabase* subjects,
    const int query_id,
    const char* fwd_query,
    const char* rev_query,
    const int query_size,
    TrimPcaList* tpca_list,
    TrimPca* tpca,
    const char** _qas,
    const char** _sas,
    int* _as_size)
{
    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_size = tpca->as_size;
    const char* qae = qas + as_size;
    const char* sae = sas + as_size;
    int qb = tpca->qoff;
    int qe = tpca->qend;
    int sb = tpca->soff;
    int se = tpca->send;

    if (!truncate_align_bad_ends_4(qas, sas, as_size, &qb, &qe, &sb, &se, &qas, &qae, &sas, &sae)) return false;
    as_size = qae - qas;

    const char* query = (tpca->qdir == FWD) ? fwd_query : rev_query;
    const u8* chr_seq = subjects->GetSequence(tpca->pca.sid);
    s_validate_aligned_string(__FILE__, __FUNCTION__, __LINE__,
        query_id, query, qb, qe, qas, 
        tpca->pca.sid, chr_seq, sb, se, sas, as_size, TRUE);

    tpca->qoff = qb;
    tpca->qend = qe;
    tpca->soff = sb;
    tpca->send = se;

    *_qas = qas;
    *_sas = sas;
    *_as_size = as_size;
    return true;
}

static void 
s_dump_secondary_pca(PoreCAlign* pca, const char* query_name, HbnUnpackedDatabase* subjects, ostringstream& os)
{
    os << "SA:Z:";
    const char* subject_name = subjects->SeqName(pca->sid);
    const char strand = (pca->qdir == FWD) ? '+' : '-';
    os << strand << ':' << pca->chain_qoff << ':' << pca->chain_qend
       << ':' << subject_name << ':' << pca->soff << ':' << pca->send << ':' << pca->pi;
}

static bam1_t*
build_one_read_bam(RestrictEnzymeLociList* reloci_list,
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList* tpca_list,
        TrimPca* tpca)
{
    uint16_t flag = 0; if (tpca->qdir == REV) flag |= 0x10;

    //const char* qas = nullptr;
    //const char* sas = nullptr;
    //int as_size = 0;
    //if (!s_truncate_pca_align_bad_ends(subjects, query_id, fwd_query, rev_query, query_size, tpca_list, tpca, &qas, &sas, &as_size)) return nullptr;

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_size = tpca->as_size;

    float pi = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    ostringstream os;
    dump_sam_cigar(tpca->qoff, tpca->qend, query_size, qas, sas, as_size, os);
    string cigar_s = os.str();
    validate_cigar_and_seq_size(cigar_s.c_str(), query_size);

    uint32_t* cigar_a = nullptr;
    size_t cigar_buffer_size = 0;
    auto num_cigar_op = sam_parse_cigar(cigar_s.c_str(), nullptr, &cigar_a, &cigar_buffer_size);
    if (num_cigar_op < 0) HBN_ERR("FAIL at parsing CIGAR");

    const char* query = (tpca->qdir == FWD) ? fwd_query : rev_query;
    const char* qv = nullptr;
    if (tpca->qdir == FWD && fwd_qv) qv = fwd_qv;
    if (tpca->qdir == REV && rev_qv) qv = rev_qv;

    bam1_t* bam = bam_init1();
    int r = bam_set1(bam, strlen(query_name), query_name, flag, 
        tpca->pca.sid, tpca->soff, tpca->pca.map_q, num_cigar_op, cigar_a,
        -1, -1, 0, query_size, query, qv, 512);
    if (r < 0) HBN_ERR("FAIL at building BAM record");
    free(cigar_a);

    char tagname[2];
    /// qid
    tagname[0] = 'q'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_id));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qdir
    tagname[0] = 'q'; tagname[1] = 'd';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.qdir));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qs
    tagname[0] = 'q'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qe
    tagname[0] = 'q'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqs
    tagname[0] = 'q'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqe
    tagname[0] = 'q'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ql
    tagname[0] = 'q'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_size));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sid
    tagname[0] = 's'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.sid));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ss
    tagname[0] = 's'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// se
    tagname[0] = 's'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vss
    tagname[0] = 'v'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vse
    tagname[0] = 'v'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sl
    tagname[0] = 's'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.ssize));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// identity
    tagname[0] = 'p'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'f', sizeof(float), (const uint8_t*)(void*)(&pi));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    if (tpca->pca.qc >= 0) {
        ostringstream os;
        s_dump_secondary_pca(all_pca_a + tpca->pca.qc, query_name, subjects, os);
        tagname[0] = 'S'; tagname[1] = 'A';
        r = bam_aux_update_str(bam, tagname, os.str().size(), os.str().c_str());
        if (r) HBN_ERR("Could not add '%s' tag to BAM record", os.str().c_str());
    }

    return bam;
}


static bam1_t*
build_one_frag_bam(RestrictEnzymeLociList* reloci_list,
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList* tpca_list,
        TrimPca* tpca)
{
    uint16_t flag = 0; if (tpca->qdir == REV) flag |= 0x10;

    //const char* qas = nullptr;
    //const char* sas = nullptr;
    //int as_size = 0;
    //if (!s_truncate_pca_align_bad_ends(subjects, query_id, fwd_query, rev_query, query_size, tpca_list, tpca, &qas, &sas, &as_size)) return nullptr;

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    int as_size = tpca->as_size;

    char frag_prefix[FRAD_ID_SIZE];
    frag_id_to_string(query_id, tpca->frag_id, tpca->pca.sid, tpca->soff, frag_prefix);
    ostringstream os;
    os << query_name << '_' << frag_prefix;
    string frag_name = os.str();

    float pi = calc_ident_perc(qas, sas, as_size, nullptr, nullptr);
    os.str("");
    dump_cigar_string(qas, sas, as_size, os);
    string cigar_s = os.str();
    validate_cigar_and_seq_size(cigar_s.c_str(), tpca->qend - tpca->qoff);

    uint32_t* cigar_a = nullptr;
    size_t cigar_buffer_size = 0;
    auto num_cigar_op = sam_parse_cigar(cigar_s.c_str(), nullptr, &cigar_a, &cigar_buffer_size);
    if (num_cigar_op < 0) HBN_ERR("FAIL at parsing CIGAR");

    const char* query = (tpca->qdir == FWD) ? fwd_query + tpca->qoff : rev_query + tpca->qoff;
    const char* qv = nullptr;
    if (tpca->qdir == FWD && fwd_qv) qv = fwd_qv + tpca->qoff;
    if (tpca->qdir == REV && rev_qv) qv = rev_qv + tpca->qoff;

    bam1_t* bam = bam_init1();
    int r = bam_set1(bam, frag_name.size(), frag_name.c_str(), flag, 
        tpca->pca.sid, tpca->soff, tpca->pca.map_q, num_cigar_op, cigar_a,
        -1, -1, 0, tpca->qend - tpca->qoff, query, qv, 512);
    if (r < 0) HBN_ERR("FAIL at building BAM record");
    free(cigar_a);

    char tagname[2];
    /// qid
    tagname[0] = 'q'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_id));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qdir
    tagname[0] = 'q'; tagname[1] = 'd';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.qdir));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qs
    tagname[0] = 'q'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// qe
    tagname[0] = 'q'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqs
    tagname[0] = 'q'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_qoff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vqe
    tagname[0] = 'q'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_qend));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ql
    tagname[0] = 'q'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&query_size));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sid
    tagname[0] = 's'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.sid));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// ss
    tagname[0] = 's'; tagname[1] = 's';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// se
    tagname[0] = 's'; tagname[1] = 'e';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vss
    tagname[0] = 'v'; tagname[1] = 'S';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_soff));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// vse
    tagname[0] = 'v'; tagname[1] = 'E';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.enzyme_send));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// sl
    tagname[0] = 's'; tagname[1] = 'l';
    r = bam_aux_append(bam, tagname, 'i', sizeof(int), (const uint8_t*)(void*)(&tpca->pca.ssize));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    /// identity
    tagname[0] = 'p'; tagname[1] = 'i';
    r = bam_aux_append(bam, tagname, 'f', sizeof(float), (const uint8_t*)(void*)(&pi));
    if (r) HBN_ERR("FAIL at appending tag %c%c to BAM record", tagname[0], tagname[1]);

    if (tpca->pca.qc >= 0) {
        ostringstream os;
        s_dump_secondary_pca(all_pca_a + tpca->pca.qc, query_name, subjects, os);
        tagname[0] = 'S'; tagname[1] = 'A';
        r = bam_aux_update_str(bam, tagname, os.str().size(), os.str().c_str());
        if (r) HBN_ERR("Could not add '%s' tag to BAM record", os.str().c_str());
    }

    return bam;
}

void
dump_one_pca(RestrictEnzymeLociList* reloci_list,
    HbnUnpackedDatabase* subjects,
    const char* query_name,
    const int query_id,
    const int query_size,
    PoreCAlign* all_pca_a,
    int all_pca_c,
    TrimPcaList* tpca_list,
    TrimPca* tpca,
    ostringstream& out)
{
    PoreCAlign* pca = &tpca->pca;
    const char* subject_name = subjects->SeqName(pca->sid);
    const int subject_size = subjects->SeqSize(pca->sid);

    out << query_name;
    out << '\t' << query_size;
    out << '\t' << tpca->fwd_qoff;
    out << '\t' << tpca->fwd_qend;
    out << '\t' << ((pca->qdir == FWD) ? '+' : '-');

    out << '\t' << subject_name;
    out << '\t' << subject_size;
    out << '\t' << tpca->soff;
    out << '\t' << tpca->send;

    const char* qas = tpca_list->align_strings.c_str() + tpca->qas_offset;
    const char* sas = tpca_list->align_strings.c_str() + tpca->sas_offset;
    const int as_size = tpca->as_size;
    int mat = 0;
    for (int i = 0; i < as_size; ++i) mat += (qas[i] == sas[i]);

    out << '\t' << mat;
    out << '\t' << as_size;
    out << '\t' << tpca->pca.map_q;

    if (pca->qdir == FWD) {
        out << '\t' << "qS:i:" << pca->enzyme_qoff;
    } else {
        int enzyme_size = reloci_list->enzyme.enzyme_size;
        int x = (pca->enzyme_qend == pca->qsize) ? 0 : (pca->qsize - pca->enzyme_qend - enzyme_size);
        out << '\t' << "qS:i:" << x;
    }

    if (pca->qdir == FWD) {
        out << '\t' << "qE:i:" << pca->enzyme_qend;
    } else {
        int enzyme_size = reloci_list->enzyme.enzyme_size;
        int x = (pca->enzyme_qoff == 0) ? pca->qsize : (pca->qsize - pca->enzyme_qoff - enzyme_size);
        out << '\t' << "qE:i:" << x;        
    }

    out << '\t' << "vS:i:" << pca->enzyme_soff;
    out << '\t' << "vE:i:" << pca->enzyme_send;

    double pi = calc_ident_perc(qas, sas, as_size, NULL, NULL);
    out << '\t' << "pi:f:" << pi;

    out << '\t' << "qi:i:" << query_id;
    
    if (pca->qc >= 0) {
        out << '\t';
        s_dump_secondary_pca(all_pca_a + pca->qc, query_name, subjects, out);
    }

    out << '\n';
}

void HbnOutputs::x_init_sam_hdr(int argc, char* argv[], HbnUnpackedDatabase* reference)
{
    M_sam_hdr = sam_hdr_init();
    init_hbn_bam_hdr(argc, argv, reference, M_sam_hdr);
    int r = sam_hdr_write(M_sam_out, M_sam_hdr);
    if (r) HBN_ERR("FAIL at writing BAM header to file");
}

void HbnOutputs::dump(RestrictEnzymeLociList* reloci_list, 
        HbnUnpackedDatabase* subjects,
        const char* query_name,
        const int query_id,
        const char* fwd_query,
        const char* rev_query,
        const char* fwd_qv,
        const char* rev_qv,
        const int query_size,
        PoreCAlign* all_pca_a,
        int all_pca_c,
        TrimPcaList& pca_list)
{
    TrimPca* a = pca_list.tpca_list.data();
    int c = pca_list.tpca_list.size();
    int n = 0;
    for (int i = 0; i < c; ++i) if (a[i].is_valid) a[n++] = a[i];
    c = n;
    for (int i = 0; i < c; ++i) {
        if (a[i].qdir == FWD) {
            a[i].fwd_qoff = a[i].qoff;
            a[i].fwd_qend = a[i].qend;
        } else {
            a[i].fwd_qoff = query_size - a[i].qend;
            a[i].fwd_qend = query_size - a[i].qoff;
        }
    }
    pdqsort(a, a + c, [](const TrimPca& x, const TrimPca& y) { return x.fwd_qoff < y.fwd_qoff; });
    for (int i = 0; i < c; ++i) a[i].frag_id = i;

    if (M_outfmt == eOutputFmt_BAM || M_outfmt == eOutputFmt_SAM) {
        vector<bam1_t*> bam_list;
        bam1_t* bam;
        for (int i = 0; i < c; ++i) {
            bam = build_one_read_bam(reloci_list, subjects, query_name, query_id,
                fwd_query, rev_query, fwd_qv, rev_qv, query_size, all_pca_a, all_pca_c, &pca_list, a + i);
            if (!bam) continue;
            bam_list.push_back(bam);
        }
        lock_guard<mutex> lg(M_out_mutex);
        for (auto b : bam_list) {
            int r = sam_write1(M_sam_out, M_sam_hdr, b);
            if (r < 0) HBN_ERR("FAIL at writing BAM record");
            bam_destroy1(b);
        }
    } else if (M_outfmt == eOutputFmt_FragBAM || M_outfmt == eOutputFmt_FragSAM) {
        vector<bam1_t*> bam_list;
        bam1_t* bam;
        for (int i = 0; i < c; ++i) {
            bam = build_one_frag_bam(reloci_list, subjects, query_name, query_id,
                fwd_query, rev_query, fwd_qv, rev_qv, query_size, all_pca_a, all_pca_c, &pca_list, a + i);
            if (!bam) continue;
            bam_list.push_back(bam);
        }
        lock_guard<mutex> lg(M_out_mutex);
        for (auto b : bam_list) {
            int r = sam_write1(M_sam_out, M_sam_hdr, b);
            if (r < 0) HBN_ERR("FAIL at writing BAM record");
            bam_destroy1(b);
        }        
    } else if (M_outfmt == eOutputFmt_PAF) {
        ostringstream os;
        for (int i = 0; i < c; ++i) {
            dump_one_pca(reloci_list, subjects, query_name, query_id, query_size, all_pca_a, all_pca_c, &pca_list, a+i, os);
        }
        lock_guard<mutex> lg(M_out_mutex);
        hbn_fwrite(os.str().c_str(), 1, os.str().size(), M_out);
    }
}