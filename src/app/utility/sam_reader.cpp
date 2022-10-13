#include "sam_reader.hpp"
#include "../../algo/hbn_traceback_aux.h"

using namespace std;

static const int kReadIdSize = 10;
static const int kFragIdSize = 3;

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
split_string_by_char(const char* s, const int sl, const char delim, vector<pair<const char*, int>>& tokens)
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

static void
s_parse_cigar(const char* query, const int query_size, const char* cigar, int cigar_size, const char* md, int md_size,
    string& qas, string& sas,
    int* _qb, int* _qe, int* _sb, int* _se, double* _pi)
{
    int qb = 0, qe = 0, sb = 0, se = 0;
    int i = 0;
    while (i < cigar_size && isdigit(cigar[i])) ++i;
    hbn_assert(i < cigar_size);
    if (cigar[i] == 'S') {
        qb = atoi(cigar);
        ++i;
    } else {
        i = 0;
    }
    qe = qb;

    qas.clear();
    sas.clear();
    while (i < cigar_size) {
        hbn_assert(isdigit(cigar[i]));
        int num = atoi(cigar + i);
        while (i < cigar_size && isdigit(cigar[i])) ++i;
        hbn_assert(i < cigar_size);
        char op = cigar[i];
        ++i;
        if (op == 'M') {
            for (int k = 0; k < num; ++k) {
		hbn_assert(qe + k < query_size, "qe = %d, k = %d, query_size = %d", qe, k, query_size);
                qas += query[qe + k];
                sas += query[qe + k];
            }
            qe += num;
            se += num;
        } else if (op == 'I') {
            for (int k = 0; k < num; ++k) {
		hbn_assert(qe + k < query_size);
                qas += query[qe + k];
                sas += GAP_CHAR;
            }
            qe += num;
        } else if (op == 'D') {
            for (int k = 0; k < num; ++k) {
                qas += GAP_CHAR;
                sas += 'N';
            }
            se += num;
        } else if (op == 'S') {
            continue;
        } else {
            HBN_ERR("Illegal CIGAR operation '%c'", op);
        }
    }

    //HBN_LOG("recover from md");
    hbn_assert(qas.size() == sas.size());
    //const size_t as_size = qas.size();
    //dump_align_string(qas.c_str(), sas.c_str(), as_size, stderr);
    int md_i = 0;
    int dist = 0;
    char op = 0;
    int as_i = 0;
    while (1) {
        dist = atoi(md + md_i);
        while (md_i < md_size && isdigit(md[md_i])) ++md_i;
        if (md_i >= md_size) break;
        op = md[md_i];
        if (op == '^') ++md_i;
        int cnt = 0;
        while (cnt < dist) {
            if (sas[as_i] != GAP_CHAR) ++cnt;
            ++as_i;
        }
        while (md_i < md_size && (!isdigit(md[md_i]))) {
            if (sas[as_i] == GAP_CHAR) {
                ++as_i;
                continue;
            }
            if (op == '^') hbn_assert(sas[as_i] == 'N');
            //fprintf(stderr, "set sas[%d] from %c to %c\n", as_i, sas[as_i], md[md_i]);
            sas[as_i] = md[md_i];
            ++as_i;
            ++md_i;
        }
    }
    for (auto c : sas) hbn_assert(c != 'N');
for (size_t i = 0; i < qas.size(); ++i) {
	char c = qas[i];
	hbn_assert(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == '-');
	c = sas[i];
	hbn_assert(c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == '-');
}

    *_qb = qb;
    *_qe = qe;
    *_sb = sb;
    *_se = se;
    *_pi = calc_ident_perc(qas.c_str(), sas.c_str(), qas.size(), NULL, NULL);
}

static const char*
s_extract_md(const std::string& sam, int& md_size)
{
    auto pos = sam.find("MD:Z:");
    if (pos == string::npos) {
        fprintf(stderr, "MD string is not presented in sam record\n");
        fprintf(stderr, "%s\n", sam.c_str());
        abort();
    }
    pos += 5;
    md_size = 0;
    size_t n = sam.size();
    while (pos + md_size < n) {
        if (isspace(sam[pos + md_size])) break;
        ++md_size;
    }
    return sam.c_str() + pos;
}

//////////////////

SamLoader::SamLoader(const char* sam_path)
{
    m_line_reader = HbnLineReaderNew(sam_path);
    while (!HbnLineReaderAtEof(m_line_reader)) {
        HbnLineReaderReadOneLine(m_line_reader);
        if (ks_A(m_line_reader->line, 0) != '@') {
            HbnLineReaderUngetline(m_line_reader);
            break;
        }
        m_sam_hdr_list.emplace_back(string(ks_s(m_line_reader->line), ks_size(m_line_reader->line)));
    }
    m_unget_frag_align = false;
}

SamLoader::~SamLoader()
{
    if (m_line_reader) m_line_reader = HbnLineReaderFree(m_line_reader);
}

void SamLoader::unget_frag_align()
{
    m_unget_frag_align = true;
}

static void
s_make_rc_read(const string& src, string& dst)
{
    dst.clear();
    dst.insert(dst.end(), src.rbegin(), src.rend());
    for (auto& c : dst) {
        int e = c;
hbn_assert(e == 'A' || e == 'C' || e == 'G' || e == 'T');
        e = nst_nt4_table[e];
        e = 3 - e;
        e = DECODE_RESIDUE(e);
        c = e;
    }
}

bool SamLoader::load_next_align()
{
    if (m_unget_frag_align) {
        m_unget_frag_align = false;
        return true;
    }

    if (HbnLineReaderAtEof(m_line_reader)) return false;
    HbnLineReaderReadOneLine(m_line_reader);
    m_line.assign(ks_s(m_line_reader->line), ks_size(m_line_reader->line));
    m_columns.clear();
    split_string_by_char(m_line.c_str(), m_line.size(), '\t', m_columns);

    FragmentAlignment align;
    m_frag_align.qname.assign(m_columns[0].first, m_columns[0].second);
    int flag = atoi(m_columns[1].first);
    m_frag_align.qdir = (flag & 0x10) ? REV : FWD;
    m_frag_align.sname.assign(m_columns[2].first, m_columns[2].second);
    if (m_columns[9].first[0] != '*') {
        if (m_frag_align.qdir == FWD) {
            m_fwd_read.assign(m_columns[9].first, m_columns[9].second);
            s_make_rc_read(m_fwd_read, m_rev_read);
        } else {
            m_rev_read.assign(m_columns[9].first, m_columns[9].second);
            s_make_rc_read(m_rev_read, m_fwd_read);
        }
    }
    const string& read = (m_frag_align.qdir == FWD) ? m_fwd_read : m_rev_read;
    int md_size = 0;
    const char* md = s_extract_md(m_line, md_size);
    s_parse_cigar(read.c_str(), read.size(), m_columns[5].first, m_columns[5].second, md, md_size,
        m_qas, m_sas, &m_frag_align.qoff, &m_frag_align.qend, &m_frag_align.soff, &m_frag_align.send, &m_frag_align.pi);
    m_frag_align.qsize = read.size();
    int soff = atoi(m_columns[3].first);
    hbn_assert(soff > 0);
    --soff;
    m_frag_align.soff += soff;
    m_frag_align.send += soff;
    m_frag_align.as_size = m_qas.size();

    if (m_frag_align.qdir == REV) {
        int x = m_frag_align.qsize - m_frag_align.qend;
        int y = m_frag_align.qsize - m_frag_align.qoff;
        m_frag_align.qoff = x;
        m_frag_align.qend = y;
    }
    m_frag_align.map_q = atoi(m_columns[4].first);

#if 0
    fprintf(stderr, "[%s, %d, %d, %d, %d] x [%s, %d, %d], %g\n", m_frag_align.qname.c_str(), m_frag_align.qdir,
        m_frag_align.qoff, m_frag_align.qend, m_frag_align.qsize,
        m_frag_align.sname.c_str(), m_frag_align.soff, m_frag_align.send, m_frag_align.pi);
#endif
    return true;
}

int SamLoader::load_next_align_list(std::vector<FragmentAlignment>& frag_align_list,
            std::vector<std::string>& sam_list,
            std::string& align_string_list,
            bool& is_complete_map)
{
    frag_align_list.clear();
    sam_list.clear();
    align_string_list.clear();
    string qname;
    while (load_next_align()) {
        if (qname.empty()) {
            qname.assign(m_columns[0].first, m_columns[0].second);
        } else {
            if (strncmp(qname.c_str(), m_columns[0].first, qname.size())) {
                unget_frag_align();
                break;
            }
        }
        m_frag_align.qas_offset = align_string_list.size();
        align_string_list += m_qas;
        m_frag_align.sas_offset = align_string_list.size();
        align_string_list += m_sas;
        sam_list.push_back(m_line);
        frag_align_list.push_back(m_frag_align);
    }
    if (frag_align_list.empty()) return 0;
    is_complete_map = (sam_list.front().find("mt:Z:complete_map") != string::npos);
    return frag_align_list.size();
}

////////////////

bool SamLoader::load_next_align_suffix()
{
    if (m_unget_frag_align) {
        m_unget_frag_align = false;
        return true;
    }

    if (HbnLineReaderAtEof(m_line_reader)) return false;
    HbnLineReaderReadOneLine(m_line_reader);
    m_line.assign(ks_s(m_line_reader->line), ks_size(m_line_reader->line));
    m_columns.clear();
    split_string_by_char(m_line.c_str(), m_line.size(), '\t', m_columns);
    extract_read_frag_id_from_read_name(m_columns[0].first, m_columns[0].second, &m_frag_align.read_id, &m_frag_align.frag_id);

    {
        const char* s = m_columns[0].first;
        int x = m_columns[0].second;
        while (x) {
            --x;
            if (s[x] == '_') break;
        }
        hbn_assert(x > 0);
        m_frag_align.qname.assign(s, x);
    }

    int flag = atoi(m_columns[1].first);
    m_frag_align.qdir = (flag & 0x10) ? REV : FWD;
    m_frag_align.sname.assign(m_columns[2].first, m_columns[2].second);
    if (m_columns[9].first[0] != '*') {
        if (m_frag_align.qdir == FWD) {
            m_fwd_read.assign(m_columns[9].first, m_columns[9].second);
            s_make_rc_read(m_fwd_read, m_rev_read);
        } else {
            m_rev_read.assign(m_columns[9].first, m_columns[9].second);
            s_make_rc_read(m_rev_read, m_fwd_read);
        }
    }
    const string& read = (m_frag_align.qdir == FWD) ? m_fwd_read : m_rev_read;
    int md_size = 0;
    const char* md = s_extract_md(m_line, md_size);
    s_parse_cigar(read.c_str(), read.size(), m_columns[5].first, m_columns[5].second, md, md_size,
        m_qas, m_sas, &m_frag_align.qoff, &m_frag_align.qend, &m_frag_align.soff, &m_frag_align.send, &m_frag_align.pi);
    m_frag_align.qsize = read.size();
    int soff = atoi(m_columns[3].first);
    hbn_assert(soff > 0);
    --soff;
    m_frag_align.soff += soff;
    m_frag_align.send += soff;
    m_frag_align.as_size = m_qas.size();

    if (m_frag_align.qdir == REV) {
        int x = m_frag_align.qsize - m_frag_align.qend;
        int y = m_frag_align.qsize - m_frag_align.qoff;
        m_frag_align.qoff = x;
        m_frag_align.qend = y;
    }
    m_frag_align.map_q = atoi(m_columns[4].first);

#if 0
    fprintf(stderr, "[%s, %d, %d, %d, %d] x [%s, %d, %d], %g\n", m_frag_align.qname.c_str(), m_frag_align.qdir,
        m_frag_align.qoff, m_frag_align.qend, m_frag_align.qsize,
        m_frag_align.sname.c_str(), m_frag_align.soff, m_frag_align.send, m_frag_align.pi);
#endif
    return true;
}

int SamLoader::load_next_align_list_suffix(std::vector<FragmentAlignment>& frag_align_list,
            std::vector<std::string>& sam_list,
            std::string& align_string_list,
            bool& is_complete_map)
{
    frag_align_list.clear();
    sam_list.clear();
    align_string_list.clear();
	int qid = -1;
    while (load_next_align_suffix()) {
        if (qid == -1) {
		qid = m_frag_align.read_id;
        } else {
            if (qid != m_frag_align.read_id) {
                unget_frag_align();
                break;
            }
        }
        m_frag_align.qas_offset = align_string_list.size();
        align_string_list += m_qas;
        m_frag_align.sas_offset = align_string_list.size();
        align_string_list += m_sas;
        sam_list.push_back(m_line);
        frag_align_list.push_back(m_frag_align);
    }
    if (frag_align_list.empty()) return 0;
    is_complete_map = (sam_list.front().find("mt:Z:complete_map") != string::npos);
    return frag_align_list.size();
}

//////////////////

void
read_frag_id_to_string(const int read_id, const int frag_id, char buf[])
{
    char buf1[64], buf2[64];
    u64_to_fixed_width_string_r(read_id, buf1, kReadIdSize);
    u64_to_fixed_width_string_r(frag_id, buf2, kFragIdSize);
    sprintf(buf, "%s:%s", buf1, buf2);
}

static void
s_illegal_frag_name_format(const int line, const char* s, const int sl)
{
    int read_id = 1, frag_id = 2;
    char buf[64];
    read_frag_id_to_string(read_id, frag_id, buf);
    fprintf(stderr, "[%d] Illegal frag name format '", line);
    hbn_fwrite(s, 1, sl, stderr);
    fprintf(stderr, "'\n");
    fprintf(stderr, "A legal frag name with read_id = %d and frag_id = %d should end with '%s'\n", read_id, frag_id, buf);
    abort();
}

void
extract_read_frag_id_from_read_name(const char* s, const int sl, int* read_id, int* frag_id)
{
    int p = sl;
    while (p) {
        --p;
        if (s[p] == '_') break;
    }
    if (s[p] != '_') s_illegal_frag_name_format(__LINE__, s, sl);
    ++p;
    if (p >= sl) s_illegal_frag_name_format(__LINE__, s, sl);
    if (!isdigit(s[p])) s_illegal_frag_name_format(__LINE__, s, sl);
    if (read_id) *read_id = atoi(s + p);

    while (p < sl && s[p] != ':') ++p;
    if (s[p] != ':') s_illegal_frag_name_format(__LINE__, s, sl);
    ++p;
    if (p >= sl) s_illegal_frag_name_format(__LINE__, s, sl);
    if (!isdigit(s[p])) s_illegal_frag_name_format(__LINE__, s, sl);
    if (frag_id) *frag_id = atoi(s + p);
}

void 
load_subject_id_from_sam(const char* sam_path, std::map<std::string, int>& name2id_map)
{
    HbnLineReader* in = HbnLineReaderNew(sam_path);
    vector<pair<const char*, int>> cols;
    string name;
    int id = 0;
    string line;
    while (!HbnLineReaderAtEof(in)) {
        HbnLineReaderReadOneLine(in);
        if (ks_front(in->line) != '@') break;
        line.assign(ks_s(in->line), ks_size(in->line));
        if (strncmp(line.c_str(), "@SQ", 3)) continue;
        cols.clear();
        split_string_by_char(line.c_str(), line.size(), '\t', cols);
        for (auto& c : cols) {
            if (strncmp(c.first, "SN:", 3)) continue;
            name.assign(c.first + 3, c.second - 3);
            auto pos = name2id_map.find(name);
            hbn_assert(pos == name2id_map.end());
            name2id_map.insert(pair<string, int>(name, id));
            ++id;
            break;
        }
    }
    in = HbnLineReaderFree(in);
    fprintf(stderr, "Load %d subjects\n", id);
}
