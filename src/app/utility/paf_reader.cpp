#include "paf_reader.hpp"

#include "sam_reader.hpp"

using namespace std;

bool PAFReader::x_read_one_line()
{
    if (HbnLineReaderAtEof(m_line_reader)) return false;
    HbnLineReaderReadOneLine(m_line_reader);
    m_line.assign(ks_s(m_line_reader->line), ks_size(m_line_reader->line));
    m_cols.clear();
    split_string_by_char(m_line.c_str(), m_line.size(), '\t', m_cols);
    return true;
}

void PAFReader::set_paf_is_loaded()
{
    m_paf_is_loaded = true;
}

bool PAFReader::read_next_paf()
{
    if (m_paf_is_loaded) {
        m_paf_is_loaded = false;
        return true;
    }
    if (!x_read_one_line()) return false;

    m_paf.qdir = (m_cols[4].first[0] == '+') ? FWD : REV;
    m_paf.qsize = atoi(m_cols[1].first);
    m_paf.qoff = atoi(m_cols[2].first);
    m_paf.qend = atoi(m_cols[3].first);

    m_paf.ssize = atoi(m_cols[6].first);
    m_paf.soff = atoi(m_cols[7].first);
    m_paf.send = atoi(m_cols[8].first);

    m_paf.map_q = atoi(m_cols[11].first);

    return true;
}

PAFReader::PAFReader(const char* m4_path)
{
    m_line_reader = HbnLineReaderNew(m4_path);
    m_paf_is_loaded = false;
}

PAFReader::~PAFReader()
{
    HbnLineReaderFree(m_line_reader);
}

int PAFReader::read_next_paf_list(std::vector<PAF>& paf_list, std::string& read_name, std::vector<std::string>& ref_name_list)
{
    paf_list.clear();
    read_name.clear();
    ref_name_list.clear();

    if (!read_next_paf()) return 0;
    m_paf.sname_idx = ref_name_list.size();
    paf_list.push_back(m_paf);
    read_name.assign(m_cols[0].first, m_cols[0].second);
    ref_name_list.push_back(string(m_cols[5].first, m_cols[5].second));

    while (read_next_paf()) {
        if (read_name.size() != m_cols[0].second || strncmp(read_name.c_str(), m_cols[0].first, m_cols[0].second)) {
            set_paf_is_loaded();
            break;
        }
        m_paf.sname_idx = ref_name_list.size();
        paf_list.push_back(m_paf);
        ref_name_list.push_back(string(m_cols[5].first, m_cols[5].second));
    }
    return paf_list.size();
}

int PAFReader::read_next_paf_list(std::vector<PAF>& paf_list, std::string& read_name, std::vector<std::string>& ref_name_list, std::vector<std::string>& line_list)
{
    paf_list.clear();
    read_name.clear();
    ref_name_list.clear();
    line_list.clear();

    if (!read_next_paf()) return 0;
    m_paf.sname_idx = ref_name_list.size();
    paf_list.push_back(m_paf);
    read_name.assign(m_cols[0].first, m_cols[0].second);
    ref_name_list.push_back(string(m_cols[5].first, m_cols[5].second));
    line_list.push_back(m_line);

    while (read_next_paf()) {
        if (read_name.size() != m_cols[0].second || strncmp(read_name.c_str(), m_cols[0].first, m_cols[0].second)) {
            set_paf_is_loaded();
            break;
        }
        m_paf.sname_idx = ref_name_list.size();
        paf_list.push_back(m_paf);
        ref_name_list.push_back(string(m_cols[5].first, m_cols[5].second));
        line_list.push_back(m_line);
    }
    return paf_list.size();
}