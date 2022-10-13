#ifndef __SAM_READER_HPP
#define __SAM_READER_HPP

#include "../../corelib/hbn_aux.h"
#include "../../corelib/line_reader.h"

#include <map>
#include <string>
#include <vector>

typedef struct {
    std::string qname;
    int read_id;
    int frag_id;
    int hap_type;
    int qdir;
    int qoff;
    int qend;
    int qsize;
    std::string sname;
    int soff;
    int send;
    double pi;
    int map_q;
    size_t qas_offset;
    size_t sas_offset;
    int as_size;
} FragmentAlignment;

class SamLoader
{
public:
    SamLoader(const char* sam_path);
    ~SamLoader();
    int load_next_align_list(std::vector<FragmentAlignment>& frag_align_list,
            std::vector<std::string>& sam_list,
            std::string& align_string_list,
            bool& is_complete_map);
    int load_next_align_list_suffix(std::vector<FragmentAlignment>& frag_align_list,
            std::vector<std::string>& sam_list,
            std::string& align_string_list,
            bool& is_complete_map);
    bool load_next_align_suffix();

    FragmentAlignment& align() {
        return m_frag_align;
    }
    std::string& sam_line() {
        return m_line;
    }
    std::vector<std::string>& sam_hdr_list() {
        return m_sam_hdr_list;
    }

private:
    bool load_next_align();
    void unget_frag_align();

private:
    HbnLineReader*      m_line_reader;
    std::string         m_line;
    std::vector<std::pair<const char*, int>>
                        m_columns;
    std::string         m_fwd_read;
    std::string         m_rev_read;
    FragmentAlignment   m_frag_align;
    bool                m_unget_frag_align;
    std::string         m_qas;
    std::string         m_sas;
    std::vector<std::string>
                        m_sam_hdr_list;
};

void
split_string_by_char(const char* s, const int sl, const char delim, std::vector<std::pair<const char*, int>>& tokens);

void
read_frag_id_to_string(const int read_id, const int frag_id, char buf[]);

void
extract_read_frag_id_from_read_name(const char* s, const int sl, int* read_id, int* frag_id);

void 
load_subject_id_from_sam(const char* sam_path, std::map<std::string, int>& name2id_map);

#endif // __SAM_READER_HPP