#ifndef __PAF_READER_HPP
#define __PAF_READER_HPP

#include "../../corelib/line_reader.h"
#include "../../ncbi_blast/str_util/ncbistr.hpp"

#include <string>
#include <vector>

typedef struct {
    int qdir;
    int qoff;
    int qend;
    int qsize;
    int soff;
    int send;
    int ssize;
    int map_q;
    int sname_idx;
} PAF;

class PAFReader
{
public:
    PAFReader(const char* paf_path);
    ~PAFReader();
    int read_next_paf_list(std::vector<PAF>& paf_list, std::string& read_name, std::vector<std::string>& ref_name_list);
    int read_next_paf_list(std::vector<PAF>& paf_list, std::string& read_name, std::vector<std::string>& ref_name_list, std::vector<std::string>& line_list);

private:
    bool x_read_one_line();
    bool read_next_paf();
    void set_paf_is_loaded();

private:
    HbnLineReader*              m_line_reader;
    std::string                 m_line;
    std::vector<std::pair<const char*, int>>    m_cols;
    bool                        m_paf_is_loaded;
    PAF                         m_paf;
};

#endif // __PAF_READER_HPP