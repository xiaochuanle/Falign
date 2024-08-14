#ifndef __SEQ_NAME2ID_MAP_HPP
#define __SEQ_NAME2ID_MAP_HPP

#include "../ncbi_blast/str_util/ncbistr.hpp"

#include <map>
#include <string>
#include <vector>

class SeqName2IdMap
{
public:
    
    void clear() {
        m_name2id.clear();
        m_name_offset_list.clear();
        m_name_list.clear();
    }

    void add_one_name(const char* name, const int name_size) {
        m_name_offset_list.push_back(m_name_list.size());
        m_name_list.insert(m_name_list.end(), name, name + name_size);
        m_name_list.push_back('\0');
    }

    void build_name2id_map() {
        m_name2id.clear();
        int num_chr = NumSeqs();
        for (int i = 0; i < num_chr; ++i) {
            const char* name = m_name_list.data() + m_name_offset_list[i];
            const int name_size = strlen(name);
            NStr::CTempString cname(name, name_size);
            m_name2id.insert(std::pair<NStr::CTempString, int>(cname, i));
        }
    }

    void add_one_name(const char* name) {
        return add_one_name(name, strlen(name));
    }

    void add_one_name(const std::string& name) {
        return add_one_name(name.c_str(), name.size());
    }

    int GetIdFromNameSafe(const char* name, const int name_size) const {
        NStr::CTempString cname(name, name_size);
        auto pos = m_name2id.find(cname);
        if (pos == m_name2id.end()) {
            std::cerr << "ERROR: Sequence name " << cname << " does not exist" << std::endl;
            abort();
        }
        return pos->second;
    }

    int GetIdFromNameSafe(const char* name) const {
        return GetIdFromNameSafe(name, strlen(name));
    }

    int GetIdFromNameSafe(const std::string& name) const {
        return GetIdFromNameSafe(name.c_str(), name.size());
    }

    const char* GetSeqName(const int seq_id) const {
        const int n = NumSeqs();
        if (seq_id < 0 || seq_id >= n) {
            HBN_ERR("Sequence id %d is out of plausible range [%d, %d)", seq_id, 0, n);
        }
        return m_name_list.data() + m_name_offset_list[seq_id];
    }

    int NumSeqs() const {
        return m_name_offset_list.size();
    }

    bool NameExists(NStr::CTempString name) {
        return m_name2id.find(name) != m_name2id.end();
    }

private:
    std::map<NStr::CTempString, int>    m_name2id;
    std::vector<size_t>                 m_name_offset_list;
    std::vector<char>                   m_name_list;
};

#endif // __SEQ_NAME2ID_MAP_HPP