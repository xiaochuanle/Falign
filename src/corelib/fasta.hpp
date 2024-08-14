#ifndef __FASTA_HPP
#define __FASTA_HPP

#include "line_reader.hpp"

#include <string>

class CFastaReader
{
public:
    CFastaReader(const char* path, const bool ignore_ill_formated_seq = true):
        m_LineReader(path),
        m_IgnoreIllFormatedSequence(ignore_ill_formated_seq) {}

    ~CFastaReader() {}

    int ReadOneSeq();

    std::string& sequence() {
        return m_Sequence;
    }

    std::string& name() {
        return m_Name;
    }

    std::string& comment() {
        return m_Comment;
    }

    std::string& plus() {
        return m_Plus;
    }

    std::string& qual() {
        return m_Qual;
    }

    const std::string& sequence() const {
        return m_Sequence;
    }

    const std::string& name() const {
        return m_Name;
    }

    const std::string& comment() const {
        return m_Comment;
    }

    const std::string& plus() const {
        return m_Plus;
    }

    const std::string& qual() const {
        return m_Qual;
    }

    void ClearSeqInfo() {
        m_Name.clear();
        m_Comment.clear();
        m_Sequence.clear();
        m_Plus.clear();
        m_Qual.clear();
    }

private:
    bool ParseDefLine(CTempString line);

    bool CheckDataLine(CTempString line);

    bool ParseDataLine(CTempString line);

    void AdvanceToNextSeqHeader();

    bool ParseFastqSeq();

private:
    HbnLineReader   m_LineReader;
    std::string     m_Name;
    std::string     m_Comment;
    std::string     m_Sequence;
    std::string     m_Plus;
    std::string     m_Qual;
    bool            m_IgnoreIllFormatedSequence;
};

typedef CFastaReader HbnFastaReader;

#endif // __FASTA_HPP