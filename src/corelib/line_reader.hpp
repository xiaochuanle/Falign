#ifndef __LINE_READER_HPP
#define __LINE_READER_HPP

#include "hbn_aux.h"
#include "../ncbi_blast/str_util/ncbistr.hpp"
#include "../htslib/hts.h"

#include <string>

using NStr::CTempString;

class CBufferLineReader
{
public:
    /*
    Read from the file, "-" (but not "./") means standard input
    An explicit call to operator++ or
    ReadLine() will bbe necessaray to fetch the first line
    */
    CBufferLineReader(const char* filename) {
        M_FileName = filename;
        M_in = hts_open(filename, "r");
        hts_set_threads(M_in, 8);
        M_ungetline = false;
        ks_initialize(&M_kline);
        M_LineNumber = 0;
    }

    ~CBufferLineReader() {
        hts_close(M_in);
        ks_free(&M_kline);
    }

    void UngetLine() {
        M_ungetline = true;
        --M_LineNumber;
    }
    
    CTempString operator*() const {
        return M_Line;
    }

    size_t GetLineNumber() const {
        return M_LineNumber;
    }

    const char* GetFileName() const { return M_FileName.c_str(); }

    bool ReadOneLine() {
        if (M_ungetline) {
            M_ungetline = false;
            ++M_LineNumber;
            return true;
        }
        int r = hts_getline(M_in, '\n', &M_kline);
        if (r == -1) return false;
        if (r <= -2) HBN_ERR("FAIL at reading line from file %s", M_FileName.c_str());
        M_Line.assign(ks_c_str(&M_kline), ks_len(&M_kline));
        ++M_LineNumber;
        return true;
    }

private:
    CBufferLineReader(const CBufferLineReader&);
    CBufferLineReader& operator=(const CBufferLineReader&);
    
private:
    std::string     M_FileName;
    htsFile*        M_in;
    bool            M_ungetline;
    kstring_t       M_kline;
    CTempString     M_Line;
    size_t          M_LineNumber;
};

typedef CBufferLineReader HbnLineReader;

#endif // __LINE_READER_HPP