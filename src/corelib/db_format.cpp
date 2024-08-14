#include "db_format.h"

#include "line_reader.hpp"

#include <string>

using namespace std;

extern "C"
EDbFormat
hbn_guess_db_format(const char* path)
{
    EDbFormat fmt = eDbFormatEmptyFile;
    HbnLineReader* reader = new HbnLineReader(path);
    while (reader->ReadOneLine()) {
        NStr::CTempString line = **reader;
        NStr::TruncateSpacesInPlace(line);
        // ignore lines containing only whitespace
        if (line.empty()) continue;
        char c = line[0];
        // no content, just a comment or blank line
        if (c == '!' || c == '#' || c == ';') continue;
        if (c == '>') {
            fmt = eDbFormatFasta;
        } else if (c == '@') {
            fmt = eDbFormatFastq;
        } else {
            fmt = eDbFormatUnknown;
        }
        // we just test the first valid line (not empty line, not comment line)
        break;
    }
    delete reader;
    return fmt;
}