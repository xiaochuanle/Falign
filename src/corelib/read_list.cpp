#include "read_list.hpp"
#include "db_format.h"
#include "hbn_aux.h"
#include "line_reader.hpp"

#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>

using namespace std;

static BOOL 
s_extract_read_list_from_dir(const char* dir_path, std::vector<std::string>& read_list)
{
    struct stat statbuf;
    if (stat(dir_path, &statbuf)) {
        HBN_ERR("Fail to stat path '%s': %s", dir_path, strerror(errno));
    }
    if (!S_ISDIR(statbuf.st_mode)) return FALSE;

    DIR* dp = NULL;
    struct dirent* entry = NULL;
    if ((dp = opendir(dir_path)) == NULL) {
        HBN_ERR("Fail to open directory '%s': %s", dir_path, strerror(errno));
    }
    struct stat file_stat;
    char path[HBN_MAX_PATH_LEN];
    while ((entry = readdir(dp)) != NULL) {
        sprintf(path, "%s/%s", dir_path, entry->d_name);
        if (stat(path, &file_stat)) HBN_ERR("Fail to stat path '%s': %s", path, strerror(errno));
        if (S_ISDIR(file_stat.st_mode)) continue;
        read_list.push_back(path);
    }
    closedir(dp);
    return TRUE;
}

static void
s_extract_read_list_from_list_file(const char* file_list_path, std::vector<std::string>& read_list)
{
    HbnLineReader* line_reader = new HbnLineReader(file_list_path);
    string sline;
    while (line_reader->ReadOneLine()) {
        NStr::CTempString line = **line_reader;
        sline.assign(line.data(), line.size());
        if (access(sline.c_str(), F_OK)) {
            HBN_LOG("Sequence file '%s' does not exist\n", sline.c_str());
            abort();
        }
        read_list.push_back(sline);
    }
    delete line_reader;
}

void make_read_list(int argc, char* argv[], std::vector<std::string>& read_list)
{
    for (int i = 0; i < argc; ++i) {
        if (s_extract_read_list_from_dir(argv[i], read_list)) continue;

        if (access(argv[i], F_OK)) {
            HBN_LOG("Sequence file '%s' does not exist.\n", argv[i]);
            abort();
        }

        EDbFormat format = hbn_guess_db_format(argv[i]);
        if (format == eDbFormatEmptyFile) {
            HBN_LOG("File '%s' contains no sequence data, skip it.\n", argv[i]);
            continue;
        }

        if (format == eDbFormatFasta || format == eDbFormatFastq) {
            read_list.push_back(argv[i]);
        } else if (format == eDbFormatUnknown) {
            s_extract_read_list_from_list_file(argv[i], read_list);
        }
    }
}