#ifndef __ARG_PARSE_HPP
#define __ARG_PARSE_HPP

#include <cstring>

#include "../ncbi_blast/str_util/ncbistr.hpp"

static inline bool parse_bool_arg_value(int argc, char* argv[], int& i, const char* arg_name, bool& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    s = true;
    i += 1;

    return true;
}

static inline bool parse_int_arg_value(int argc, char* argv[], int& i, const char* arg_name, int& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = NStr::StringToInt8(argv[i+1]);
    i += 2;

    return true;
}

static inline bool parse_real_arg_value(int argc, char* argv[], int& i, const char* arg_name, double& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = NStr::StringToDouble(argv[i+1]);
    i += 2;

    return true;
}

static inline bool parse_data_size_arg_value(int argc, char* argv[], int& i, const char* arg_name, size_t& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = NStr::StringToUInt8_DataSize(argv[i+1]);
    i += 2;

    return true;
}

static inline bool parse_string_arg_value(int argc, char* argv[], int& i, const char* arg_name, const char*& s)
{
    if (strcmp(argv[i], arg_name)) return false;

    if (i == argc) {
        fprintf(stderr, "ERROR: value to argument '%s' is missing\n", argv[i]);
        exit(1);
    }
    s = argv[i+1];
    i += 2;

    return true;
}

#endif // __ARG_PARSE_HPP