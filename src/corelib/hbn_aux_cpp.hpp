#ifndef __HBN_AUX_CPP_HPP
#define __HBN_AUX_CPP_HPP

#include "hbn_aux.h"

#include <fstream>
#include <iostream>
#include <sstream>

template <class stream_type>
bool open_stream(const char* path, stream_type& stream)
{
    if (stream.is_open()) stream.close();
    stream.open(path);
    if (!stream) {
        HBN_ERR("Fail to open file '%s': %s", path, strerror(errno));
    }
    return true;
}

template <class stream_type, typename mode_type>
bool open_stream(const char* path, stream_type& stream, mode_type mode)
{
    if (stream.is_open()) stream.close();
    stream.open(path, mode);
    if (!stream) {
        HBN_ERR("Fail to open file '%s': %s", path, strerror(errno));
    }
    return true;
}

#define open_ifstream(stream__, stream_path__) open_stream<std::ifstream>(stream__, stream_path__)
#define open_ofstream(stream__, stream_path__) open_stream<std::ofstream>(stream__, stream_path__)

#define d_open_stream(stream_type__, stream__, stream_path__) stream_type__ stream__; open_stream(stream_path__, stream__)

#define dm_open_stream(stream_type__, stream__, stream_path__, mode__) stream_type__ stream__; open_stream(stream_path__, stream__, mode__)

#define d_ifstream(stream__, stream_path__) d_open_stream(std::ifstream, stream__, stream_path__)

#define d_ofstream(stream__, stream_path__) d_open_stream(std::ofstream, stream__, stream_path__)

#endif // __HBN_AUX_CPP_HPP