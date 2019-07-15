#pragma once

#include <filesystem/path.h>

#include <cstdint>
#include <fstream>
#include <vector>


namespace gpucf
{

template<typename R, class T>
std::vector<T> read(filesystem::path f)
{
    std::ifstream in(f.str(), std::ios::binary);
    
    uint64_t size;
    in.read(reinterpret_cast<char *>(&size), sizeof(uint64_t));

    std::vector<R> raw(size);
    in.read(reinterpret_cast<char *>(raw.data()), raw.size() * sizeof(R));

    std::vector<T> data;
    data.reserve(raw.size());

    for (const R &r : raw)
    {
        data.emplace_back(r);
    }

    return data;
}

template<typename R>
std::vector<R> read(filesystem::path f)
{
    std::ifstream in(f.str(), std::ios::binary);
    
    uint64_t size;
    in.read(reinterpret_cast<char *>(&size), sizeof(uint64_t));

    std::vector<R> raw(size);
    in.read(reinterpret_cast<char *>(raw.data()), raw.size() * sizeof(R));

    return raw;
}

} // namespace gpucf

// vim: set ts=4 sw=4 sts=4 expandtab:
