#pragma once
#include <memory>

struct FileChunk{
    std::shared_ptr<char[]> data;
    size_t size;
};
