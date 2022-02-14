#pragma once
#include <string>
#include <vector>
#include "channel.h"
#include "read_file.h"

int read_folder(std::string &filepath, std::vector<std::string> &chrnames, Channel<FileChunk> &content, std::vector<uint64_t> &chrpos);

