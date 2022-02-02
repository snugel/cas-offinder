#pragma once
#include <string>
#include <vector>
#include <fstream>
#include "channel.h"
#include "read_file.h"

int read_twobit(std::string &filepath, std::vector<std::string> &chrnames, Channel<FileChunk> &content, std::vector<uint64_t> &chrpos);
