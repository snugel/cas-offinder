#pragma once
#include <string>
#include <vector>
#include "read_file.h"
#include "channel.h"

int read_fasta(std::string &filepath, std::vector<std::string> &chrnames, Channel<FileChunk> &content, std::vector<uint64_t> &chrpos);

