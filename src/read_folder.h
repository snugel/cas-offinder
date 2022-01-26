#pragma once
#include <string>
#include <vector>
#include "channel.h"

int read_folder(std::string &filepath, std::vector<std::string> &chrnames, Channel<std::string> &content, std::vector<uint64_t> &chrpos);

