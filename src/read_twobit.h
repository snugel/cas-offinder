#pragma once
#include <string>
#include <vector>
#include <fstream>
#include "channel.h"

int read_twobit(std::string &filepath, std::vector<std::string> &chrnames, Channel<std::string> &content, std::vector<uint64_t> &chrpos);
