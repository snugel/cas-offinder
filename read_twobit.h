#pragma once
#include <string>
#include <vector>
#include <fstream>

using namespace std;

int read_twobit(string &filepath, vector<string> &chrnames, string &content, vector<uint64_t> &chrpos);
