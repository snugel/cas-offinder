#pragma once

#include <array>
#include <vector>

extern std::array<char, 256> to4bitmap;

inline char to4bit(char c){
    return to4bitmap[c];
}

std::vector<uint64_t> make4bitpackedint64(const std::string & genome);
std::vector<uint32_t> make4bitpackedint32(const std::string & genome);
void clean_bogus(std::string & genome);

std::vector<uint64_t> bit64tobit32(const std::vector<uint32_t> & bit32);
  
