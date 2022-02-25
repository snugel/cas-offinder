#pragma once
#include <inttypes.h>

void str2bit4(uint8_t* dest,
              const char* src,
              int64_t write_offset,
              uint64_t n_chrs);
void str2bit4pattern(uint8_t* dest,
                     const char* src,
                     int64_t write_offset,
                     uint64_t n_chrs);
void bit42str(char* dest,
              const uint8_t* src,
              int64_t read_offset,
              uint64_t n_chrs);
void twobit2bit4(uint8_t* dest, const uint8_t* src, uint64_t n_chrs);
void memsetbit4(uint8_t* dest, uint8_t bit4val, uint64_t start, uint64_t end);
bool is_mixedbase(const char* src, uint64_t n_chrs);
bool is_match(char nucl, char pattern);

inline uint8_t to_upper(uint8_t c)
{
    return ((c)&0xdf);
}
