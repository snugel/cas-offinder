#pragma once
#include <inttypes.h>

void str2bit4(uint8_t * dest, const char * src, int64_t write_offset, uint64_t n_chrs);
void str2bit4pattern(uint8_t * dest, const char * src, int64_t write_offset, uint64_t n_chrs);
void bit42str(char * dest, const uint8_t * src, int64_t read_offset, uint64_t n_chrs);

