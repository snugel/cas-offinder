#pragma once
#include <stdint.h>
#include <stddef.h>
#include "search.h"

struct GenomeMatch
{
    const char* dna_match;
    const char* chrom_name;
    uint64_t chrom_loc;
    uint32_t pattern_idx;
    uint32_t mismatches;
};
typedef void async_match_callback(const GenomeMatch*);
void async_search(const char* genome_path,
                  DeviceType device_ty,
                  const char* compares,
                  size_t pattern_size,
                  size_t num_patterns,
                  uint32_t mismatches,
                  async_match_callback callback);

