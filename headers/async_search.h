#pragma once
#include "search.h"
#include <stddef.h>
#include <stdint.h>

struct GenomeMatch
{
    const char* dna_match;
    const char* chrom_name;
    uint64_t chrom_loc;
    uint32_t pattern_idx;
    uint32_t mismatches;
};
struct BlockConfig{
    size_t OUT_CHUNK_SIZE;
    size_t MIN_CMPS_PER_OUT;
    size_t DEFAULT_CHUNK_BYTES;
};

typedef void async_match_callback(const GenomeMatch*);
void async_search(const char* genome_path,
                  DeviceType device_ty,
                  const char* compares,
                  size_t pattern_size,
                  size_t num_patterns,
                  uint32_t mismatches,
                  const BlockConfig * config,
                  async_match_callback callback);
