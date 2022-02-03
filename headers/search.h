#pragma once
#include <inttypes.h>
struct Searcher;
struct Match {
    uint64_t loc;
    uint32_t pattern_idx;
    uint32_t mismatches;
};
Searcher* create_searcher(uint8_t* bit4patterns, uint64_t num_patterns, uint64_t pattern_size);
void search(Searcher* search, uint8_t* bit4genome, uint64_t genome_size, uint32_t max_mismatches,
            Match** match_result, uint64_t* num_matches);
void free_searcher(Searcher**);
