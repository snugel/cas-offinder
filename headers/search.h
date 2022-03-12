#pragma once
#include <inttypes.h>
struct SearchFactory;
struct Searcher;
struct Match
{
    uint64_t loc;
    uint32_t pattern_idx;
    uint32_t mismatches;
};
enum DeviceType
{
    CPU,
    GPU,
    ACCEL,
    ANY_DEV
};
SearchFactory* create_search_factory(DeviceType device);
int num_searchers_avaliable(SearchFactory* fact);
void print_device_information(SearchFactory* fact);
Searcher* create_searcher(SearchFactory* factory,
                          int searcher_idx,
                          uint64_t max_genome_size,
                          uint64_t max_matches,
                          uint8_t* bit4patterns,
                          uint64_t num_patterns,
                          uint64_t pattern_size);
void search(Searcher* search,
            uint8_t* bit4genome,
            uint64_t genome_size,
            uint32_t max_mismatches,
            Match** match_result,
            uint64_t* num_matches);

void free_searcher(Searcher**);
void free_search_factory(SearchFactory**);
