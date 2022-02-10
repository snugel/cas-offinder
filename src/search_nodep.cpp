#include "RangeIterator.h"
#include "ceildiv.h"
#include "search.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <memory>
#include <cstdlib>

using namespace std;

struct SearchFactory
{};
struct Searcher
{
    uint64_t max_matches;
    uint64_t max_genome_size;
    uint64_t num_patterns;
    uint64_t pattern_size;
    Match * out_buf;
    uint64_t * bit4blocks;
};

SearchFactory* create_search_factory(DeviceType device_ty)
{
    if (device_ty != CPU) {
        cerr << "Nodep search only supports CPU!\n";
        exit(1);
    }
    return new SearchFactory();
}
int num_searchers_avaliable(SearchFactory*)
{
    return 1;
}
Searcher* create_searcher(SearchFactory*,
                          int,
                          uint64_t max_matches,
                          uint64_t max_genome_size,
                          uint8_t* bit4patterns,
                          uint64_t num_patterns,
                          uint64_t pattern_size)
{

    int old_pattern_size = cdiv(pattern_size, 2);
    int new_pattern_blocks = cdiv(pattern_size, 2*sizeof(uint64_t));
    int new_pattern_size = new_pattern_blocks* sizeof(uint64_t);
    uint64_t * data = (uint64_t*)calloc(new_pattern_blocks*num_patterns,sizeof(uint64_t));
    uint8_t * padded_patterns = (uint8_t*)(data);
    for (size_t j : range(num_patterns)) {
        for (size_t i : range(old_pattern_size)) {
            padded_patterns[new_pattern_size * j + i] =
              bit4patterns[old_pattern_size * j + i];
        }
    }
    return new Searcher{
        .max_matches = max_matches,
        .max_genome_size = max_genome_size,
        .num_patterns = num_patterns,
        .pattern_size = pattern_size,
        .out_buf = (Match*)malloc(sizeof(Match)*max_matches),
        .bit4blocks= data,
    };
}
void search(Searcher* searcher,
            uint8_t* bit4genome,
            uint64_t genome_size,
            uint32_t max_mismatches,
            Match** match_result,
            uint64_t* num_matches)
{

    uint64_t* bit4blocks = searcher->bit4blocks;
    uint64_t num_patterns = searcher->num_patterns;
    uint64_t pattern_size = searcher->pattern_size;
    constexpr size_t block_size = 2*sizeof(uint64_t);
    uint64_t* genome_block_data = (uint64_t*)(bit4genome);

    uint64_t genome_blocks = cdiv(genome_size, block_size);
    uint32_t pattern_blocks = cdiv(pattern_size, block_size);
    uint64_t out_count = 0;
    for (size_t i : range(genome_blocks - pattern_blocks + 1)) {
        for (size_t l : range(block_size)) {
            for (size_t j : range(num_patterns)) {
                uint32_t num_mismatches = pattern_size;
                for (size_t k : range(pattern_blocks)) {
                    uint64_t prev = genome_block_data[i + k];
                    uint64_t next = genome_block_data[i + k + 1];
                    uint64_t cur = l == 0 ? prev : ((prev >> (4 * l))) |
                                  ((next << (4 * (block_size - l))));
                    num_mismatches -= __builtin_popcountll(
                      cur & bit4blocks[j * pattern_blocks + k]);
                }
                if (num_mismatches <= max_mismatches) {
                    searcher->out_buf[out_count] =
                      Match{ .loc = uint64_t(i * block_size + l),
                             .pattern_idx = uint32_t(j),
                             .mismatches = uint32_t(num_mismatches) };
                    out_count += 1;
                }
            }
        }
    }
    *match_result = (Match*)malloc(out_count * sizeof(Match));
    memcpy(*match_result, searcher->out_buf, out_count * sizeof(Match));
    *num_matches = out_count;
}
void free_searcher(Searcher** sptr) {
    free((*sptr)->out_buf);
    free((*sptr)->bit4blocks);
    delete *sptr;
    *sptr = nullptr;
}
void free_search_factory(SearchFactory** fact) {
    delete *fact;
    *fact = nullptr;
}
