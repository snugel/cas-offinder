#include "RangeIterator.h"
#include "ceildiv.h"
#include "search.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

struct SearchFactory
{};
struct Searcher
{
    uint64_t max_matches;
    uint64_t max_genome_size;
    uint8_t* bit4patterns;
    uint64_t num_patterns;
    uint64_t pattern_size;
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
    return new Searcher{
        .max_matches = max_matches,
        .max_genome_size = max_genome_size,
        .bit4patterns = bit4patterns,
        .num_patterns = num_patterns,
        .pattern_size = pattern_size,
    };
}
void search(Searcher* searcher,
            uint8_t* bit4genome,
            uint64_t genome_size,
            uint32_t max_mismatches,
            Match** match_result,
            uint64_t* num_matches)
{

    uint8_t* bit4patterns = searcher->bit4patterns;
    uint64_t num_patterns = searcher->num_patterns;
    uint64_t pattern_size = searcher->pattern_size;

    uint64_t genome_blocks = cdiv(genome_size, 2);
    uint32_t pattern_blocks = cdiv(pattern_size, 2);
    uint32_t buf_size = 1 << 16;
    *match_result = (Match*)malloc(buf_size * sizeof(Match));
    *num_matches = 0;
    for (size_t i : range(genome_blocks - pattern_blocks + 1)) {
        for (size_t l : range(2)) {
            for (size_t j : range(num_patterns)) {
                uint32_t num_mismatches = pattern_size;
                for (size_t k : range(pattern_blocks)) {
                    uint16_t prev = bit4genome[i + k];
                    uint16_t next = bit4genome[i + k + 1];
                    uint8_t cur = ((prev >> (4 * l)) & 0xff) |
                                  ((next << (4 * (2 - l))) & 0xff);
                    num_mismatches -= __builtin_popcount(
                      cur & bit4patterns[j * pattern_blocks + k]);
                }
                if (num_mismatches <= max_mismatches) {
                    if (*num_matches == buf_size) {
                        buf_size *= 2;
                        *match_result = (Match*)realloc(
                          *match_result, buf_size * sizeof(Match));
                    }
                    (*match_result)[*num_matches] =
                      Match{ .loc = uint64_t(i * 2 + l),
                             .pattern_idx = uint32_t(j),
                             .mismatches = uint32_t(num_mismatches) };
                    *num_matches += 1;
                }
            }
        }
    }
}
void free_searcher(Searcher**) {}
void free_search_factory(SearchFactory**) {}
