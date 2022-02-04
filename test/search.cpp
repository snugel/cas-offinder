#include "search.h"
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include <iostream>
#include "RangeIterator.h"
#include "bit4ops.h"
#include "ceildiv.h"
#include "test_framework.h"

void gold_search(uint8_t *bit4genome, uint64_t genome_size, uint8_t *bit4patterns,
                 uint64_t num_patterns, uint64_t pattern_size, uint32_t max_mismatches,
                 Match **match_result, uint64_t *num_matches) {
    uint64_t genome_blocks = cdiv(genome_size, 2);
    uint32_t pattern_blocks = cdiv(pattern_size, 2);
    uint32_t buf_size = 1 << 16;
    *match_result = (Match *)malloc(buf_size * sizeof(Match));
    *num_matches = 0;
    for (size_t i : range(genome_blocks - pattern_blocks + 1)) {
        for (size_t l : range(2)) {
            for (size_t j : range(num_patterns)) {
                uint32_t num_mismatches = pattern_size;
                for (size_t k : range(pattern_blocks)) {
                    uint16_t prev = bit4genome[i + k];
                    uint16_t next = bit4genome[i + k + 1];
                    uint8_t cur = ((prev >> (4 * l)) & 0xff) | ((next << (4 * (2 - l))) & 0xff);
                    num_mismatches -=
                        __builtin_popcount(cur & bit4patterns[j * pattern_blocks + k]);
                }
                if (num_mismatches <= max_mismatches) {
                    if (*num_matches == buf_size) {
                        buf_size *= 2;
                        *match_result = (Match *)realloc(*match_result, buf_size * sizeof(Match));
                    }
                    (*match_result)[*num_matches] = Match{.loc = uint64_t(i * 2 + l),
                                                          .pattern_idx = uint32_t(j),
                                                          .mismatches = uint32_t(num_mismatches)};
                    *num_matches += 1;
                }
            }
        }
    }
}

void sort_matches(Match *matches, size_t size) {
    std::sort(matches, matches + size, [](Match &a, Match &b) {
        if (a.loc < b.loc) return true;
        if (a.loc > b.loc) return false;
        if (a.pattern_idx < b.pattern_idx) return true;
        if (a.pattern_idx > b.pattern_idx) return false;
        return false;
    });
}
TEST(test_gold_search) {
    char genome[] = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    char patterns[][10] = {"GCGTAGACG", "CGTAGCTAG", "CAGTCGATC"};
    Match expected[] = {
        Match{
            .loc = 2,
            .pattern_idx = 0,
            .mismatches = 0,
        },
        Match{
            .loc = 5,
            .pattern_idx = 2,
            .mismatches = 2,
        },
        Match{
            .loc = 13,
            .pattern_idx = 2,
            .mismatches = 0,
        },
        Match{
            .loc = 21,
            .pattern_idx = 1,
            .mismatches = 0,
        },
    };
    Match *actual;
    uint64_t num_actual;
    constexpr int genome_size = sizeof(genome) - 1;
    constexpr int pattern_size = 9;
    constexpr int num_patterns = 3;
    uint8_t genome4bit[cdiv(genome_size, 2)] = {0};
    constexpr int pattern_blocks = cdiv(pattern_size, 2);
    uint8_t patterndata[num_patterns * pattern_blocks] = {0};
    str2bit4(genome4bit, genome, 0, genome_size);
    for (int i : range(num_patterns)) {
        uint8_t *cur_data = patterndata + i * pattern_blocks;
        str2bit4pattern(cur_data, patterns[i], 0, pattern_size);
    }
    int max_mismatches = 3;
    gold_search(genome4bit, genome_size, patterndata, num_patterns, pattern_size, max_mismatches,
                &actual, &num_actual);
    //    std::cerr << num_actual << "\n";
    //    for (size_t i : range(num_actual)) {
    //        std::cerr << actual[i].loc << "\t" << actual[i].pattern_idx << "\t" <<
    //        actual[i].mismatches
    //                  << "\n";
    //    }
    if (num_actual != 4) {
        return false;
    }
    bool result = !memcmp(expected, actual, sizeof(expected));
    free(actual);
    return result;
}

TEST(test_search_against_gold) {
    srand(42);
    constexpr size_t genome_size = 1 << 10;
    char genome[genome_size] = {0};
    char genome_chars[] = {'A', 'G', 'C', 'T'};
    for (size_t i : range(genome_size)) {
        genome[i] = genome_chars[rand() % sizeof(genome_chars)];
    }
    constexpr size_t genome_blocks = cdiv(genome_size, 2);
    uint8_t genome_data[genome_blocks] = {0};
    str2bit4(genome_data, genome, 0, genome_size);
    constexpr size_t num_patterns = 8;
    constexpr size_t pattern_size = 91;
    constexpr size_t pattern_blocks = cdiv(pattern_size, 2) + 1;
    uint8_t pattern_data[num_patterns * pattern_blocks] = {0};
    char pattern_chars[] = {'A', 'G', 'C', 'T', 'R', 'N', 'V'};
    for (size_t j : range(num_patterns)) {
        char pattern[pattern_size];
        for (size_t i : range(pattern_size)) {
            pattern[i] = pattern_chars[rand() % sizeof(pattern_chars)];
        }
        str2bit4pattern(&pattern_data[j * pattern_blocks], pattern, 0, pattern_size);
    }
    const int max_mismatches = 35;
    Match *gold_result;
    uint64_t num_gold_results;
    Match *actual_result;
    uint64_t num_actual_results;

    gold_search(genome_data, genome_size, pattern_data, num_patterns, pattern_size, max_mismatches,
                &gold_result, &num_gold_results);
    SearchFactory *fact = create_search_factory(CPU);
    int dev_idx = 0;
    Searcher *searcher =
        create_searcher(fact, dev_idx, 100000, 1000000, pattern_data, num_patterns, pattern_size);
    search(searcher, genome_data, genome_size, max_mismatches, &actual_result, &num_actual_results);
    free_searcher(&searcher);
    free_search_factory(&fact);
    sort_matches(actual_result, num_actual_results);
    //    std::cerr << num_gold_results << "\t" << num_actual_results << "\n";
    //    for (int i : range(10)) {
    //        std::cerr << gold_result[i].loc << "\t" << gold_result[i].pattern_idx << "\t"
    //                  << gold_result[i].mismatches << "\t";
    //        std::cerr << actual_result[i].loc << "\t" << actual_result[i].pattern_idx << "\t"
    //                  << actual_result[i].mismatches << "\n";
    //    }
    return num_gold_results == num_actual_results &&
           !memcmp(gold_result, actual_result, num_gold_results * sizeof(Match));
}
