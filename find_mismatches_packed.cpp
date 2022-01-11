#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>

#define LOCAL_BLOCK_SIZE 4

std::vector<uint64_t> make4bitpackedints(std::string genome)
{
    std::vector<uint64_t> result((genome.size() + 15) / 16);
    for (size_t i = 0; i < genome.size(); i++)
    {
        int ridx = i / 16;
        int shift = 15 - i % 16;
        result[ridx] |= uint64_t(to4bit(genome[i])) << (shift * 4);
    }
    return result;
}

void find_matches_packed_helper(
    uint64_t *genome,
    size_t genome_idx,
    uint64_t *pattern_blocks,
    size_t num_patterns,
    size_t blocks_per_pattern,
    size_t pattern_block_idx,
    std::vector<match> &matches,
    int max_mismatches,
    int max_matches
){
    // genome is expected to be at least BIGGER than the genome_size
    constexpr size_t blocks_avail = 16;
    uint16_t pattern_counts[LOCAL_BLOCK_SIZE] = {0};
    for (size_t k = 0; k < blocks_avail; k++)
    {
        size_t local_pattern_size = std::min(num_patterns - pattern_block_idx, size_t(LOCAL_BLOCK_SIZE));
        for (size_t x = 0; x < local_pattern_size; x++)
        {
            pattern_counts[x] = 0;
        }
        for (size_t l = 0; l < blocks_per_pattern; l++)
        {
            uint64_t prev = genome[genome_idx + l];
            uint64_t next = genome[genome_idx + l + 1];
            uint64_t cur = k == 0 ? prev : (prev << (k * 4)) | (next >> ((blocks_avail - k) * 4));
            for (size_t x = 0; x < local_pattern_size; x++)
            {
                pattern_counts[x] += __builtin_popcountl(cur & pattern_blocks[(pattern_block_idx + x)*blocks_per_pattern + l]);
            }
        }
        for (size_t x = 0; x < local_pattern_size; x++)
        {
            int mismatches = max_matches - pattern_counts[x];
            if (mismatches <= max_mismatches)
            {
                matches.push_back(match{
                    .loc = genome_idx * blocks_avail + k,
                    .mismatches = mismatches,
                    .pattern_idx = pattern_block_idx + x,
                });
            }
        }
    }
}

std::vector<match> find_matches_packed(std::string genome, std::vector<std::string> patterns, int max_mismatches)
{
    if (patterns.size() == 0)
    {
        return {{}};
    }
    size_t pattern_size = patterns[0].size();
    for (std::string &p : patterns)
    {
        assert(p.size() == pattern_size);
    }
    size_t blocks_per_pattern = (pattern_size + 15) / 16;

    std::vector<uint64_t> pattern_blocks(patterns.size() * blocks_per_pattern);
    for (size_t i = 0; i < patterns.size(); i++)
    {
        std::vector<uint64_t> b4pattern = make4bitpackedints(patterns[i]);
        assert(b4pattern.size() == blocks_per_pattern);
        for (size_t j = 0; j < blocks_per_pattern; j++)
        {
            pattern_blocks[i * blocks_per_pattern + j] = b4pattern[j];
        }
    }
    std::vector<match> matches;

    std::vector<uint64_t> b4genome = make4bitpackedints(genome);
    b4genome.push_back(0);

    size_t genome_size = b4genome.size() - blocks_per_pattern;

    for (size_t genome_idx = 0; genome_idx < genome_size; genome_idx += 1) {
        for (size_t pattern_block_idx = 0; pattern_block_idx < patterns.size(); pattern_block_idx += LOCAL_BLOCK_SIZE) {
            find_matches_packed_helper(
                b4genome.data(),                      // genome
                genome_idx,                           // genome_idx
                pattern_blocks.data(),                // pattern_blocks
                patterns.size(),                      // num_patterns
                blocks_per_pattern,                   // blocks_per_pattern
                pattern_block_idx,                    // pattern_block_idx
                matches,                              // matches
                max_mismatches,                       // max_mismatches
                pattern_size                          // max_matches
            );
        }
    }

    return matches;
}

TEST(test_make4bitpackedints)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint64_t> expected = {0x4212184142148241, 0x8214821841284182, 0x8148100000000000};
    std::vector<uint64_t> actual = make4bitpackedints(genome);
    return expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin());
}

TEST(test_find_matches_packed)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACGGCGTAGACG",
        "CGTAGCTAGCGTAGCTAG",
        "GATCGACTGGATCGACTG",
    };
    int mismatches = 9;
    std::vector<match> expected = find_matches_gold(genome, patterns, mismatches);
    std::vector<match> actual = find_matches_packed(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    for (match m : actual)
    {
        atomic_print_match(m);
    }
    std::cout << "expected\n";
    for (match m : expected)
    {
        atomic_print_match(m);
    }
    return matches_equal(actual, expected);
}

TEST(find_mismatches_packed_perf)
{
    std::vector<std::string> patterns(25, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for (int i : range(1000000))
    {
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 5;
    int start = clock();
    std::vector<match> actual = find_matches_packed(genome, patterns, mismatches);
    int end = clock();
    std::cout << "time: " << (end - start) / double(CLOCKS_PER_SEC) << "\n";
    return true;
}
