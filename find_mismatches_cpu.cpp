#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include "bit4ops.h"
#include "timing.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <sstream>
#include <ctime>

void find_matches_packed_helper(
    uint64_t *genome,
    size_t genome_idx,
    uint64_t *pattern_blocks,
    size_t num_patterns,
    size_t blocks_per_pattern,
    size_t pattern_idx,
    int max_mismatches,
    int max_matches,
    match * match_buffer,
    int & match_idx
){
    // genome is expected to be at least BIGGER than the genome_size
    constexpr size_t blocks_avail = 16;
    for (size_t k = 0; k < blocks_avail; k++)
    {
        uint64_t pattern_count = 0;
        for (size_t l = 0; l < blocks_per_pattern; l++)
        {
            uint64_t prev = genome[genome_idx + l];
            uint64_t next = genome[genome_idx + l + 1];
            uint64_t cur = k == 0 ? prev : (prev << (k * 4)) | (next >> ((blocks_avail - k) * 4));
            pattern_count += __builtin_popcountl(cur & pattern_blocks[pattern_idx*blocks_per_pattern + l]);
        }
       // std::cout <<int(pattern_counts[0]) << std::endl;
        uint32_t mismatches = max_matches - pattern_count;
        if (mismatches <= max_mismatches)
        {
            int next_idx = match_idx++;
            match_buffer[next_idx] = match{
                .loc = genome_idx * blocks_avail + k,
                .mismatches = mismatches,
                .pattern_idx = uint32_t(pattern_idx),
            };
        }
    }
}

std::vector<match> find_matches_gold(std::string & genome,std::vector<std::string> & patterns, int max_mismatches){
    std::vector<uint32_t> genomeb4 = make4bitpackedint32(genome);
    std::vector<match> matches;
    find_matches_gold(genomeb4, patterns, max_mismatches, [&](match m){
        matches.push_back(m);
    });
    return matches;
}
void find_matches_gold(std::vector<uint32_t> & genomeb4, std::vector<std::string> & patterns, int max_mismatches, std::function<void(match)> func)
{
    if (patterns.size() == 0)
    {
        return;
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
        std::vector<uint64_t> b4pattern = make4bitpackedint64(patterns[i]);
        assert(b4pattern.size() == blocks_per_pattern);
        for (size_t j = 0; j < blocks_per_pattern; j++)
        {
            pattern_blocks[i * blocks_per_pattern + j] = b4pattern[j];
        }
    }
    std::vector<match> matches;
    std::vector<uint64_t> b4genome = bit64tobit32(genomeb4);// ((genomeb4.size()+1)/2 + 1);
    //std::copy(genomeb4.begin(), genomeb4.end(), (uint32_t*)&b4genome[0]);
    b4genome.push_back(0);    

    size_t genome_size = b4genome.size() - blocks_per_pattern;

    #pragma omp parallel for
    for (size_t genome_idx = 0; genome_idx < genome_size; genome_idx += 1) {
        match match_buffer[1<<8];
        int match_idx = 0;
        for (size_t pattern_idx = 0; pattern_idx < patterns.size(); pattern_idx += 1) {
            find_matches_packed_helper(
                    b4genome.data(),       // genome
                    genome_idx,                           // genome_idx
                    pattern_blocks.data(),                // pattern_blocks
                    patterns.size(),                      // num_patterns
                    blocks_per_pattern,                   // blocks_per_pattern
                    pattern_idx,                    // pattern_idx
                    max_mismatches,                       // max_mismatches
                    pattern_size,                         // max_matches
                    match_buffer,
                    match_idx
                    );
        }
        
        #pragma omp critical 
        {
            for(size_t i : range(match_idx)){
                func(match_buffer[i]);
            }
        }
    }

}
