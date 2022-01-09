#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>

std::vector<uint64_t> make4bitpackedints(std::string genome){
    std::vector<uint64_t> result((genome.size() + 15) / 16);
    for(size_t i = 0; i < genome.size(); i++){
        int ridx = i / 16;
        int shift = 15 - i % 16;
        result[ridx] |= to4bit(genome[i]) << (shift*4);
    }
    return result;
}

void find_mismatches_packed_helper(
        uint64_t * genome, 
        int genome_size, 
        uint64_t * pattern_blocks, 
        int num_patterns, 
        int blocks_per_pattern,
        std::vector<match> & matches, 
        int max_mismatches,
        int max_matches
){
    // genome is expected to be at least blocks_per_pattern+1 BIGGER than the genome_size
    constexpr size_t blocks_avail = 16;
    constexpr size_t pattern_count_arr_size = 4;
    uint16_t pattern_counts[pattern_count_arr_size];
    for(size_t i = 0; i < genome_size; i++){
        for(size_t j = 0; j < num_patterns; j += pattern_count_arr_size){
            size_t local_pattern_size = std::min(num_patterns-j, pattern_count_arr_size);
            for(size_t x = 0; x < local_pattern_size; x++){
                pattern_counts[x] = 0;
            }
            for(size_t l = 0; l < blocks_per_pattern; l++){
                for(size_t k = 0; k < blocks_avail; k++){
                    uint64_t prev = genome[i+l];
                    uint64_t next = genome[i+l+1];
                    uint64_t cur = k == 0 ? prev : (prev << (k*4)) | (next >> ((blocks_avail - k)*4));
                    for(size_t x = 0; x < local_pattern_size; x++){
                        pattern_counts[x] += __builtin_popcountl(cur & pattern_blocks[j+x]);
                    }
                }
            }
            for(size_t x = 0; x < local_pattern_size; x++){
                int mismatches = max_matches - pattern_counts[x];
                if(mismatches >= max_mismatches){
                    matches.push_back(match{
                        .loc=i,
                        .mismatches=mismatches,
                        .pattern_idx=j+x
                    });
                }
            }
        }
    }
}

std::vector<std::vector<uint16_t>> find_mismatches_packed(std::string genome, std::vector<std::string> patterns, int max_mismatches){
    if(patterns.size() == 0){
        return {{}};
    }
    size_t pattern_size = patterns[0].size();
    for(std::string & p : patterns){
        assert(p.size() == pattern_size);
    }
    size_t blocks_per_pattern = (pattern_size + 15) / 16; 
    size_t num_pattern_blocks = patterns.size() * blocks_per_pattern;
    
    std::vector<uint64_t> pattern_blocks(num_pattern_blocks);
    std::vector<uint16_t*> block_counts(num_pattern_blocks);
    for(size_t i = 0; i < patterns.size(); i++){
        std::vector<uint64_t> b4pattern = make4bitpackedints(patterns[i]);
        assert(b4pattern.size() == blocks_per_pattern);
        for(size_t j = 0; j < blocks_per_pattern; j++){
            block_counts[i*blocks_per_pattern + j] = &padded_pattern_mismatches[i][(blocks_per_pattern-j-1)*16];
            pattern_blocks[i*blocks_per_pattern + j] = b4pattern[j];
        }
    }
    std::vector<uint64_t> b4genome = make4bitpackedints(genome);
    find_mismatches_packed_helper(b4genome, pattern_blocks, block_counts, max_mismatches);

    //removes memory pad on patterns
    for(auto & v : padded_pattern_mismatches){
        v.erase(v.begin(), v.begin() + (blocks_per_pattern-1)*16);
        v.resize(genome.size()-pattern_size+1);
    }
    // returns unpadded patterns
    return padded_pattern_mismatches;
}

std::vector<std::vector<uint16_t>> find_mismatches_packed(std::string genome,  std::vector<std::string> patterns, int max_mismatches){
    if(patterns.size() == 0){
        return {{}};
    }
    size_t pattern_size = patterns[0].size();
    for(std::string & p : patterns){
        assert(p.size() == pattern_size);
    }
    
    std::vector<std::vector<uint16_t>> mismatches(patterns.size(), std::vector<uint16_t>(genome.size() - pattern_size + 1));
    for(char & c : genome){
        c = to4bit(c);
    }
    for(std::string & pat : patterns){
        for(char & c : pat){
            c = to4bit(c);
        }
    }
    for(int j = 0; j < patterns.size(); j++){
        for(int i = 0; i < genome.size() - pattern_size + 1; i++){
            for(int k = 0; k < pattern_size; k++){
                if(!(genome[i+k] & patterns[j][k])){
                    mismatches[j][i] += 1;
                }
            }
        }
    }
    return mismatches;
}

TEST(test_to_4bit_vec){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint64_t> expected = {0x4212184142148241, 0x8214821841284182, 0x8148100000000000};
    std::vector<uint64_t> actual = make4bitpackedints(genome);
    return expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin());
}

TEST(test_find_matches_gold){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACG",
        "CGTAGCTAG",
        "GATCGACTG"
    };
    int mismatches = 3;
    std::vector<match> expected = find_matches_gold(genome, patterns, mismatches);
    std::vector<match> actual = find_mismatches_packed(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    for(match m : actual){
        atomic_print_match(m);
    }
    return matches_equal(actual, expected);
}
TEST(find_mismatches_perf){
    std::vector<std::string> patterns(25, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for(int i : range(10000)){
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 5;
    int start = clock();
    std::vector<match> actual = find_matches_packed(genome, patterns, mismatches);
    int end = clock();
    std::cout << "time: " << (end - start)/double(CLOCKS_PER_SEC) << "\n";
    return true;
}
