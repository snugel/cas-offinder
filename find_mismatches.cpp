#include "find_mismatches.h"
#include "test/test_framework.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>




uint64_t to4bit(char c){
    constexpr uint64_t G = 0x1;
    constexpr uint64_t C = 0x2;
    constexpr uint64_t A = 0x4;
    constexpr uint64_t T = 0x8;
    switch (c)
    {
    case 'G':
        return G;
    case 'C':
        return C;
    case 'A':
        return A;
    case 'T':
        return T;
    case 'R':
        return A | G;
    case 'Y':
        return C | T;
    case 'S':
        return G | C;
    case 'W':
        return A | T;
    case 'K':
        return G | T;
    case 'M':
        return A | C;
    case 'B':
        return C | G | T;
    case 'D':
        return A | G | T;
    case 'H':
        return A | C | T;
    case 'V':
        return A | C | G;
    case 'N':
        return A | C | G | T;
    default:
        throw std::runtime_error("got unexpected letter, input must be mixed base");
    }
}

std::vector<uint64_t> to4bit(std::string genome){
    std::vector<uint64_t> result((genome.size() + 15) / 16);
    for(size_t i = 0; i < genome.size(); i++){
        int ridx = i / 16;
        int shift = 15 - i % 16;
        result[ridx] |= to4bit(genome[i]) << (shift*4);
    }
    return result;
}

void find_mismatches(std::vector<uint64_t> & genome, std::vector<uint64_t> & pattern_blocks, std::vector<uint16_t*> & b4_counts){
    constexpr int blocks_avail = (128-64)/4;
    for(size_t i = 0; i < genome.size(); i++){
        uint64_t prev = genome[i];
        uint64_t next = i+1 < genome.size() ? genome[i+1] : 0;
        for(size_t k = 0; k < blocks_avail; k++){
            uint64_t cur = (prev << (k*4)) | (next >> ((blocks_avail - k)*4));
            for(size_t j = 0; j < pattern_blocks.size(); j++){
                b4_counts[j][i*16+k] -= __builtin_popcountll(cur & pattern_blocks[j]);
            }
        }
    }
}


std::vector<std::vector<uint16_t>> find_mismatches_packed(std::string genome, std::vector<std::string> patterns){
    if(patterns.size() == 0){
        return {{}};
    }
    size_t pattern_size = patterns[0].size();
    for(std::string & p : patterns){
        assert(p.size() == pattern_size);
    }
    size_t blocks_per_pattern = (pattern_size + 15) / 16; 
    size_t num_pattern_blocks = patterns.size() * blocks_per_pattern;
    
    std::vector<std::vector<uint16_t>> padded_pattern_mismatches(patterns.size(), std::vector<uint16_t>(genome.size() + (blocks_per_pattern+1)*16, uint16_t(pattern_size)));

    std::vector<uint64_t> pattern_blocks(num_pattern_blocks);
    std::vector<uint16_t*> block_counts(num_pattern_blocks);
    for(size_t i = 0; i < patterns.size(); i++){
        std::vector<uint64_t> b4pattern = to4bit(patterns[i]);
        assert(b4pattern.size() == blocks_per_pattern);
        for(size_t j = 0; j < blocks_per_pattern; j++){
            block_counts[i*blocks_per_pattern + j] = &padded_pattern_mismatches[i][(blocks_per_pattern-j-1)*16];
            pattern_blocks[i*blocks_per_pattern + j] = b4pattern[j];
        }
    }
    std::vector<uint64_t> b4genome = to4bit(genome);
    int start = clock();
    find_mismatches(b4genome, pattern_blocks, block_counts);
    int end = clock();
    std::cout << "innertime: " << (end - start)/double(CLOCKS_PER_SEC) << "\n";

    //removes memory pad on patterns
    for(auto & v : padded_pattern_mismatches){
        v.erase(v.begin(), v.begin() + (blocks_per_pattern-1)*16);
        v.resize(genome.size()-pattern_size+1);
    }
    // returns unpadded patterns
    return padded_pattern_mismatches;
}

struct match_info{
    int x;
};
std::vector<std::vector<uint16_t>> find_mismatches(std::string genome,  std::vector<std::string> patterns){
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
                    if(mismatches[j][i] > 5){
                        break;
                    }
                }
            }
        }
    }
    return mismatches;
}

std::vector<std::vector<uint16_t>> find_mismatches_gold(std::string genome, std::vector<std::string> patterns){
    if(patterns.size() == 0){
        return {{}};
    }
    size_t pattern_size = patterns[0].size();
    for(std::string & p : patterns){
        assert(p.size() == pattern_size);
    }
    std::vector<std::vector<uint16_t>> mismatches(patterns.size(), std::vector<uint16_t>(genome.size() - pattern_size + 1, uint16_t(pattern_size)));
    for(int j = 0; j < patterns.size(); j++){
        for(int i = 0; i < genome.size() - pattern_size + 1; i++){
            for(int k = 0; k < pattern_size; k++){
                if(to4bit(genome.at(i+k)) & to4bit(patterns.at(j).at(k))){
                    mismatches.at(j).at(i) -= 1;
                }
            }
        }
    }
    return mismatches;
}

TEST(find_mismatches_gold){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACG",
        "CGTAGCTAG"
    };
    std::vector<std::vector<uint16_t>> expected = {
        //6, 9, 0, 9, 7, 6, 8, 6, 6, 8, 4, 9, 8, 5, 7, 7, 6, 6, 7, 8, 4, 8, 7, 8, 5, 6, 9, 7, 6,
        { 6, 9, 0, 9, 7, 6, 8, 6, 6, 8, 4, 9, 8, 5, 7, 7, 7, 6, 7, 8, 4, 8, 7, 8, 5, 6, 9, 7, 6, },
        //6, 6, 8, 4, 7, 8, 6, 7, 7, 7, 6, 7, 7, 8, 4, 9, 5, 4, 7, 9, 9, 0, 9, 9, 7, 5, 6, 8, 6,
        { 8, 6, 8, 4, 7, 8, 6, 7, 7, 7, 6, 7, 7, 8, 4, 9, 7, 4, 7, 9, 9, 0, 9, 9, 7, 5, 6, 8, 6, },
    };
    std::vector<std::vector<uint16_t>> actual = find_mismatches(genome, patterns);
    for(auto & v : actual){
        for(uint16_t i : v){
            std::cout << i << ", ";
        }
        std::cout << "\n";
    }
    if(actual.size() != expected.size() || actual[0].size() != expected[0].size()){
        return false;
    }
    for(int i = 0; i < actual.size(); i++){
        for(int j = 0; j < actual[0].size(); j++){
            if(actual.at(i).at(j) != expected.at(i).at(j)){
                return false;
            }
        }
    }
    return true;
}
TEST(find_mismatches_perf){
    std::vector<std::string> patterns(25, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for(int i = 0; i < 100000; i++){
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int start = clock();
    find_mismatches(genome, patterns);
    int end = clock();
    std::cout << "time: " << (end - start)/double(CLOCKS_PER_SEC) << "\n";
}


TEST(test_to_4bit_vec){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint64_t> expected = {0x4212184142148241, 0x8214821841284182, 0x8148100000000000};
    std::vector<uint64_t> actual = to4bit(genome);
    return expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin());
}