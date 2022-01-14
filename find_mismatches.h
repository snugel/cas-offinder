#pragma once
#include<vector>
#include<string>


struct match{
    uint64_t loc;
    uint32_t mismatches;
    uint32_t pattern_idx;
};
std::vector<match> find_matches(std::string genome, std::vector<std::string> patterns, int max_mismatches);
void sort_matches(std::vector<match> & matches);
void atomic_print_match(match & m);
bool matches_equal(std::vector<match> & m1, std::vector<match> & m2);


// internal facing, use externally at your own risk
std::vector<match> find_matches_gold(std::string genome, std::vector<std::string> patterns, int max_mismatches);
