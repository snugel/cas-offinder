#pragma once
#include<vector>
#include<functional>
#include<string>
#include "channel.h"


struct GenomeInput{
    std::shared_ptr<uint32_t> data; // block of genome memory 
    size_t size;    // size of genome to search. Must be pattern_size+1 smaller than the actual data 
    uint64_t idx; // start index of genome chunk
};
struct match{
    uint64_t loc;
    uint32_t mismatches;
    uint32_t pattern_idx;
};
struct WorkerOutput{
    GenomeInput input;
    std::shared_ptr<match> matches;
    size_t num_matches;
};
void sort_matches(std::vector<match> & matches);
void atomic_print_match(match & m);
bool matches_equal(std::vector<match> & m1, std::vector<match> & m2);


std::vector<match> find_matches(std::string & genome, std::vector<std::string> & patterns, int max_mismatches);
std::vector<match> find_matches_gold(std::string & genome, std::vector<std::string> & patterns, int max_mismatches);

void find_matches_worker(
    Channel<GenomeInput> * genome_data_stream, 
    std::vector<std::string> patterns, 
    int max_mismatches, 
    Channel<WorkerOutput> * out_stream
);
