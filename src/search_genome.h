#pragma once
#include <string>
#include <memory>
#include "chromloc.h"
#include "find_mismatches.h"
#include "channel.h"

struct GenomeMatch{
    std::string dna_match;
    std::string cromosome;
    uint64_t chrom_loc;
    uint32_t pattern_idx;
    uint32_t mismatches;
};
void search_genome(
    std::vector<std::string> patterns, 
    std::string path,
    int mismatches, 
    Channel<GenomeMatch> * out_channel
);
