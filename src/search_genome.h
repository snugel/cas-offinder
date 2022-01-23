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
    std::vector<chromloc> chrom_locs, 
    int mismatches, 
    Channel<GenomeInput> * inp_channel, 
    Channel<GenomeMatch> * out_channel
);
