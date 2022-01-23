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

std::vector<match> find_matches_gold(std::string & genome,std::vector<std::string> & patterns, int max_mismatches){
    std::vector<match> matches;
    size_t pattern_size = patterns.at(0).size();
    for(size_t gi : range(genome.size() - pattern_size + 1)){
        for(size_t pi : range(patterns.size())){
            uint32_t mismatches = 0;
            for(size_t i : range(pattern_size)){
                mismatches += !(to4bit(genome[gi+i]) & to4bit(patterns[pi][i]));
                if(mismatches > max_mismatches){
                    break;
                }
            }
            if(mismatches <= max_mismatches){
                matches.push_back(match{
                    .loc=gi,
                    .mismatches=mismatches,
                    .pattern_idx=uint32_t(pi),
                });
            }
        }
    }
    return matches;
}