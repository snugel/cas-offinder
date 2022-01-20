#include "postprocess.h"
#include "RangeIterator.h"
#include "test/test_framework.h"
#include <vector>
#include <iostream>
#include <algorithm>

void get_chrom_info(
    const std::vector<match> & matches, 
    const std::vector<uint64_t> & chrom_poses,
    std::vector<uint64_t> & rel_locs,
    std::vector<uint64_t> & chrom_idxs
){
    rel_locs.resize(matches.size());
    chrom_idxs.resize(matches.size());
    for(size_t i : range(matches.size())){
        auto less_eq_el = std::upper_bound(chrom_poses.begin(), chrom_poses.end(), matches[i].loc) - 1;
        size_t pos_idx = less_eq_el - chrom_poses.begin();
        chrom_idxs[i] = pos_idx;
        rel_locs[i] = matches[i].loc - *less_eq_el;
    }
}

void filter_chrom_walls(
    std::vector<match> & matches, 
    const std::vector<uint64_t> & chrom_poses,
    size_t query_size
){
    matches.erase(std::remove_if(matches.begin(), matches.end(), [&](match m){
        auto greater_than_el = std::upper_bound(chrom_poses.begin(), chrom_poses.end(), m.loc);
        return  (*greater_than_el < m.loc + query_size);
    }), matches.end());
}
