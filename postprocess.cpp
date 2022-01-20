#include "postprocess.h"
#include "RangeIterator.h"
#include "test/test_framework.h"
#include <vector>
#include <iostream>
#include <algorithm>


chrom_info get_chrom_info(
    uint64_t loc, 
    const std::vector<uint64_t> & chrom_poses
){
    auto less_eq_el = std::upper_bound(chrom_poses.begin(), chrom_poses.end(), loc) - 1;
    size_t pos_idx = less_eq_el - chrom_poses.begin();
    return chrom_info{
        .rel_loc=loc - *less_eq_el,
        .chrom_idx=pos_idx,
    };
}

bool crosses_chrom_wall(
    uint64_t loc, 
    const std::vector<uint64_t> & chrom_poses,
    size_t query_size
){
    auto greater_than_el = std::upper_bound(chrom_poses.begin(), chrom_poses.end(), loc);
    return  (*greater_than_el < loc + query_size);
}

