#include <vector>
#include "find_mismatches.h"

void get_chrom_info(
    const std::vector<match> & matches, 
    const std::vector<uint64_t> & chrom_poses,
    std::vector<uint64_t> & rel_locs,
    std::vector<uint64_t> & chrom_idxs
);

void filter_chrom_walls(
    std::vector<match> & matches, 
    const std::vector<uint64_t> & chrom_poses,
    size_t query_size
);