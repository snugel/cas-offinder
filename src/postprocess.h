#include <vector>
#include "find_mismatches.h"

struct chrom_info{
    uint64_t rel_loc;
    uint64_t chrom_idx;
};

chrom_info get_chrom_info(
    uint64_t loc, 
    const std::vector<uint64_t> & chrom_poses
);

bool crosses_chrom_wall(
    uint64_t loc, 
    const std::vector<uint64_t> & chrom_poses,
    size_t query_size
);
