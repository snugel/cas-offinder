#include "postprocess.h"
#include "test/test_framework.h"
#include <iostream>


TEST(test_get_chrom_info){
    std::vector<match> expected = {
        match{ .loc=2,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=5,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=10,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=13,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=18,  .mismatches=0,  .pattern_idx=0,  },
    };
    std::vector<uint64_t> chrom_posses = {
        0,
        5,
        7,
        11,
        13,
        22
    };
    std::vector<uint64_t> rel_locs;
    std::vector<uint64_t> chrom_idxs;
    get_chrom_info(expected, chrom_posses, rel_locs, chrom_idxs);
    
    std::vector<uint64_t> expected_rel_locs = { 2, 0, 3, 0, 5 };
    std::vector<uint64_t> expected_chrom_idxs = { 0, 1, 2, 4, 4 };
    return rel_locs.size() == expected_rel_locs.size() && std::equal(rel_locs.begin(), rel_locs.end(), expected_rel_locs.begin()) && 
          chrom_idxs.size() == expected_chrom_idxs.size() && std::equal(chrom_idxs.begin(), chrom_idxs.end(), expected_chrom_idxs.begin());
}


TEST(test_filter_chrom_walls){
    std::vector<match> input = {
        match{ .loc=2,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=5,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=10,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=13,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=18,  .mismatches=0,  .pattern_idx=0,  },
    };
    std::vector<uint64_t> chrom_posses = {
        0,
        5,
        7,
        11,
        13,
        22
    };
    int query_size = 2;
    filter_chrom_walls(input, chrom_posses, query_size);
    
    std::vector<match> expected = {
        match{ .loc=2,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=5,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=13,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=18,  .mismatches=0,  .pattern_idx=0,  },
    };
    // for(match m : input){
    //     atomic_print_match(m);
    // }
    return matches_equal(input, expected);
}