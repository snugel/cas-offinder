#include "test/test_framework.h"
#include "find_mismatches.h"
#include <iostream>
#include "RangeIterator.h"
#include "timing.h"

TEST(test_find_matches_gold){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACG",
        "CGTAGCTAG",
        "CAGTCGATC"
    };
    int mismatches = 3;
    std::vector<match> expected = {
        match{ .loc=2,  .mismatches=0,  .pattern_idx=0,  },
        match{ .loc=21,  .mismatches=0,  .pattern_idx=1,  },
        match{ .loc=5,  .mismatches=2,  .pattern_idx=2,  },
        match{ .loc=13,  .mismatches=0,  .pattern_idx=2,  },
    };
    std::vector<match> actual = find_matches_gold(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    for(match m : actual){
        atomic_print_match(m);
    }
    return matches_equal(actual, expected);
}

TEST(find_mismatches_packed_perf)
{
    std::vector<std::string> patterns(25, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for (int i : range(10000))
    {
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 5;
    std::cout << "time: " << time_spent([&](){
    find_matches_gold(genome, patterns, mismatches);
    }) << std::endl;
    return true;
}




TEST(test_find_matches_opencl_small)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATGACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACGGCGTAGACG",
        "CGTAGCTAGCGTAGCTAG",
        "GATCGACTGGATCGACTG",
    };
    int mismatches = 12;
    std::vector<match> expected = find_matches_gold(genome, patterns, mismatches);
    std::vector<match> actual = find_matches(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    for (match m : actual)
    {
        atomic_print_match(m);
    }
    std::cout << "expected\n";
    for (match m : expected)
    {
        atomic_print_match(m);
    }
    return matches_equal(actual, expected);
}


TEST(test_find_matches_opencl_large)
{
    std::string genome;
    for (int i : range(1000000))
    {
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    std::vector<std::string> patterns = {
        "GCGTAGACGGCGTAGACG",
        "CGTAGCTAGCGTAGCTAG",
        "GATCGACTGGATCGACTG",
    };
    int mismatches = 7;
    std::vector<match> expected = find_matches_gold(genome, patterns, mismatches);
    std::vector<match> actual = find_matches(genome, patterns, mismatches);
    sort_matches(actual);
    sort_matches(expected);
    std::cout << "large output size: " << std::dec << actual.size() << "\n";
    return matches_equal(actual, expected);
}

TEST(find_mismatches_opencl_perf)
{
    std::vector<std::string> patterns(50, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for (int i : range(10000))
    {
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 10;
    std::cout << "time: " << time_spent([&](){
    find_matches(genome, patterns, mismatches);
    }) << std::endl;
    return true;
}
