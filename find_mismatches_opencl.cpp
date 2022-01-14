#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include "opencl_executor.h"
#include "bit4ops.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>

#define LOCAL_BLOCK_SIZE 1
using block_ty = uint32_t;
constexpr size_t bit4_c = sizeof(block_ty) * 8 / 4;

std::vector<match> find_matches_opencl(std::string genome, std::vector<std::string> patterns, int max_mismatches)
{
    if (patterns.size() == 0)
    {
        return {{}};
    }
    size_t num_patterns = patterns.size();
    int pattern_size = patterns[0].size();
    for (std::string &p : patterns)
    {
        assert(p.size() == pattern_size);
    }
    size_t blocks_per_pattern = (pattern_size + bit4_c - 1) / bit4_c;

    std::vector<block_ty> pattern_blocks(patterns.size() * blocks_per_pattern);
    for (size_t i = 0; i < patterns.size(); i++)
    {
        std::vector<block_ty> b4pattern = make4bitpackedint32(patterns[i]);
        assert(b4pattern.size() == blocks_per_pattern);
        for (size_t j = 0; j < blocks_per_pattern; j++)
        {
            pattern_blocks[i * blocks_per_pattern + j] = b4pattern[j];
        }
    }
    std::vector<match> matches;


    OpenCLPlatform plat;
    OpenCLExecutor executor(
        "find_mismatches_packed.cl", 
        plat.get_platform(), 
        plat.get_first_device()
    );
    std::vector<block_ty> b4genome = make4bitpackedint32(genome);

    const int OUT_BUF_SIZE = 1<<24;
    CLBuffer<match> output_buf = executor.new_clbuffer<match>(OUT_BUF_SIZE);
    CLBuffer<int> output_count = executor.new_clbuffer<int>(1);
    CLBuffer<block_ty> genome_buf = executor.new_clbuffer<block_ty>(b4genome.size()+2);
    CLBuffer<block_ty> pattern_buf = executor.new_clbuffer<block_ty>(pattern_blocks.size());

    size_t genome_size = b4genome.size() - blocks_per_pattern + 1;
    size_t num_pattern_blocks = (patterns.size() + LOCAL_BLOCK_SIZE - 1) / LOCAL_BLOCK_SIZE;
    CLKernel compute_results = executor.new_clkernel(
        "find_matches_packed_helper",
        CL_NDRange(genome_size, num_pattern_blocks),
        CL_NDRange(),
        CL_NDRange(),
        {
            genome_buf.k_arg(),
            pattern_buf.k_arg(),
            make_arg(num_patterns),
            make_arg(blocks_per_pattern),
            make_arg(max_mismatches),
            make_arg(pattern_size),
            output_buf.k_arg(),
            output_count.k_arg(),
        });
        int start = clock();
    genome_buf.clear_buffer();
    pattern_buf.clear_buffer();
    output_count.clear_buffer();
    output_buf.clear_buffer();
    genome_buf.write_buffer(&b4genome[0], b4genome.size());
    pattern_buf.write_buffer(pattern_blocks);
    compute_results.run();
    int out_count = 0;
    output_count.read_buffer(&out_count, 1);

    if(out_count > 0){
        matches.resize(out_count);
        output_buf.read_buffer(&matches[0], out_count);
    }
    int end = clock();
    std::cout << "inner time: " << (end - start) / double(CLOCKS_PER_SEC) << "\n";
    std::cout << "out count" << out_count << "\n";
    return matches;
}



TEST(test_find_matches_opencl)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATGACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACGGCGTAGACG",
        "CGTAGCTAGCGTAGCTAG",
        "GATCGACTGGATCGACTG",
    };
    int mismatches = 12;
    std::vector<match> expected = find_matches_gold(genome, patterns, mismatches);
    std::vector<match> actual = find_matches_opencl(genome, patterns, mismatches);
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

TEST(find_mismatches_opencl_perf)
{
    std::vector<std::string> patterns(50, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for (int i : range(100000000))
    {
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 10;
    int start = clock();
    std::vector<match> actual = find_matches_opencl(genome, patterns, mismatches);
    int end = clock();
    std::cout << "time: " << (end - start) / double(CLOCKS_PER_SEC) << "\n";
    return true;
}
