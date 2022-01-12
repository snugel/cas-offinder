#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include "opencl_executor.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>

#define LOCAL_BLOCK_SIZE 4

std::vector<uint64_t> make4bitpackedints(std::string genome);


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
    size_t blocks_per_pattern = (pattern_size + 15) / 16;

    std::vector<uint64_t> pattern_blocks(patterns.size() * blocks_per_pattern);
    for (size_t i = 0; i < patterns.size(); i++)
    {
        std::vector<uint64_t> b4pattern = make4bitpackedints(patterns[i]);
        assert(b4pattern.size() == blocks_per_pattern);
        for (size_t j = 0; j < blocks_per_pattern; j++)
        {
            pattern_blocks[i * blocks_per_pattern + j] = b4pattern[j];
        }
    }
    std::vector<match> matches;

    std::vector<uint64_t> b4genome = make4bitpackedints(genome);

    OpenCLPlatform plat;
    OpenCLExecutor executor(
        "find_mismatches_packed.cl", 
        plat.get_platform(), 
        plat.get_first_device()
    );
    const int OUT_BUF_SIZE = 1<<20;
    CLBuffer<match> output_buf = executor.new_clbuffer<match>(OUT_BUF_SIZE);
    CLBuffer<int> output_count = executor.new_clbuffer<int>(1);
    CLBuffer<uint64_t> genome_buf = executor.new_clbuffer<uint64_t>(b4genome.size()+2);
    CLBuffer<uint64_t> pattern_buf = executor.new_clbuffer<uint64_t>(pattern_blocks.size());

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

    genome_buf.clear_buffer();
    pattern_buf.clear_buffer();
    output_count.clear_buffer();
    output_buf.clear_buffer();
    genome_buf.write_buffer(&b4genome[0], b4genome.size());
    pattern_buf.write_buffer(pattern_blocks);
    compute_results.run();
    int out_count = 0;
    output_count.read_buffer(&out_count, 1);
    matches.resize(out_count);
    output_buf.read_buffer(&matches[0], out_count);
    std::cout << "out count" << out_count << "\n";
    return matches;
    
    const size_t CHUNK_SIZE = 256;
    for (size_t genome_block = 0; genome_block < genome_size; genome_block += CHUNK_SIZE) {
        match match_buffer[1024];
        int match_idx = 0;
        for (size_t genome_idx = 0; genome_idx < std::min(genome_size-genome_block,CHUNK_SIZE); genome_idx += 1) {
            for (size_t pattern_block_idx = 0; pattern_block_idx < patterns.size(); pattern_block_idx += LOCAL_BLOCK_SIZE) {
                // find_matches_packed_helper(
                //     b4genome.data() + genome_block,       // genome
                //     genome_idx,                           // genome_idx
                //     pattern_blocks.data(),                // pattern_blocks
                //     patterns.size(),                      // num_patterns
                //     blocks_per_pattern,                   // blocks_per_pattern
                //     pattern_block_idx,                    // pattern_block_idx
                //     max_mismatches,                       // max_mismatches
                //     pattern_size,                         // max_matches
                //     match_buffer,
                //     match_idx
                // );
            }
        }
        matches.insert(matches.end(), match_buffer, match_buffer + match_idx);
    }
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