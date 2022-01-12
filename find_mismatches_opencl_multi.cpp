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



std::vector<match> find_matches_opencl_multi(std::string genome, std::vector<std::string> patterns, int max_mismatches)
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
    size_t next_genome_idx = 0;
    const size_t CHUNK_SIZE = 1<<16;
    #pragma omp parallel for
    for(size_t device_idx = 0; device_idx < plat.num_cl_devices(); device_idx++){
        OpenCLExecutor executor(
            "find_mismatches_packed.cl", 
            plat.get_platform(), 
            plat.get_first_device()
        );

        const int OUT_BUF_SIZE = 1<<24;
        CLBuffer<match> output_buf = executor.new_clbuffer<match>(OUT_BUF_SIZE);
        CLBuffer<int> output_count = executor.new_clbuffer<int>(1);
        CLBuffer<uint64_t> genome_buf = executor.new_clbuffer<uint64_t>(CHUNK_SIZE+3);
        CLBuffer<uint64_t> pattern_buf = executor.new_clbuffer<uint64_t>(pattern_blocks.size());

        size_t genome_size = CHUNK_SIZE - blocks_per_pattern + 1;
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
            }
        );
        pattern_buf.clear_buffer();
        pattern_buf.write_buffer(pattern_blocks);
        output_buf.clear_buffer();
        while(true){
            size_t genome_idx;
            #pragma omp critical 
            {
                genome_idx = next_genome_idx;
                next_genome_idx += CHUNK_SIZE;
            }
            if(genome_idx >= b4genome.size()){
                break;
            }
        
            genome_buf.clear_buffer();
            output_count.clear_buffer();
            genome_buf.write_buffer(&b4genome[genome_idx], std::min(CHUNK_SIZE, b4genome.size() - genome_idx));
            compute_results.run();
            int out_count = 0;
            output_count.read_buffer(&out_count, 1);

            if(out_count > 0){
                #pragma omp critical 
                {
                    size_t old_size = matches.size();
                    matches.resize(old_size + out_count);
                    output_buf.read_buffer(&matches[old_size], out_count);
                }
            }
        }
    }
    return matches;
}



TEST(test_find_matches_opencl_multi)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATGACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATGTCTGATG";
    std::vector<std::string> patterns = {
        "GCGTAGACGGCGTAGACG",
        "CGTAGCTAGCGTAGCTAG",
        "GATCGACTGGATCGACTG",
    };
    int mismatches = 12;
    std::vector<match> expected = find_matches_gold(genome, patterns, mismatches);
    std::vector<match> actual = find_matches_opencl_multi(genome, patterns, mismatches);
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

TEST(find_mismatches_opencl_multi_perf)
{
    std::vector<std::string> patterns(50, "GCGTAGACGGCGTAGACGGCGTANNRGR");
    std::string genome;
    for (int i : range(100000000))
    {
        genome += "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    }
    int mismatches = 10;
    int start = clock();
    std::vector<match> actual = find_matches_opencl_multi(genome, patterns, mismatches);
    int end = clock();
    std::cout << "time: " << (end - start) / double(CLOCKS_PER_SEC) << "\n";
    return true;
}
