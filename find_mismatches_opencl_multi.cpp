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

    OpenCLPlatform plat;
    size_t next_genome_idx = 0;
    constexpr size_t GENOME_CHUNK_SIZE = 1<<24-1;
    constexpr size_t B4_CHUNK_SIZE = (GENOME_CHUNK_SIZE + 15) / 16 + 1;
    #pragma omp parallel for
    for(size_t device_idx = 0; device_idx < plat.num_cl_devices(); device_idx++){
        OpenCLExecutor executor(
            "find_mismatches_packed.cl", 
            plat.get_platform(), 
            plat.get_device_ids()[device_idx]
        );

        const int OUT_BUF_SIZE = 1<<24;
        CLBuffer<match> output_buf = executor.new_clbuffer<match>(OUT_BUF_SIZE);
        CLBuffer<int> output_count = executor.new_clbuffer<int>(1);
        CLBuffer<uint64_t> genome_buf = executor.new_clbuffer<uint64_t>(B4_CHUNK_SIZE+2 + blocks_per_pattern);
        CLBuffer<uint64_t> pattern_buf = executor.new_clbuffer<uint64_t>(pattern_blocks.size());

        size_t num_pattern_blocks = (patterns.size() + LOCAL_BLOCK_SIZE - 1) / LOCAL_BLOCK_SIZE;
        CLKernel compute_results = executor.new_clkernel(
            "find_matches_packed_helper",
            CL_NDRange(B4_CHUNK_SIZE, num_pattern_blocks),
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
                next_genome_idx += GENOME_CHUNK_SIZE;
            }
            // std::cout << genome_idx << std::endl;
            if(genome_idx >= genome.size()){
                break;
            }
            std::string cur_genome = genome.substr(genome_idx, GENOME_CHUNK_SIZE+pattern_size);
            // std::cout << cur_genome << std::endl;
            std::vector<uint64_t> b4genome = make4bitpackedints(cur_genome);

            genome_buf.clear_buffer();
            output_count.clear_buffer();
            genome_buf.write_buffer(&b4genome[0], b4genome.size());
            compute_results.run();
            int out_count;
            output_count.read_buffer(&out_count, 1);
            
            if(out_count > 0){
                std::vector<match> new_matches(out_count);
                output_buf.read_buffer(&new_matches[0], out_count);
                std::vector<match> pruned_matches;
                for(match m : new_matches){
                    if(m.loc < GENOME_CHUNK_SIZE){
                        m.loc += genome_idx;
                        pruned_matches.push_back(m);
                    }
                }
                #pragma omp critical 
                {
                    matches.insert(matches.end(), pruned_matches.begin(), pruned_matches.end());
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
    std::cout << "time: " << time_spent([&](){
    find_matches_opencl_multi(genome, patterns, mismatches);
    }) << std::endl;
    return true;
}
