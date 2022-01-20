#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include "opencl_executor.h"
#include "bit4ops.h"
#include "timing.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <cassert>
#include <immintrin.h>
#include <ctime>

using block_ty = uint32_t;
constexpr size_t bit4_c = sizeof(block_ty) * 8 / 4;


std::vector<match> find_matches(std::string & genome, std::vector<std::string> & patterns, int max_mismatches){
    std::vector<uint32_t> genomeb4 = make4bitpackedint32(genome);
    std::vector<match> matches;
    find_matches_gold(genomeb4, patterns, max_mismatches, [&](match m){
        #pragma omp critical 
        {
        matches.push_back(m);
        }
    });
    return matches;
}

void find_matches(std::vector<uint32_t> & genomeb4, std::vector<std::string> & patterns, int max_mismatches, std::function<void(match)> func)
{
    if (patterns.size() == 0)
    {
        return;
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

    OpenCLPlatform plat;
    volatile size_t next_genome_idx = 0;
    constexpr size_t GENOME_CHUNK_SIZE = 1<<24;
    #pragma omp parallel for
    for(size_t device_idx = 0; device_idx < plat.num_cl_devices(); device_idx++){
        OpenCLExecutor executor(
            "find_mismatches_packed.cl", 
            plat.get_platform(), 
            plat.get_device_ids()[device_idx]
        );

        const int OUT_BUF_SIZE = 1<<22;
        CLBuffer<match> output_buf = executor.new_clbuffer<match>(OUT_BUF_SIZE);
        CLBuffer<int> output_count = executor.new_clbuffer<int>(1);
        CLBuffer<block_ty> genome_buf = executor.new_clbuffer<block_ty>(GENOME_CHUNK_SIZE+2 + blocks_per_pattern);
        CLBuffer<block_ty> pattern_buf = executor.new_clbuffer<block_ty>(pattern_blocks.size());

        size_t num_pattern_blocks = patterns.size();
        CLKernel compute_results = executor.new_clkernel(
            "find_matches_packed_helper",
            CL_NDRange(GENOME_CHUNK_SIZE, num_pattern_blocks),
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
            if(genome_idx >= genomeb4.size()){
                break;
            }
            // std::cout << cur_genome << std::endl;

            genome_buf.clear_buffer();
            output_count.clear_buffer();
            genome_buf.write_buffer(&genomeb4[genome_idx], std::min(genomeb4.size() - blocks_per_pattern - genome_idx, GENOME_CHUNK_SIZE) + blocks_per_pattern);
            compute_results.run();
            int out_count;
            output_count.read_buffer(&out_count, 1);
            
            if(out_count > 0){
                std::vector<match> new_matches(out_count);
                output_buf.read_buffer(&new_matches[0], out_count);
                #pragma omp critical 
                {
                // filter out and output matches
                for(match m : new_matches){
                    if(m.loc < GENOME_CHUNK_SIZE * bit4_c){
                        m.loc += genome_idx * bit4_c;
                        func(m);
                    }
                }
                }
            }
        }
    }
}
