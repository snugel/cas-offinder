#include "test/test_framework.h"
#include "find_mismatches.h"
#include "RangeIterator.h"
#include "opencl_executor.h"
#include "oclkernels.h"
#include "bit4ops.h"
#include "timing.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <thread>
#include <cassert>
#include <immintrin.h>
#include <ctime>

using block_ty = uint32_t;
constexpr size_t bit4_c = sizeof(block_ty) * 8 / 4;


std::vector<match> find_matches(std::string & genome, std::vector<std::string> & patterns, int max_mismatches){
    std::vector<uint32_t> genomeb4 = make4bitpackedint32(genome);
    std::shared_ptr<uint32_t> ptr(new uint32_t[genomeb4.size() + patterns.at(0).size() + 1]());
    std::copy(genomeb4.begin(), genomeb4.end(), ptr.get());

    Channel<GenomeInput> input_channel;
    Channel<WorkerOutput> out_channel;
    input_channel.send(GenomeInput{.data=ptr,.size=genomeb4.size(),.idx=0});
    input_channel.terminate();

    std::thread find_matches_thread(find_matches_worker, &input_channel, patterns, max_mismatches, &out_channel);

    std::vector<match> matches;
    WorkerOutput outs;
    while(out_channel.receive(outs)){
        // filter out and output matches
        for(size_t i : range(outs.num_matches)){
            match m = outs.matches.get()[i];
            if(m.loc < outs.input.size * bit4_c){
                m.loc += outs.input.idx * bit4_c;
                matches.push_back(m);
            }
        }
    }

    find_matches_thread.join();

    return matches;
}
void find_matches_device_worker(
    Channel<GenomeInput> * genome_data_stream, 
    std::vector<block_ty> pattern_blocks, 
    OpenCLPlatform * plat,
    size_t device_idx,
    uint32_t pattern_size, 
    uint32_t num_patterns, 
    uint32_t blocks_per_pattern ,
    int max_mismatches, 
    Channel<WorkerOutput> * out_stream){

    constexpr size_t MAX_GENOME_CHUNK_SIZE = 1<<24;
    
    std::string arguments = "-Dblocks_per_pattern=" + std::to_string(blocks_per_pattern);
    OpenCLExecutor executor(
        program_src, 
        plat->get_platform(), 
        plat->get_device_ids()[device_idx],
        arguments
    );

    const int OUT_BUF_SIZE = 1<<22;
    CLBuffer<match> output_buf = executor.new_clbuffer<match>(OUT_BUF_SIZE);
    CLBuffer<int> output_count = executor.new_clbuffer<int>(1);
    CLBuffer<block_ty> genome_buf = executor.new_clbuffer<block_ty>(MAX_GENOME_CHUNK_SIZE+2 + blocks_per_pattern);
    CLBuffer<block_ty> pattern_buf = executor.new_clbuffer<block_ty>(pattern_blocks.size());

    size_t num_pattern_blocks = num_patterns;

    CLKernel compute_results = executor.new_clkernel(
        "find_matches_packed_helper",
        {
            genome_buf.k_arg(),
            pattern_buf.k_arg(),
            make_arg(max_mismatches),
            make_arg(pattern_size),
            output_buf.k_arg(),
            output_count.k_arg(),
        }
    );
    pattern_buf.clear_buffer();
    pattern_buf.write_buffer(pattern_blocks);
    output_buf.clear_buffer();
    GenomeInput input;
    while(genome_data_stream->receive(input)){
        // std::cout << cur_genome << std::endl;

        genome_buf.clear_buffer();// TODO: check if removing this improves performance 
        output_count.clear_buffer();
        genome_buf.write_buffer(input.data.get(),input.size+blocks_per_pattern+1);
        compute_results.run(
            CL_NDRange(input.size, num_pattern_blocks),
            CL_NDRange(),
            CL_NDRange()
        );
        int out_count;
        output_count.read_buffer(&out_count, 1);
        
        if(out_count > 0){
            std::shared_ptr<match> out_matches(new match[out_count]);
            output_buf.read_buffer(out_matches.get(), out_count);
            out_stream->send(WorkerOutput{
                .input=input,
                .matches=out_matches,
                .num_matches=size_t(out_count),
            });
        }
    }

}

void find_matches_worker(Channel<GenomeInput> * genome_data_stream, std::vector<std::string> patterns, int max_mismatches, Channel<WorkerOutput> * out_stream)
{
    uint32_t pattern_size = patterns.at(0).size();
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
    std::vector<std::thread> threads;
    for(size_t dev_idx : range(plat.num_cl_devices())){
        threads.emplace_back(
            find_matches_device_worker,
            genome_data_stream,
            pattern_blocks,
            &plat,
            dev_idx,
            pattern_size,
            patterns.size(),
            blocks_per_pattern,
            max_mismatches,
            out_stream
        );
    }

    for(size_t dev_idx : range(plat.num_cl_devices())){
        threads[dev_idx].join();
    }

    out_stream->terminate();
}
