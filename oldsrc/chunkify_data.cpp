#include "chunkify_data.h"
#include "timing.h"
#include "bit4ops.h"
#include <cassert>
#include <algorithm>
#include <iostream>

double total_time_g = 0;
void add_chunk(std::string & cur_chunk, Channel<GenomeInput> * genome_input, size_t & position){
    std::vector<uint32_t> b32_data;
    total_time_g += time_spent([&](){
    clean_bogus(cur_chunk);
    b32_data = make4bitpackedint32(cur_chunk);
    });
    std::shared_ptr<uint32_t> inpt(new uint32_t[b32_data.size()]);
    std::copy(b32_data.begin(), b32_data.end(), inpt.get());
    genome_input->send(GenomeInput{
        .data=inpt,
        .size=b32_data.size(),
        .idx=position,
    });
    std::cerr << "created chunk:\n";
}

void chunkify_data(Channel<FileChunk> * input_channel, Channel<GenomeInput> * genome_input, int pattern_size, int chunk_size){
    std::string buffer;
    std::string next;
    size_t position = 0;
    assert(chunk_size%8 == 0);
    double total_time = 0;
    while(input_channel->receive(next)){
        total_time += time_spent([&](){
	buffer += next;
        int64_t bidx = 0;
        for(;bidx <= int64_t(buffer.size() - pattern_size - chunk_size); bidx += chunk_size){
            std::string cur_chunk(buffer.begin()+bidx, buffer.begin()+bidx + chunk_size+pattern_size);
            add_chunk(cur_chunk, genome_input, position);
            position += chunk_size;
        }
        // bidx -= chunk_size;
        if(bidx > 0){
            buffer.erase(buffer.begin(), buffer.begin() + bidx);
        }});
    }
    if(buffer.size() > 0){
        int64_t bidx = 0;
        for(;bidx < int64_t(buffer.size() - pattern_size); bidx += chunk_size){
            size_t end = std::min(size_t(bidx + chunk_size), buffer.size() - pattern_size);
            std::string cur_chunk(buffer.begin()+bidx, buffer.begin()+end + pattern_size);
            add_chunk(cur_chunk, genome_input, position);
            position += end - bidx;
        }
    }
    std::cerr << "chunkify inner_time: " << total_time_g << "\n";
    std::cerr << "chunkify time: " << total_time << "\n";
    genome_input->terminate();
}
