#include "chunkify_data.h"
#include "bit4ops.h"
#include "test/test_framework.h"
#include "RangeIterator.h"
#include <thread>
#include <iostream>

using namespace std;

constexpr size_t bit4_c = sizeof(uint32_t)*2;

static FileChunk to_chunk(string nextchunk){
    std::shared_ptr<char[]> arr(new char[nextchunk.size()]);
    std::copy(nextchunk.begin(), nextchunk.end(), arr.get()); 
    return FileChunk{.data=arr, .size=nextchunk.size()};
}


TEST(test_chunkify_data){
    std::vector<std::string> input_strings = {
        "ATCGAC",
        "CTAATCGACTGCAGTCGATCGATAT",
        "ATATCAGTCTGAGTCGACC",
        "ATATCAGTCTGAGTCGACC",
    };
    std::string concat_data;// = "ATCGACCTAATCGACTGCAGTCGATCGATATATATCAGTCTGAGTCGACC";
    for(string s : input_strings){
        concat_data += s;
    }
    int pattern_size = 3;
    int chunk_size = 32;
    vector<vector<uint32_t>> expected_datas = {
        make4bitpackedint32(concat_data.substr(0,(chunk_size+pattern_size))),
        make4bitpackedint32(concat_data.substr(chunk_size,(2*chunk_size+pattern_size))),
        make4bitpackedint32(concat_data.substr(2*chunk_size)),
    };
    vector<uint32_t> posses = {0, 32, 64};
    Channel<FileChunk> input_channel;
    Channel<GenomeInput> out_channel;
    thread chunkfy_t(chunkify_data, &input_channel, &out_channel, pattern_size, chunk_size);
    for(string s : input_strings){
        input_channel.send(to_chunk(s));
    }
    input_channel.terminate();
    chunkfy_t.join();
    GenomeInput input;
    size_t i = 0;
    bool result = true;
    while(out_channel.receive(input)){

        for(auto x : expected_datas[i]){
            std::cout << x <<" ";
        }
        cout << "\n";
        for(auto i : range(input.size)){
            std::cout << input.data.get()[i] <<" ";
        }
        cout << "\n";
        
        if(input.size != expected_datas[i].size() || !std::equal(expected_datas[i].begin(), expected_datas[i].end(), input.data.get())){
            std::cout << "false\n";
            result = result && false;
        }
        else{
            std::cout << "false\n";

        }
        // std::cout << posses[i] << "\t" << expected_datas[i] << std::endl;
        // std::cout << input.idx << "\t" << std::string(input.data.get(),input.data.get()+input.size) << std::endl;
        i += 1;
    }

    return result;
}
