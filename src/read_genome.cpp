#include "read_genome.h"
#include <iostream>
#include <fstream>
#include <thread>
#include <cassert>

using namespace std;

void read_genome(std::string data_folder, size_t pattern_size, Channel<GenomeInput> * data_output_p){
	std::cerr << "Reading genome..." << std::endl;
	constexpr size_t block_size = 1<<24;
	size_t buffer_padding = (pattern_size+7)/8 + 1;


	ifstream bit4_file(data_folder + "genome.4bit", ios::binary);
	const auto begin = bit4_file.tellg();
	bit4_file.seekg (0, ios::end);
	const auto end = bit4_file.tellg();
	const auto fsize = (end-begin);
	bit4_file.seekg (0, ios::beg);

	std::this_thread::sleep_for(std::chrono::milliseconds(500));

    std::cerr << "Aprox genome size: " << fsize*2 << std::endl;

	assert(fsize % 4 == 0);
	uint64_t num_blocks = fsize/4;
	std::vector<uint32_t> prev_padding(buffer_padding);
	bit4_file.read((char*)prev_padding.data(), buffer_padding*sizeof(uint32_t));

	for(size_t i = 0; i < num_blocks; i += block_size){
		size_t this_block_size = std::min(block_size, num_blocks-i);
		size_t padded_size = this_block_size + buffer_padding;
		std::shared_ptr<uint32_t> data(new uint32_t[padded_size]());
		std::copy(prev_padding.begin(), prev_padding.end(), data.get());	
		bit4_file.read((char*)data.get() + buffer_padding, this_block_size*sizeof(uint32_t));
		data_output_p->send(GenomeInput{.data=data,.size=this_block_size,.idx=i});
		std::copy(data.get()+this_block_size, data.get() + padded_size, prev_padding.begin());	
		// bit4_file.seekg (i*sizeof(uint32_t), ios::beg);
		std::cerr << "read block\n";
	}
	data_output_p->terminate();
}