#include "find_mismatches.h"
#include "read_fasta.h"
#include "timing.h"
#include "chromloc.h"
#include "postprocess.h"
#include "bit4ops.h"
#include "RangeIterator.h"
#include "search_genome.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <thread>

#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdarg>

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

int main(int argc, char *argv[]) {
	constexpr size_t num_pos_args = 5;
	if(argc < num_pos_args){
		std::cerr << "Needs at least 4 command line arguments: data_folder, mismatches, device, queries\n";
		exit(1);
	}
	double calc_time = time_spent([&](){

	std::string data_folder(argv[1]);
	int mismatches = atoi(argv[2]);
	char device = *argv[3];
	char ** pattern_chrs = &argv[4];
	int num_patterns = argc - num_pos_args+1;

	std::cerr << "Number of patterns: " << num_patterns << std::endl;
	std::vector<std::string> patterns(num_patterns); 
	for(size_t i : range(num_patterns)){
		patterns[i] = std::string(pattern_chrs[i]);
	}
	
	if(data_folder.back() != '/'){
		data_folder.push_back('/');
	}

	std::ifstream chrom_locs_file(data_folder + "chrom_locs.csv");
	std::vector<chromloc> chrom_locs = parse_chromloc_file(chrom_locs_file);

	Channel<GenomeInput> data_input_channel(6);
	Channel<GenomeMatch> data_output_channel;
	std::thread reader_thread(read_genome, data_folder, patterns.at(0).size(), &data_input_channel);
	std::thread searcher_thread(search_genome, patterns, chrom_locs, mismatches, &data_input_channel, &data_output_channel);
	// search_genome(patterns, mismatches, data_folder);

    std::cout << "Idx\tChromosome\tLocation\tDNA\tMismatches\n";

	GenomeMatch m;
	while(data_output_channel.receive(m)){
		std::cout 
			<< m.pattern_idx << "\t"
			<< m.cromosome << "\t"
			<< m.chrom_loc << "\t"
			<< m.dna_match << "\t"
			<< m.mismatches << "\n";
	}
	reader_thread.join();
	searcher_thread.join();
	});
	std::cerr << "Time spent: " << calc_time << std::endl;
	return 0;
}
