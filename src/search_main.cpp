#include "find_mismatches.h"
#include "read_fasta.h"
#include "timing.h"
#include "chromloc.h"
#include "postprocess.h"
#include "bit4ops.h"
#include "RangeIterator.h"
#include "search_genome.h"
#include "read_genome.h"

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
	

	Channel<GenomeMatch> data_output_channel;
	std::thread searcher_thread(search_genome, patterns, data_folder, mismatches, &data_output_channel);
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
	searcher_thread.join();
	});
	std::cerr << "Time spent: " << calc_time << std::endl;
	return 0;
}
