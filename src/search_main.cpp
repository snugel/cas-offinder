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

#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>

#include <cstdlib>
#include <cstring>
#include <cassert>
#include <cstdarg>

using namespace std;


int main(int argc, char *argv[]) {
	constexpr size_t num_pos_args = 4;
	if(argc > num_pos_args){
		std::cerr << "Needs at least 5 command line arguments: data_folder, mismatches, device, queries\n";
		exit(1);
	}
	double calc_time = time_spent([&](){

	char * data_folder(argv[1]);
	int mismatches = atoi(argv[2]);
	char device = *argv[3];
	char ** pattern_chrs = &argv[4];
	int num_patterns = argc - num_pos_args;

	std::vector<std::string> patterns(num_patterns); 
	for(size_t i : range(num_patterns)){
		patterns[i] = std::string(pattern_chrs[num_patterns]);
	}
	search_genome(patterns, mismatches, data_folder);

	});
	std::cerr << "Time spent: " << calc_time << std::endl;
	return 0;
}
