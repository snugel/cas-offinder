#include "find_mismatches.h"
#include "read_fasta.h"
#include "timing.h"

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
	if(argc != 5){
		std::cerr << "Needs 4 command line arguments: 4bit_path, pattern_path, mismatches, device\n";
		exit(1);
	}
	double calc_time = time_spent([&](){
	std::string bit4_path(argv[1]);
	std::string pattern_path(argv[2]);
	int mismatches = atoi(argv[3]);
	char device = *argv[4];

	std::cerr << "Reading patterns..." << std::endl;

	std::vector<std::string> patterns; 
	std::string line;
	std::ifstream pattern_file(pattern_path);
	while (std::getline(pattern_file, line))
	{
		if(line.size()){
			patterns.push_back(line);
		}
	}

	std::cerr << "Reading genome..." << std::endl;

	ifstream bit4_file(bit4_path, ios::binary);
	const auto begin = bit4_file.tellg();
	bit4_file.seekg (0, ios::end);
	const auto end = bit4_file.tellg();
	const auto fsize = (end-begin);
	bit4_file.seekg (0, ios::beg);

	assert(fsize % 4 == 0);
	vector<uint32_t> data(fsize/4);
	bit4_file.read((char*)&data[0], fsize);
    std::cerr << "Aprox genome size: " << fsize*2 << std::endl;
	
	std::cerr << "Searching genome..." << std::endl;
    
    std::cout << "Idx\tLocation\tMismatches\n";
	find_matches(data, patterns, mismatches, [](match m){
        atomic_print_match(m);
    });

	});
	std::cerr << "Time spent: " << calc_time << std::endl;
	return 0;
}
