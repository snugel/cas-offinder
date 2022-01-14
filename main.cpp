#include "find_mismatches.h"
#include "read_fasta.h"
#include "read_twobit.h"
#include "timing.h"
#include "postprocess.h"

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
#include <cstdarg>

using namespace std;

int main(int argc, char *argv[]) {
	if(argc != 5){
		std::cerr << "Needs 4 command line arguments: 2bit_path, pattern_path, mismatches, device\n";
		exit(1);
	}
	double calc_time = time_spent([&](){
	std::string bit2_path(argv[1]);
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

	vector<string> m_chrnames;
	string m_chrdata;
	vector<uint64_t> m_chrpos;
	if (read_twobit(bit2_path, m_chrnames, m_chrdata, m_chrpos)) {
		cerr << "Non-acceptable file: " << bit2_path << endl;
		exit(1);
	}
	
	std::cerr << "Genome size: " << m_chrdata.size() << std::endl;
	std::cerr << "Searching genome..." << std::endl;
	std::vector<match> matches = find_matches(m_chrdata, patterns, mismatches);
	std::cerr << "Postprocessing..." << std::endl;
	filter_chrom_walls(matches, m_chrpos, patterns[0].size());
	sort_matches(matches);

	std::vector<uint64_t> rel_locs;
	std::vector<uint64_t> chrom_idxs;
	get_chrom_info(
		matches, 
		m_chrpos,
		rel_locs,
		chrom_idxs
	);
    std::cout << "Idx\tLocation\tMismatches\n";
    for(match m : matches){
        std::cout << m.pattern_idx << '\t' << m.loc << '\t' << m.mismatches << '\n';
    }
	});
	std::cerr << "Time spent: " << calc_time << std::endl;
	return 0;
}
