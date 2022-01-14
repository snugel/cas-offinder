#include "find_mismatches.h"
#include "read_fasta.h"
#include "read_twobit.h"

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
		return 1;
	}
	std::string bit2_path(argv[1]);
	std::string pattern_path(argv[2]);
	int mismatches = atoi(argv[3]);
	char device = *argv[4];

	// parse patterns
	std::vector<std::string> patterns; 
	std::string line;
	std::ifstream pattern_file(pattern_path);
	while (std::getline(pattern_file, line))
	{
		if(line.size()){
			patterns.push_back(line);
		}
	}
		std::cerr << "p" << std::endl;
	for(auto p : patterns){
		std::cerr << p << std::endl;
	}

		std::cerr << "p" << std::endl;
	vector<string> m_chrnames;
	string m_chrdata;
	vector<unsigned long long> m_chrpos;
	
	// parse genome
	if (read_twobit(bit2_path, m_chrnames, m_chrdata, m_chrpos)) {
		cerr << "Non-acceptable file: " << bit2_path << endl;
		return 1;
	}
	std::vector<match> matches = find_matches(m_chrdata, patterns, mismatches);
    std::cout << "Idx\tLocation\tMismatches\n";
    for(match m : matches){
        std::cout << m.pattern_idx << '\t' << m.loc << '\t' << m.mismatches << '\n';
    }
}
