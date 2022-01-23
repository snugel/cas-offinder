
#include "search_genome.h"
#include "RangeIterator.h"
#include "chromloc.h"
#include "postprocess.h"
#include "find_mismatches.h"
#include "bit4ops.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include <vector>

using namespace std;

void search_genome(std::vector<std::string> & patterns, int mismatches, std::string data_folder){
    
	if(data_folder.back() != '/'){
		data_folder.push_back('/');
	}
    size_t pattern_size = patterns.at(0).size();

	std::cerr << "Reading genome..." << std::endl;

	ifstream bit4_file(data_folder + "genome.4bit", ios::binary);
	const auto begin = bit4_file.tellg();
	bit4_file.seekg (0, ios::end);
	const auto end = bit4_file.tellg();
	const auto fsize = (end-begin);
	bit4_file.seekg (0, ios::beg);

	assert(fsize % 4 == 0);
	vector<uint32_t> data(fsize/4);
	bit4_file.read((char*)&data[0], fsize);
    std::cerr << "Aprox genome size: " << fsize*2 << std::endl;

    std::ifstream chrom_locs_file(data_folder + "chrom_locs.csv");
	std::vector<chromloc> chrom_locs = parse_chromloc_file(chrom_locs_file);
	std::vector<uint64_t> chrom_poses(chrom_locs.size());
	for(size_t i : range(chrom_locs.size())){
		chrom_poses[i] = chrom_locs[i].loc;
	}

	std::cerr << "Searching genome..." << std::endl;
    
    std::cout << "Idx\tChromosome\tLocation\tDNA\tMismatches\n";
	find_matches(data, patterns, mismatches, [&](match m){
		if(!crosses_chrom_wall(m.loc, chrom_poses, pattern_size)){
			chrom_info chr_info = get_chrom_info(m.loc, chrom_poses);
			chromloc chrloc = chrom_locs.at(chr_info.chrom_idx);
	        std::string sstr = (
				std::to_string(m.pattern_idx) + "\t" + 
				chrloc.name + "\t" + 
				std::to_string(chr_info.rel_loc) + "\t" + 
				bit4tostr(data, m.loc, m.loc + pattern_size) + "\t" +
				std::to_string(m.mismatches) + "\n" 
			);
			std::cout << sstr;
		}
    });

}
