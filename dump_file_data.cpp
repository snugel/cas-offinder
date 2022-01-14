#include "read_twobit.h"
#include "timing.h"
#include "RangeIterator.h"
#include "bit4ops.h"

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cassert>

#include <sys/time.h>
#include <unistd.h>
#include <sys/stat.h>

#include <cstdlib>
#include <cstring>
#include <cstdarg>

using namespace std;

int main(int argc, char *argv[]) {
	if(argc != 3){
		std::cerr << "Needs 2 command line arguments: 2bit_path, out_folder\n";
		exit(1);
	}
	std::string bit2_path(argv[1]);
	std::string out_folder(argv[2]);

	std::cerr << "Reading genome..." << std::endl;

	vector<string> m_chrnames;
	string m_chrdata;
	vector<uint64_t> m_chrpos;
	if (read_twobit(bit2_path, m_chrnames, m_chrdata, m_chrpos)) {
		cerr << "Non-acceptable file: " << bit2_path << endl;
		exit(1);
	}
    // for(size_t i : range(m_chrnames.size())){
    // 	std::cerr << m_chrnames[i] << "\t" << m_chrpos[i] << std::endl;
    // }
	
	std::cerr << "Genome size: " << m_chrdata.size() << std::endl;

	std::cerr << "Cleaning genome data..." << std::endl;

    clean_bogus(m_chrdata);
    std::vector<uint32_t> bit4_data = make4bitpackedint32(m_chrdata);

	std::cerr << "Creating output directory..." << std::endl;

    if (out_folder.back() != '/'){
        out_folder.push_back('/');
    }
    if (mkdir(out_folder.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == -1)
    {
        if( errno == EEXIST ) {
            // alredy exists
        } else {
            // something else
            std::cerr << "cannot create folder error:" << strerror(errno) << std::endl;
            throw std::runtime_error( strerror(errno) );
        }
    }
    assert(m_chrpos.size() == m_chrnames.size());
    std::ofstream chrom_locs(out_folder + "chrom_locs.csv");
    chrom_locs << "Name\tLoc\n";
    for(size_t i : range(m_chrpos.size())){
        chrom_locs << m_chrnames[i] << "\t" << m_chrpos[i] << "\n";
    }

	std::cerr << "Dumping genome data..." << std::endl;
    
    std::ofstream genome_out(out_folder + "genome.4bit", std::ios::binary);
    genome_out.write((char*)&bit4_data[0], bit4_data.size() * sizeof(bit4_data[0]));

	return 0;
}
