#include "read_folder.h"
#include "timing.h"
#include "RangeIterator.h"
#include "bit4ops.h"
#include "chromloc.h"

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
		std::cerr << "Needs 2 command line arguments: path, out_folder\n";
		exit(1);
	}
	std::string path(argv[1]);
	std::string out_folder(argv[2]);

	std::cerr << "Reading genome..." << std::endl;

	vector<string> m_chrnames;
	string m_chrdata;
	vector<uint64_t> m_chrpos;
	int cnum = 0;
	if (read_folder(path, m_chrnames, m_chrdata, m_chrpos)) {
        cerr << "Non-acceptable file/folder: " << path << endl;
		exit(1);
	}
	
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
	std::vector<chromloc> locs(m_chrpos.size());
	for(size_t i : range(m_chrpos.size())){
		locs[i] = chromloc{
			.name=m_chrnames[i],
			.loc=m_chrpos[i],
		};
	}
    std::ofstream chrom_locs(out_folder + "chrom_locs.csv");
	write_chromloc_file(locs, chrom_locs);

	std::cerr << "Dumping genome data..." << std::endl;
    
    std::ofstream genome_out(out_folder + "genome.4bit", std::ios::binary);
    genome_out.write((char*)&bit4_data[0], bit4_data.size() * sizeof(bit4_data[0]));

	return 0;
}
