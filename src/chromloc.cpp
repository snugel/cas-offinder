#include "chromloc.h"
#include "RangeIterator.h"
#include <algorithm>

std::vector<chromloc> parse_chromloc_file(std::istream & stream){
    std::vector<chromloc> locs;
    std::vector<std::string> lines;
    std::string line;
	while (getline(stream, line)){
        lines.push_back(line);
    }
    if(lines.size() == 0 || lines[0] != "Name\tLoc"){
        throw std::runtime_error("not a valid chrom loc file");
    }
    uint64_t last_loc = 0;
    for(size_t i : range(1, lines.size())){
        std::string line = lines[i];
        if (std::count(line.begin(), line.end(), '\t') != 1){
            throw std::runtime_error("not a valid chrom loc file, needs one tab per line");        
        }
        size_t sep = line.find_first_of('\t');
        std::string name = line.substr(0, sep);
        uint64_t loc = stoll(line.substr(sep+1));
        if(loc < last_loc){
            throw std::runtime_error("not a valid chrom loc file, locations out of order");        
        }
        last_loc = loc;
        locs.push_back(chromloc{
            .name=name,
            .loc=loc,
        });
    }
	return locs;
}
void write_chromloc_file(std::vector<chromloc> & locs, std::ostream & chrom_locs){
    chrom_locs << "Name\tLoc\n";
    for(chromloc loc : locs){
        chrom_locs << loc.name << "\t" << loc.loc << "\n";
    }
}
