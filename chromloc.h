#include <string>
#include <istream>
#include <ostream>
#include <vector>
#include <inttypes.h>

struct chromloc{
	std::string name;
	uint64_t loc;
};
std::vector<chromloc> parse_chromloc_file(std::istream & stream);
void write_chromloc_file(std::vector<chromloc> & locs, std::ostream & stream);
