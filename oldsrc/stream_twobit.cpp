#include "stream_twobit.h"

#include <iostream>
#include <utility>

#include <cstdlib>
#include <cstring> // memset

#ifdef _WIN32
#include <assert.h>
#else
#include <unistd.h> // readlink
#include <limits.h> // PATH_MAX
#endif

using namespace std;

inline uint32_t read_uint(FILE *input) {
	uint32_t val;
	fread((char*)&val, 4, 1, input);
	return val;
}


int get_twobit_chromlocs(std::string filepath, std::vector<chromloc> & out_locs, uint64_t & cur_pos){
	FILE *input = fopen(filepath.c_str(), "rb");

	if (read_uint(input) != 440477507) { // Magic
		fclose(input);
		return 1;
	}
	if (read_uint(input) != 0) { // Version should be 0
		fclose(input);
		return 1;
	}

	uint32_t chrcnt = read_uint(input);
	fseek(input, 4, SEEK_CUR); // Reserved

    char len_chrname;
    size_t readsize;
    char chrname[256+1];
    std::vector<std::string> chrnames;
    std::vector<uint64_t> chrpos;
    std::vector<uint64_t> filepos;
	for (uint32_t i=0; i<chrcnt; i++) {
		readsize = fread(&len_chrname, 1, 1, input);
		readsize = fread(chrname, 1, len_chrname, input);
		chrname[len_chrname] = 0;
		chrnames.push_back(string(chrname));
        filepos.push_back(read_uint(input)); // Absolute position of each sequence
	}

	chrpos.push_back(cur_pos);
	for (uint32_t i=0; i<chrcnt; i++) {
    	fseek(input, filepos[i], SEEK_CUR); // go to start of chromosome
		uint32_t chrlen = read_uint(input);
        cur_pos += chrlen;
        if(i != chrcnt-1){
            chrpos.push_back(chrlen);
        }
	}
	fclose(input);
	return 0;
}
