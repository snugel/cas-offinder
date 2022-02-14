#include "read_fasta.h"
#include <iostream>
#include <algorithm>
#include <fstream>
#include <string>

#ifndef _WIN32
  #include <unistd.h> // readlink
  #include <limits.h> // PATH_MAX
  #include <cstring> // memset
#endif

using namespace std;

int read_fasta(std::string &filepath, std::vector<std::string> &chrnames, Channel<std::string> &content, std::vector<uint64_t> &chrpos){
	string line, name;
	ifstream input;
#ifdef _WIN32
	input.open(filepath.c_str());
#else
	char path_buf[PATH_MAX+1]; memset(path_buf, 0, PATH_MAX + 1);
	int path_cnt = readlink(filepath.c_str(), path_buf, PATH_MAX);
	if (path_cnt >= 0)
		input.open(path_buf);
	else
		input.open(filepath.c_str());
#endif
	char c; input.get(c);
	if (c != '>') {
		input.close();
		return 1;
	}
	input.seekg(0, input.beg);
	string nextchunk;
	while (getline(input, line)){
		if (!line.empty()) {
			if (line[line.length()-1] == '\r')
				line = line.substr(0, line.length()-1);
			if (line[0] == '>') { // Identifier marker
				name = line.substr(1);
				//if (((cnt++) % 10000) == 0) cerr << "Reading " << name << endl;
				chrnames.push_back(name);
				//if (chrpos.size() != 0) content += i"; // seperator
				if(chrpos.size() == 0){
					chrpos.push_back(0);
				}
				else{
					chrpos.push_back(chrpos.back() + nextchunk.size());
				}
				nextchunk.resize(0);
			}
			else {
				transform(line.begin(), line.end(), line.begin(), ::toupper);
				try {
					nextchunk += line;
				}
				catch (const std::bad_alloc&) {
					cerr << "File is too big!" << endl << "Split FASTA file, or try 64bit version of Cas-OFFinder." << endl;
					exit(1);
				}
			}
		}
	}
	input.close();
	return 0;
}
