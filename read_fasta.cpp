#include "read_fasta.h"

#ifndef _WIN32
  #include <unistd.h> // readlink
  #include <limits.h> // PATH_MAX
  #include <cstring> // memset
#endif

using namespace std;

int read_fasta(string &filepath, vector<string> &chrnames, string &content, vector<unsigned long long> &chrpos) {
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
	chrnames.clear();
	content.clear();
	chrpos.clear();
	char c; input.get(c);
	if (c != '>') {
		input.close();
		return 1;
	}
	input.seekg(0, input.beg);
	while (getline(input, line)){
		if (!line.empty()) {
			if (line[line.length()-1] == '\r')
				line = line.substr(0, line.length()-1);
			if (line[0] == '>') { // Identifier marker
				name = line.substr(1);
				//if (((cnt++) % 10000) == 0) cerr << "Reading " << name << endl;
				chrnames.push_back(name);
				if (chrpos.size() != 0) content += ";"; // seperator
				chrpos.push_back(content.size());
			}
			else {
				transform(line.begin(), line.end(), line.begin(), ::toupper);
				try {
					content += line;
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
