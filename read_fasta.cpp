#include "read_fasta.h"

using namespace std;

int read_fasta(string &filepath, vector<string> &chrnames, string &content, vector<unsigned long long> &chrpos) {
	string line, name;
	ifstream input;
	input.open(filepath.c_str());
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
				//if (((cnt++) % 10000) == 0) cout << "Reading " << name << endl;
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
					cout << "File is too big!" << endl << "Split FASTA file, or try 64bit version of Cas-OFFinder." << endl;
					exit(1);
				}
			}
		}
	}
	input.close();
	return 0;
}
