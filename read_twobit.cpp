#include "read_twobit.h"

using namespace std;

unsigned int read_uint(ifstream *input) {
	unsigned int val = 0;
	input->read(reinterpret_cast<char*>(val), sizeof(unsigned int));
	return val;
}

int read_twobit(string filepath, vector<string> *chrnames, string *content, vector<unsigned long long> *chrpos) {
	ifstream input;
	input.open(filepath.c_str());

	if (read_uint(&input) != 440477507) return 1;
	if (read_uint(&input) != 0) return 1;

	unsigned int chrcnt = read_uint(&input);
	input.seekg(4, input.cur); // Reserved

	return 1;
}