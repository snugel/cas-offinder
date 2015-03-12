#include "read_twobit.h"

using namespace std;

typedef union {
	char bytes[4];
	unsigned int num;
} record;

typedef union {
	unsigned char bytes[4];
	unsigned int num;
} RECORD4;

unsigned int read_uint(ifstream &input) {
	static RECORD4 val4;
	input.read((char*)val4.bytes, 4);
	return val4.num;
}

char bit_to_seq(unsigned char b) {
	switch(b) {
	case 0:
		return 'T';
	case 1:
		return 'C';
	case 2:
		return 'A';
	case 3:
		return 'G';
}

int read_twobit(string &filepath, vector<string> &chrnames, string &content, vector<unsigned long long> &chrpos) {
	ifstream input(filepath.c_str(), ios::in|ios::binary);
	if (read_uint(input) != 440477507) {
		input.close();
		return 1;
	}
	if (read_uint(input) != 0) {
		input.close();
		return 1;
	}

	unsigned int chrcnt = read_uint(input);
	input.seekg(4, input.cur); // Reserved
	input.close();
	return 1;
}
