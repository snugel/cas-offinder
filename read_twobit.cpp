#include "read_twobit.h"

using namespace std;

typedef union {
	char bytes[4];
	unsigned int num;
} record;

unsigned int swapEndian(unsigned int num) {
	unsigned int b0, b1, b2, b3;
	b0 = (num & 0x000000ff) << 24u;
	b1 = (num & 0x0000ff00) << 8u;
	b2 = (num & 0x00ff0000) >> 8u;
	b3 = (num & 0xff000000) >> 24u;
	return (b0 | b1 | b2 | b3);
}

unsigned int read_uint(ifstream &input) {
	record val;
	input >> val.bytes;
	return swapEndian(val.num);
}

int read_twobit(string &filepath, vector<string> *chrnames, string *content, vector<unsigned long long> *chrpos) {
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