#include "read_twobit.h"

#include <iostream>
#include <utility>

#include <cstdlib>

using namespace std;

typedef union {
	unsigned char bytes[4];
	unsigned int num;
} RECORD4;

unsigned int read_uint(FILE *input) {
	static RECORD4 val4;
	fread((char*)val4.bytes, 4, 1, input);
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
}

int read_twobit(string &filepath, vector<string> &chrnames, string &content, vector<unsigned long long> &chrpos) {
	unsigned int i, j, k, chrcnt, chrlen, nblockcnt, maskblockcnt, rem;
	unsigned char achar;
	char len_chrname;
	char chrname[256];
	vector<unsigned int> nblockstarts;
	vector<unsigned int> nblocksizes;
	
	FILE *input = fopen(filepath.c_str(), "r");
	if (read_uint(input) != 440477507) { // Magic
		fclose(input);
		return 1;
	}
	if (read_uint(input) != 0) { // Version should be 0
		fclose(input);
		return 1;
	}

	chrcnt = read_uint(input);
	fseek(input, 4, SEEK_CUR); // Reserved

	for (i=0; i<chrcnt; i++) {
		fread(&len_chrname, 1, 1, input);
		fread(chrname, 1, len_chrname, input);
		chrname[len_chrname] = 0;
		chrnames.push_back(string(chrname));
		fseek(input, 4, SEEK_CUR); // Absolute position of each sequence
	}

	chrpos.push_back(0);
	for (i=0; i<chrcnt; i++) {
		chrlen = read_uint(input);
		nblockcnt = read_uint(input);
		for (j=0; j<nblockcnt; j++) nblockstarts.push_back(read_uint(input));
		for (j=0; j<nblockcnt; j++) nblocksizes.push_back(read_uint(input));
		maskblockcnt = read_uint(input);
		fseek(input, maskblockcnt*8 + 4, SEEK_CUR);
		char *chrbuf = new char[chrlen+1]; chrbuf[chrlen] = 0;
		rem = chrlen&3;
		for (j=0; j<chrlen-rem; j+=4) {
			fread((char*)&achar, 1, 1, input);
			for (k=0; k<4; k++)
				chrbuf[j+k] = bit_to_seq((achar>>((3-k)*2))&0x3);
		}
		fread((char*)&achar, 1, 1, input);
		for (j=0; j<rem; j++)
			chrbuf[chrlen-rem+j] = bit_to_seq((achar>>((3-j)*2))&0x3);
		for (j=0; j<nblockcnt; j++)
			for (k=nblockstarts[j]; k<nblockstarts[j]+nblocksizes[j]; k++)
				chrbuf[k] = 'N';
		content += chrbuf;
		delete [] chrbuf;
		if (i < chrcnt-1) {
			content += ';';
			chrpos.push_back(content.size());
		}
	}
	fclose(input);
	return 0;
}
