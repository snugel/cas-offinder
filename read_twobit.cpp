#include "read_twobit.h"

#include <iostream>
#include <utility>

#include <cstdlib>
#include <cstring> // memset

#ifndef _WIN32
#include <unistd.h> // readlink
#include <limits.h> // PATH_MAX
#endif

using namespace std;

typedef union {
	unsigned char bytes[4];
	unsigned int num;
} RECORD4;

inline unsigned int read_uint(FILE *input) {
	static RECORD4 val4;
	static size_t readsize;
	readsize = fread((char*)val4.bytes, 4, 1, input);
	return val4.num;
}

inline char bit_to_seq(unsigned char b) {
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
	return 0;
}

int read_twobit(string &filepath, vector<string> &chrnames, string &content, vector<unsigned long long> &chrpos) {
	unsigned int i, k, chrcnt, chrlen, nblockcnt, maskblockcnt, rawlen, rem, cnt;
#ifdef _MSC_VER
	int j;
#else
	unsigned int j;
#endif
	size_t readsize;
	char len_chrname;
	char chrname[256];
	vector<unsigned int> nblockstarts;
	vector<unsigned int> nblocksizes;

#ifdef _WIN32
	FILE *input = fopen(filepath.c_str(), "rb");
#else
	FILE *input;
	char path_buf[PATH_MAX + 1]; memset(path_buf, 0, PATH_MAX + 1);
	int path_cnt = readlink(filepath.c_str(), path_buf, PATH_MAX);
	if (path_cnt >= 0)
		input = fopen(path_buf, "rb");
	else
		input = fopen(filepath.c_str(), "rb");
#endif

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
		readsize = fread(&len_chrname, 1, 1, input);
		readsize = fread(chrname, 1, len_chrname, input);
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

		rem = chrlen&3;
		rawlen = chrlen/4+(rem==0?0:1);
		char *chrbuf = new char[chrlen+1]; chrbuf[chrlen] = 0;
		char *raw_chrbuf = new char[rawlen+1];
		readsize = fread(raw_chrbuf, 1, rawlen, input);
		cnt = 0;

		#pragma omp parallel for private(j, k)
		for (j=0; j<chrlen/4; j++) {
			for (k=0; k<4; k++) {
				chrbuf[j*4+k] = bit_to_seq((raw_chrbuf[j]>>((3-k)*2))&0x3);
			}
		}

		if (rem) {
			for (j=0; j<rem; j++)
				chrbuf[chrlen-rem+j] = bit_to_seq((raw_chrbuf[rawlen-1]>>((3-j)*2))&0x3);
		}

		for (j=0; j<nblockcnt; j++) {
			memset(chrbuf + nblockstarts[j], 'N', nblocksizes[j]);
			//for (k=nblockstarts[j]; k<nblockstarts[j]+nblocksizes[j]; k++) {
			//	chrbuf[k] = 'N';
			//}
		}
		content += chrbuf;
		delete [] chrbuf;
		nblockstarts.clear();
		nblocksizes.clear();
		if (i < chrcnt-1) {
			content += ';';
			chrpos.push_back(content.size());
		}
	}
	fclose(input);
	return 0;
}
