#include "bit4ops.h"
#include "ceildiv.h"
#include "read_genome.h"
#include <cassert>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

struct TwoBitReader
{
    FILE* file;
    vector<string> chrnames;
    uint64_t chrcnt;
    uint64_t chridx;
};

static unsigned int read_uint(FILE* input)
{
    unsigned int num = 0;
    size_t readlen = fread((char*)&num, 4, 1, input);
    if (readlen != 1)
        throw runtime_error("invalid 2bit file");

    return num;
}
TwoBitReader* create_2bit_reader(const char* path)
{
    FILE* input = fopen(path, "rb");

    if (!input) {
        return nullptr;
    }
    unsigned int val = read_uint(input);
    if (val != 0x1A412743) { // Magic
        fclose(input);
        return nullptr;
    }
    if (read_uint(input) != 0) { // Version should be 0
        fclose(input);
        return nullptr;
    }
    vector<string> chrnames;
    unsigned int i, chrcnt;

    size_t readsize;
    char len_chrname;
    char chrname[256];
    chrcnt = read_uint(input);
    fseek(input, 4, SEEK_CUR); // Reserved

    for (i = 0; i < chrcnt; i++) {
        readsize = fread(&len_chrname, 1, 1, input);
        readsize = fread(chrname, 1, len_chrname, input);
        chrname[int(len_chrname)] = 0;
        chrnames.push_back(string(chrname));
        fseek(input, 4, SEEK_CUR); // Absolute position of each sequence
    }

    return new TwoBitReader{
        .file = input,
        .chrnames = chrnames,
        .chrcnt = chrcnt,
        .chridx = 0,
    };
}
static char base_mapping[] = { 'T', 'C', 'A', 'G' };
static char bit_to_seq(unsigned char b)
{
    return base_mapping[b];
}

ChromData read_next_2bit(TwoBitReader* reader)
{
    if (reader->chridx >= reader->chrcnt) {
        return ChromData{};
    }
    FILE* input = reader->file;

    unsigned int j, k, chrlen, nblockcnt, maskblockcnt, rawlen, rem;
    int jj;
    vector<unsigned int> nblockstarts;
    vector<unsigned int> nblocksizes;

    size_t readsize;

    chrlen = read_uint(input);
    nblockcnt = read_uint(input);
    for (j = 0; j < nblockcnt; j++)
        nblockstarts.push_back(read_uint(input));
    for (j = 0; j < nblockcnt; j++)
        nblocksizes.push_back(read_uint(input));
    maskblockcnt = read_uint(input);
    fseek(input, maskblockcnt * 8 + 4, SEEK_CUR);

    rem = chrlen & 3;
    rawlen = chrlen / 4 + (rem == 0 ? 0 : 1);
    char* chrbuf = new char[chrlen + 1];
    chrbuf[chrlen] = 0;
    char* raw_chrbuf = new char[rawlen + 1];
    readsize = fread(raw_chrbuf, 1, rawlen, input);

    for (jj = 0; jj < (int)chrlen / 4; jj++) {
        for (k = 0; k < 4; k++) {
            chrbuf[jj * 4 + k] =
              bit_to_seq((raw_chrbuf[jj] >> ((3 - k) * 2)) & 0x3);
        }
    }

    if (rem) {
        for (j = 0; j < rem; j++)
            chrbuf[chrlen - rem + j] =
              bit_to_seq((raw_chrbuf[rawlen - 1] >> ((3 - j) * 2)) & 0x3);
    }

    for (j = 0; j < nblockcnt; j++) {
        memset(chrbuf + nblockstarts[j], 'N', nblocksizes[j]);
    }
    uint8_t* data = (uint8_t*)malloc(cdiv(chrlen, 2));
    str2bit4(data, chrbuf, 0, chrlen);
    delete[] chrbuf;
    nblockstarts.clear();
    nblocksizes.clear();

    string s = reader->chrnames.at(reader->chridx);
    char* name = (char*)malloc(s.size() + 1);
    strcpy(name, s.c_str());
    reader->chridx += 1;
    return ChromData{
        .name = name,
        .bit4data = data,
        .n_nucl = chrlen,
    };
}
void free_2bit_reader(TwoBitReader** readptr)
{
    TwoBitReader* reader = *readptr;
    if (reader->file) {
        fclose(reader->file);
    }
    delete reader;
    *readptr = nullptr;
}
