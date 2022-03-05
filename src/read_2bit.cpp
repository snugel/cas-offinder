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
void sfread(char* ptr, size_t dsize, size_t dcount, FILE* file)
{
    size_t readlen = fread(ptr, dsize, dcount, file);
    if (readlen != dcount)
        throw runtime_error("invalid 2bit file");
}
static unsigned int read_uint(FILE* input)
{
    unsigned int num = 0;
    sfread((char*)&num, 4, 1, input);
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

    char len_chrname;
    char chrname[256];
    chrcnt = read_uint(input);
    fseek(input, 4, SEEK_CUR); // Reserved

    for (i = 0; i < chrcnt; i++) {
        sfread(&len_chrname, 1, 1, input);
        sfread(chrname, 1, len_chrname, input);
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

ChromData read_next_2bit(TwoBitReader* reader)
{
    if (reader->chridx >= reader->chrcnt) {
        return ChromData{};
    }
    FILE* input = reader->file;

    unsigned int j, chrlen, nblockcnt, maskblockcnt, rawlen;

    chrlen = read_uint(input);
    nblockcnt = read_uint(input);

    vector<unsigned int> nblockstarts(nblockcnt);
    vector<unsigned int> nblocksizes(nblockcnt);
    sfread((char*)&nblockstarts[0], sizeof(nblockstarts[0]), nblockcnt, input);
    sfread((char*)&nblocksizes[0], sizeof(nblocksizes[0]), nblockcnt, input);

    maskblockcnt = read_uint(input);
    fseek(input, maskblockcnt * 8 + 4, SEEK_CUR);

    rawlen = cdiv(chrlen, 4);
    vector<uint8_t> raw_buf(rawlen);
    sfread((char*)raw_buf.data(), 1, rawlen, input);

    uint8_t* data = (uint8_t*)malloc(roundup(chrlen, 4) / 2);
    twobit2bit4(data, raw_buf.data(), roundup(chrlen, 4));
    memsetbit4(data, 0, chrlen, roundup(chrlen, 4));

    for (j = 0; j < nblockcnt; j++) {
        memsetbit4(data, 0, nblockstarts[j], nblockstarts[j] + nblocksizes[j]);
    }
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
