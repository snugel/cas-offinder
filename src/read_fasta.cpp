#include "bit4ops.h"
#include "ceildiv.h"
#include "read_genome.h"
#include <cstring>
#include <fstream>

struct FastaReader
{
    char* data;
    size_t size;
    size_t idx;
};

FastaReader* create_fasta_reader(const char* path)
{
    std::ifstream is(path);
    if (!is) {
        return nullptr;
    }
    // Determine the file length
    is.seekg(0, std::ios_base::end);
    std::size_t size = is.tellg();
    is.seekg(0, std::ios_base::beg);
    // Create a vector to store the data
    char* data = new char[size];
    // Load the data
    is.read(data, size);
    if (data[0] != '>') {
        delete[] data;
        return nullptr;
    }
    return new FastaReader{
        .data = data,
        .size = size,
        .idx = 0,
    };
}
ChromData read_next_fasta(FastaReader* reader)
{
    char chrname[1 << 10] = { 0 };
    reader->idx++; // skip >
    if (reader->idx >= reader->size) {
        return ChromData{
            .name = nullptr,
            .bit4data = nullptr,
            .n_nucl = 0,
        };
    }
    size_t namepos = 0;
    for (; reader->idx < reader->size && reader->data[reader->idx] != '\n';
         reader->idx++, namepos++) {
        chrname[namepos] = reader->data[reader->idx];
    }
    size_t widx = 0;
    for (; reader->idx < reader->size && reader->data[reader->idx] != '>';
         reader->idx++) {
        for (;
             reader->idx < reader->size && reader->data[reader->idx] != '\n' &&
             reader->data[reader->idx] != '\r';
             reader->idx++, widx++) {
            reader->data[widx] = to_upper(reader->data[reader->idx]);
        }
    }
    uint8_t* outdata = new uint8_t[cdiv(widx, 2)]();
    str2bit4(outdata, reader->data, 0, widx);
    char* name = (char*)malloc(namepos + 1);
    memcpy(name, chrname, namepos + 1);
    return ChromData{
        .name = name,
        .bit4data = outdata,
        .n_nucl = widx,
    };
}
void free_fasta_reader(FastaReader** rptr)
{
    delete[](*rptr)->data;
    delete *rptr;
    *rptr = nullptr;
}
