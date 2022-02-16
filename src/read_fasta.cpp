#include "bit4ops.h"
#include "ceildiv.h"
#include "read_genome.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>

constexpr size_t MAX_LINE_LEN = 1 << 10;
struct FastaReader
{
    FILE* file;
    char linedata[MAX_LINE_LEN];
};

FastaReader* create_fasta_reader(const char* path)
{
    FILE* file = fopen(path, "r");
    if (!file) {
        return nullptr;
    }
    FastaReader* reader = new FastaReader{
        .file = file,
    };
    // file needs to start with a chromosome name
    fgets(reader->linedata, MAX_LINE_LEN, reader->file);
    if (reader->linedata[0] != '>') {
        fclose(file);
        return nullptr;
    }
    return reader;
}
ChromData read_next_fasta(FastaReader* reader)
{
    // if file is empty:
    if (feof(reader->file)) {
        return ChromData{
            .name = nullptr,
            .bit4data = nullptr,
            .n_nucl = 0,
        };
    }
    char* linedata = reader->linedata;
    FILE* file = reader->file;
    int namelen = strlen(linedata + 1);
    char* name = (char*)malloc(namelen + 1);
    memcpy(name, linedata + 1, namelen + 1);

    int data_size = 1 << 16;
    uint8_t* outdata = (uint8_t*)malloc(data_size);
    int chromsize = 0;
    while (fgets(linedata, MAX_LINE_LEN, reader->file) && linedata[0] != '>') {
        int linelen = strlen(linedata);
        for (; linelen > 0 &&
               (linedata[linelen - 1] == '\r' || linedata[linelen - 1] == '\n');
             linelen--)
            ;
        if (linelen + chromsize > data_size * 2) {
            data_size *= 2;
            outdata = (uint8_t*)realloc(outdata, data_size);
        }
        for (size_t i = 0; i < linelen; i++) {
            linedata[i] = to_upper(linedata[i]);
        }
        str2bit4(outdata, linedata, chromsize, linelen);
        chromsize += linelen;
    }
    return ChromData{
        .name = name,
        .bit4data = outdata,
        .n_nucl = chromsize,
    };
}
void free_fasta_reader(FastaReader** rptr)
{
    fclose((*rptr)->file);
    delete *rptr;
    *rptr = nullptr;
}
