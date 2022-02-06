#include "RangeIterator.h"
#include "bit4ops.h"
#include "ceildiv.h"
#include "parse_input.h"
#include "read_genome.h"
#include "search.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

static DeviceType get_dev_ty(char c)
{
    switch (c) {
        case 'C':
            return CPU;
        case 'G':
            return GPU;
        case 'A':
            return ACCEL;
    }
    throw runtime_error(
      "not a valid device name, must be one of 'G', 'C', 'A'");
}

int main(int argc, char** argv)
{
    assert(argc == 4 && "requires 3 CLI arguments, in_file, device, out_file");
    //    assert(argc =)
    const char* in_fname = argv[1];
    assert(strlen(argv[2]) == 1 && "device should be 'C','G', or 'A'");
    char device_chr = argv[2][0];
    const char* out_fname = argv[3];

    InFileInfo input = read_file(in_fname);
    DeviceType device_ty = get_dev_ty(device_chr);
    ofstream out_file(out_fname);

    FolderReader* reader = create_folder_reader(input.genome_path);
    if (!reader) {
        throw runtime_error(
          "Input file's first line does not contain a valid genome folder");
    }

    size_t pattern_bytes = cdiv(input.pattern_size, 2);
    uint8_t* patterns_bit4 =
      (uint8_t*)malloc(pattern_bytes * input.num_patterns);
    for (size_t i : range(input.num_patterns)) {
        str2bit4(patterns_bit4 + i * pattern_bytes,
                 input.pattern + input.pattern_size * i,
                 0,
                 input.pattern_size);
    }
    SearchFactory* fact = create_search_factory(device_ty);
    int num_searchers = num_searchers_avaliable(fact);
    constexpr size_t CHUNK_SIZE = 1 << 24;
    constexpr size_t OUT_CHUNK_SIZE = 1 << 22;
    constexpr size_t MAX_BLOCK_SIZE =
      2 * 8 *
      4; // 2 nucl per byte, 8 bytes per 64 bit word, 4 64 bit words per vector?
    const size_t chunk_pad_size =
      cdiv(input.pattern_size, MAX_BLOCK_SIZE) * MAX_BLOCK_SIZE;
    const size_t PAD_CHUNK_BYTES = cdiv(CHUNK_SIZE + chunk_pad_size, 2);
    Searcher* searcher1 = create_searcher(fact,
                                          0,
                                          PAD_CHUNK_BYTES,
                                          OUT_CHUNK_SIZE,
                                          patterns_bit4,
                                          input.num_patterns,
                                          input.pattern_size);
    struct Chunk
    {
        uint8_t* buffer;
        size_t end;
        size_t startidx;
        vector<size_t> chrlocs;
        vector<char*> chrnames;
    };
    Chunk chunk{ .buffer = (uint8_t*)malloc(PAD_CHUNK_BYTES) };
    ChromData data = read_next_folder(reader);
    size_t chrloc = 0;
    while (data.name) {
        size_t cpysize = min(PAD_CHUNK_BYTES - chunk.end, data.n_nucl - chrloc);
        memcpy(data.bit4data, chunk.buffer, 3);
        data = read_next_folder(reader);
        chrloc = 0;
    }
    cerr << "successfully finished\n";
    return 0;
}
