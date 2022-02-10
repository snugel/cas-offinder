#include "async_search.h"
#include "RangeIterator.h"
#include "bit4ops.h"
#include "blockify.h"
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

void async_search(const char* genome_path,
                  DeviceType device_ty,
                  const char* compares,
                  size_t pattern_size,
                  size_t num_patterns,
                  uint32_t mismatches,
                  async_match_callback callback)
{

    FolderReader* reader = create_folder_reader(genome_path);
    if (!reader) {
        throw runtime_error(
          "Input file's first line does not contain a valid genome folder");
    }

    size_t pattern_bytes = cdiv(pattern_size, 2);
    uint8_t* patterns_bit4 =
      (uint8_t*)malloc(pattern_bytes * num_patterns);
    for (size_t i : range(num_patterns)) {
        str2bit4pattern(patterns_bit4 + i * pattern_bytes,
                        compares + pattern_size * i,
                        0,
                        pattern_size);
    }
    SearchFactory* fact = create_search_factory(device_ty);
    int num_searchers = num_searchers_avaliable(fact);
    const size_t OUT_CHUNK_SIZE = 1 << 22;
    const size_t MIN_CMPS_PER_OUT = 1 << 14;
    const size_t MAX_BLOCK_SIZE = 2 * 8 * 4;
    const size_t MIN_CHUNK_BYTES =
      roundup((OUT_CHUNK_SIZE * MIN_CMPS_PER_OUT) / (num_patterns * 2),
              MAX_BLOCK_SIZE);
    const size_t CHUNK_BYTES = min(size_t(1 << 26), MIN_CHUNK_BYTES);
    // 2 nucl per byte, 8 bytes per 64 bit word, 4 64 bit words per vector?
    const size_t CHUNK_PAD_BYTES =
      roundup(cdiv(pattern_size, 2), MAX_BLOCK_SIZE);
    const size_t PADDED_CHUNK_BYTES = CHUNK_BYTES + CHUNK_PAD_BYTES;
    Searcher* searcher1 = create_searcher(fact,
                                          0,
                                          PADDED_CHUNK_BYTES,
                                          OUT_CHUNK_SIZE,
                                          patterns_bit4,
                                          num_patterns,
                                          pattern_size);
    struct Metadata
    {
        char* name;
        size_t chromsize;
    };
    Blockifier* gen =
      create_blockifier(CHUNK_BYTES, CHUNK_PAD_BYTES, sizeof(Metadata));
    ChromData data = read_next_folder(reader);
    vector<Block> blocks;
    while (data.name) {
        Metadata m{
            .name = data.name,
            .chromsize = data.n_nucl,
        };
        Chunk c{
            .data = data.bit4data,
            .size = cdiv(data.n_nucl, 2),
            .metadata = &m,
        };
        add_chunk(gen, &c);
        while (is_block_ready(gen)) {
            blocks.push_back(pop_block(gen));
        }
        data = read_next_folder(reader);
    }
    blocks.push_back(pop_block(gen));
    if (blocks.back().end == 0) {
        blocks.pop_back();
    }
    vector<char> dna_match_buf(pattern_size+1);
    for (Block b : blocks) {
        Match* result;
        uint64_t num_result;
        search(
          searcher1, b.buf, b.end * 2, mismatches, &result, &num_result);
        for (Match m : each(result, result + num_result)) {
            size_t byte_loc = m.loc / 2;
            size_t byte_offset = m.loc % 2;
            ChunkInfo cinfo = get_chunk_info(gen, &b, byte_loc);
            Metadata metainfo = *((Metadata*)cinfo.metadata);
            size_t actual_loc = cinfo.dataidx * 2 + byte_offset;
            if (actual_loc <= metainfo.chromsize - pattern_size) {
                bit42str(dna_match_buf.data(),b.buf,actual_loc,pattern_size);
                GenomeMatch gm{
                    .dna_match=dna_match_buf.data(),
                    .chrom_name=metainfo.name,
                    .chrom_loc=actual_loc,
                    .pattern_idx=m.pattern_idx,
                    .mismatches=m.mismatches,
                };
                callback(&gm);
            }
        }
    }
}
