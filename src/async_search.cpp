#include "async_search.h"
#include "RangeIterator.h"
#include "bit4ops.h"
#include "blockify.h"
#include "ceildiv.h"
#include "parse_input.h"
#include "read_genome.h"
#include "channel.h"
#include "search.h"
#include <cassert>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

using namespace std;

struct Metadata
{
    char* name;
    size_t chromsize;
};

struct BlockOutput{
    Block b;
    Match * matches;
    size_t num_matches;
};

void reader_thread(const char* genome_path, Channel<Block> * block_stream, size_t CHUNK_BYTES, size_t CHUNK_PAD_BYTES){
    FolderReader* reader = create_folder_reader(genome_path);
    if (!reader) {
        throw runtime_error(
          "Input file's first line does not contain a valid genome folder");
    }
    ChromData data = read_next_folder(reader);
    Blockifier* gen =
      create_blockifier(CHUNK_BYTES, CHUNK_PAD_BYTES, sizeof(Metadata));

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
            block_stream->send(pop_block(gen));
        }
        data = read_next_folder(reader);
    }
    Block last_block = pop_block(gen);
    if (last_block.end > 0) {
        block_stream->send(last_block);
    }
    block_stream->terminate();
}
void searcher_thread(
        SearchFactory* fact,
        Channel<Block> * block_in_stream,
        Channel<BlockOutput> *block_out_stream,
        size_t device_idx,
        size_t OUT_CHUNK_SIZE,
        size_t PADDED_CHUNK_BYTES,
        uint8_t * patterns_bit4,
        size_t num_patterns,
        size_t pattern_size,
        size_t mismatches
        ){

    Searcher* searcher = create_searcher(fact,
                                          device_idx,
                                          OUT_CHUNK_SIZE,
                                          PADDED_CHUNK_BYTES * 2,
                                          patterns_bit4,
                                          num_patterns,
                                         pattern_size);
    Block input;
    while(block_in_stream->receive(input)){
        assert(input.end <= PADDED_CHUNK_BYTES);
        Match * matches;
        size_t num_matches;
        search(searcher, input.buf, input.end * 2, mismatches, &matches, &num_matches);
        block_out_stream->send(BlockOutput{
                                   .b=input,
                                   .matches=matches,
                                   .num_matches=num_matches,
                               });
    }
    block_out_stream->terminate();
}

void async_search(const char* genome_path,
                  DeviceType device_ty,
                  const char* compares,
                  size_t pattern_size,
                  size_t num_patterns,
                  uint32_t mismatches,
                  async_match_callback callback)
{
    const size_t OUT_CHUNK_SIZE = 1 << 22;
    const size_t MIN_CMPS_PER_OUT = 1 << 14;
    const size_t MAX_BLOCK_SIZE = 2 * 8 * 4;
    const size_t MIN_CHUNK_BYTES = roundup(
      (OUT_CHUNK_SIZE * MIN_CMPS_PER_OUT) / (num_patterns * 2), MAX_BLOCK_SIZE);
    const size_t CHUNK_BYTES = min(size_t(1 << 24), MIN_CHUNK_BYTES);
    // 2 nucl per byte, 8 bytes per 64 bit word, 4 64 bit words per vector?
    const size_t CHUNK_PAD_BYTES =
      roundup(cdiv(pattern_size, 2), MAX_BLOCK_SIZE);
    const size_t PADDED_CHUNK_BYTES = CHUNK_BYTES + CHUNK_PAD_BYTES;

    Channel<Block> block_channel(8);
    std::thread reader_thread_obj(reader_thread,genome_path, &block_channel,CHUNK_BYTES, CHUNK_PAD_BYTES);

    size_t pattern_bytes = cdiv(pattern_size, 2);
    uint8_t* patterns_bit4 = (uint8_t*)malloc(pattern_bytes * num_patterns);
    for (size_t i : range(num_patterns)) {
        str2bit4pattern(patterns_bit4 + i * pattern_bytes,
                        compares + pattern_size * i,
                        0,
                        pattern_size);
    }
    SearchFactory* fact = create_search_factory(device_ty);
    int num_searchers = num_searchers_avaliable(fact);
    vector<thread> searcher_threads;
    Channel<BlockOutput> output_stream(num_searchers*2+4, num_searchers);
    for(size_t i : range(num_searchers)){
        searcher_threads.emplace_back(
                    searcher_thread,
            fact,
            &block_channel,
                &output_stream,
                i,
                OUT_CHUNK_SIZE,
                PADDED_CHUNK_BYTES,
                patterns_bit4,
                num_patterns,
                pattern_size,
                mismatches
                );
    }
    vector<char> dna_match_buf(pattern_size + 1);
    BlockOutput output;
    while(output_stream.receive(output)){
        Match* result = output.matches;
        uint64_t num_result = output.num_matches;
        Block b = output.b;
        assert(b.end <= PADDED_CHUNK_BYTES);
        for (Match m : each(result, result + num_result)) {
            size_t byte_loc = m.loc / 2;
            size_t byte_offset = m.loc % 2;
            ChunkInfo cinfo = get_chunk_info(&b, byte_loc);
            if (cinfo.metadata) {
                Metadata metainfo = *((Metadata*)cinfo.metadata);
                size_t actual_loc = cinfo.dataidx * 2 + byte_offset;
                if (actual_loc <= metainfo.chromsize - pattern_size) {
                    bit42str(dna_match_buf.data(), b.buf, m.loc, pattern_size);
                    GenomeMatch gm{
                        .dna_match = dna_match_buf.data(),
                        .chrom_name = metainfo.name,
                        .chrom_loc = actual_loc,
                        .pattern_idx = m.pattern_idx,
                        .mismatches = m.mismatches,
                    };
                    callback(&gm);
                }
            }
        }
    }
    reader_thread_obj.join();
    for(thread & t : searcher_threads){
        t.join();
    }
}
