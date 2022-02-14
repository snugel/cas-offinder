#include "blockify.h"
#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <vector>

using namespace std;

struct Blockifier
{
    Block b;
    Chunk c;
    size_t cloc;
    size_t b_size;
    size_t b_padding;
    size_t metadata_size;
    size_t meta_buf_cap;
};
constexpr size_t DEF_BUF_SIZE = 3;
static Block new_block(Blockifier* gen)
{
    gen->meta_buf_cap = DEF_BUF_SIZE;
    return Block{
        .buf = (uint8_t*)(malloc(gen->b_size + gen->b_padding)),
        .end = 0,
        .start_chunk_loc = 0,
        .num_chunks = 0,
        .chunk_poses = (size_t*)malloc(sizeof(size_t) * DEF_BUF_SIZE),
        .metadatas = malloc(gen->metadata_size * DEF_BUF_SIZE),
        .metadata_bytes = gen->metadata_size,
        .buf_size = gen->b_size,
    };
}
Blockifier* create_blockifier(size_t block_size,
                              size_t block_padding,
                              size_t metadata_size)
{
    Blockifier* gen = new Blockifier{
        .b = Block{},
        .c = Chunk{ .data = nullptr, .size = 0, .metadata = nullptr },
        .cloc = 0,
        .b_size = block_size,
        .b_padding = block_padding,
        .metadata_size = metadata_size,
        .meta_buf_cap = DEF_BUF_SIZE,
    };
    gen->b = new_block(gen);
    return gen;
}
bool is_block_ready(Blockifier* gen)
{
    return gen->b.end >= gen->b_size;
}
static void add_metadata(Blockifier* gen)
{
    size_t new_idx = gen->b.num_chunks;
    // new_idx+1 is to leave some room for padding
    if (new_idx + 1 >= gen->meta_buf_cap) {
        gen->meta_buf_cap *= 2;
        gen->b.chunk_poses = (size_t*)realloc(
          gen->b.chunk_poses, gen->meta_buf_cap * sizeof(size_t));
        gen->b.metadatas = (size_t*)realloc(
          gen->b.metadatas, gen->meta_buf_cap * gen->metadata_size);
    }
    gen->b.chunk_poses[new_idx] = gen->b.end;
    memcpy((char*)(gen->b.metadatas) + new_idx * gen->metadata_size,
           gen->c.metadata,
           gen->metadata_size);
    gen->b.num_chunks++;
}
static void populate_block(Blockifier* gen)
{
    size_t padded_write_len =
      min(gen->b_size + gen->b_padding - gen->b.end, gen->c.size - gen->cloc);
    size_t write_len = min(gen->b_size - gen->b.end, gen->c.size - gen->cloc);
    assert(int(write_len) >= 0 && int(padded_write_len) >= 0);
    if (padded_write_len) {
        add_metadata(gen);
        // if there was a previous chunk and a new block
        if (gen->cloc != 0 && gen->cloc < gen->c.size && gen->b.end == 0) {
            gen->b.start_chunk_loc = gen->cloc;
        }
        const uint8_t* start_src = gen->c.data + gen->cloc;
        uint8_t* start_dest = gen->b.buf + gen->b.end;
        memcpy(start_dest, start_src, padded_write_len);
        gen->b.end += padded_write_len;
        gen->cloc += write_len;
    }
}
void add_chunk(Blockifier* gen, const Chunk* chunk)
{
    assert(!is_block_ready(gen) && "cannot call `add_chunk` if a block is "
                                   "ready, must call `pop_block` first");
    gen->c = *chunk;
    gen->cloc = 0;
    populate_block(gen);
}

Block pop_block(Blockifier* gen)
{
    Block next_b = gen->b;
    // fill remaining unfilled padding with zeros
    if (int(gen->b_size + gen->b_padding - next_b.end) > 0)
        memset(next_b.buf + next_b.end,
               0,
               gen->b_size + gen->b_padding - next_b.end);
    next_b.chunk_poses[next_b.num_chunks] = next_b.end;
    gen->b = new_block(gen);
    populate_block(gen);
    return next_b;
}
void free_blockifier(Blockifier**) {}

// bool in_chunk(const Blockifier* gen, const Block* b, size_t bidx)
//{
//     return b->end > bidx && bidx < 0;
// }
ChunkInfo get_chunk_info(const Block* b, size_t bidx)
{
    if (bidx >= b->end || bidx >= b->buf_size) {
        return ChunkInfo{ .dataidx = size_t(-1), .metadata = 0 };
    }
    size_t midx =
      (std::upper_bound(b->chunk_poses, b->chunk_poses + b->num_chunks, bidx) -
       b->chunk_poses) -
      1;
    size_t dataidx =
      midx == 0 ? b->start_chunk_loc + bidx : bidx - b->chunk_poses[midx];
    void* metadata = (char*)(b->metadatas) + (midx * b->metadata_bytes);

    return ChunkInfo{ .dataidx = dataidx, .metadata = metadata };
}
