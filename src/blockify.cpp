#include "blockify.h"
#include <cstdlib>
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
Blockifier* create_blockifier(size_t block_size,
                              size_t block_padding,
                              size_t metadata_size)
{
    const size_t DEF_BUF_SIZE = 2;
    return new Blockifier{
        .b =
          Block{
            .buf = (uint8_t*)(malloc(block_size + block_padding)),
            .end = 0,
            .start_chunk_loc = 0,
            .num_chunks = 0,
            .chunk_poses = (size_t*)malloc(sizeof(size_t) * DEF_BUF_SIZE),
            .metadatas = malloc(metadata_size * DEF_BUF_SIZE),
          },
        .c = Chunk{ .data = nullptr, .size = 0, .metadata = nullptr },
        .cloc = 0,
        .b_size = block_size,
        .b_padding = block_padding,
        .metadata_size = metadata_size,
        .meta_buf_cap = DEF_BUF_SIZE,
    };
}
bool is_block_ready(Blockifier* gen)
{
    return false;
}
void add_chunk(Blockifier* gen, Chunk* chunk) {}
Block get_block(Blockifier* gen)
{
    return Block{};
}
void free_blockifier(Blockifier**) {}
