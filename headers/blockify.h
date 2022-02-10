#pragma once
#include <stddef.h>
#include <stdint.h>

struct Block
{
    uint8_t* buf;
    size_t end;
    size_t start_chunk_loc;
    size_t num_chunks;
    size_t* chunk_poses;
    void* metadatas;
//    size_t metadata_bytes;
};
struct Chunk
{
    uint8_t* data;
    size_t size;
    void* metadata;
};
struct Blockifier;
Blockifier* create_blockifier(size_t block_size,
                              size_t block_padding,
                              size_t metadata_size);
void add_chunk(Blockifier* gen, const Chunk* chunk);
bool is_block_ready(Blockifier* gen);
Block pop_block(Blockifier* gen);
void free_blockifier(Blockifier**);

// block reading utilities
struct ChunkInfo
{
    size_t dataidx;
    void* metadata;
};
ChunkInfo get_chunk_info(const Blockifier * gen, const Block* b, size_t bidx);
