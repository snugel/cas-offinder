#include "blockify.h"
#include "RangeIterator.h"
#include "test/test_framework.h"
#include <cstring>
#include <vector>

using namespace std;
static vector<Block> convert(size_t block_size,
                             size_t block_padding,
                             size_t metadata_size,
                             vector<Chunk> chunks)
{
    vector<Block> res;
    Blockifier* gen =
      create_blockifier(block_size, block_padding, metadata_size);
    for (Chunk c : chunks) {
        add_chunk(gen, &c);

        while (is_block_ready(gen)) {
            res.push_back(get_block(gen));
        }
    }
    return res;
}

TEST(test_blockify_small_chunks)
{
    int metadata = 42;
    uint8_t data[] = { 0x41, 0x21, 0x01 };
    vector<Chunk> small_chunks = {
        Chunk{ .data = data, .size = 3, .metadata = &metadata },
        Chunk{ .data = data, .size = 2, .metadata = &metadata },
        Chunk{ .data = data, .size = 3, .metadata = &metadata },
    };
    vector<Block> actual_l = convert(8, 2, sizeof(metadata), small_chunks);
    t_assert(actual_l.size() == 1);
    Block actual = actual_l.at(0);
    uint8_t expected[] = { 0x41, 0x21, 0x01, 0x41, 0x21, 0x41, 0x21, 0x01 };
    size_t expected_poses[] = { 0, 3, 5 };
    int expected_meatadata[] = { 42, 42, 42 };
    t_check(actual.end == 8);
    t_check(!memcmp(actual.buf, expected, sizeof(expected)));
    t_check(actual.start_chunk_loc == 0);
    t_check(actual.num_chunks == 3);
    t_check(           !memcmp(
                           actual.chunk_poses, expected_poses, sizeof(expected_poses)));
    t_check(           !memcmp(
                           actual.metadatas, expected_meatadata, sizeof(expected_meatadata)));
    return true;
}

TEST(test_blockify_small_blocks)
{
    int metadata1 = 1;
    int metadata2 = 2;
    uint8_t data[] = { 0x41, 0x21, 0x01, 0x89, 0x72, 0x49, 0x29 };
    vector<Chunk> chunks = {
        Chunk{ .data = data, .size = 6, .metadata = &metadata1 },
        Chunk{ .data = data, .size = 7, .metadata = &metadata2 },
    };
    vector<Block> actual_l = convert(5, 3, sizeof(metadata1), chunks);
    uint8_t expected_c[][5]{
        { 0x41, 0x21, 0x01, 0x89, 0x72 },
        { 0x01, 0x89, 0x72, 0x49, 0x41 },
        { 0x41, 0x21, 0x01, 0x89, 0x72 },
        { 0x01, 0x89, 0x72, 0x49, 0x29 },
    };
    t_assert(actual_l.size() == 4);
    for (size_t i : range(4)) {
        t_check(!memcmp(expected_c[i], actual_l[i].buf, sizeof(expected_c[0])));
    }
    size_t expected_poses_2[] = { 0, 1 };
    int expected_meatadata_2[] = { 1, 2 };
    t_check(!memcmp(expected_poses_2, actual_l[2].chunk_poses, sizeof(expected_poses_2)));
    t_check(!memcmp(expected_meatadata_2, actual_l[2].metadatas, sizeof(expected_meatadata_2)));
    return true;
}
