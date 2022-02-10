#include "blockify.h"
#include "RangeIterator.h"
#include "test/test_framework.h"
#include <cstring>
#include <iostream>
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
            res.push_back(pop_block(gen));
        }
    }
    res.push_back(pop_block(gen));
    if (res.back().end == 0) {
        res.pop_back();
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
    t_check(
      !memcmp(actual.chunk_poses, expected_poses, sizeof(expected_poses)));
    t_check(!memcmp(
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
    vector<Block> actual_l = convert(2, 3, sizeof(metadata1), chunks);
    uint8_t expected_c[][5]{
        { 0x41, 0x21, 0x01, 0x89, 0x72 },
        { 0x01, 0x89, 0x72, 0x49, 0 },
        { 0x72, 0x49, 0,0,0},
        { 0x41, 0x21, 0x01, 0x89, 0x72 },
        { 0x01, 0x89, 0x72, 0x49, 0x29 },
        { 0x72, 0x49, 0x29, 0, 0 },
        { 0x29, 0, 0 , 0, 0},
    };
    constexpr size_t expected_blocks = sizeof(expected_c)/sizeof(expected_c[0]);
    t_assert(actual_l.size() == expected_blocks);
    for (size_t i : range(expected_blocks)) {
        //        for (size_t j : range(5)) {
        //            cerr << hex << int(expected_c[i][j]) << dec << "\t";
        //        }
        //        cerr << "\n";
        //        for (size_t j : range(5)) {
        //            cerr << hex << int(actual_l[i].buf[j]) << dec << "\t";
        //        }
        //        cerr << "\n" << i << "\n";
        t_check(!memcmp(expected_c[i], actual_l[i].buf, sizeof(expected_c[0])));
    }
    size_t expected_poses_1[] = { 0 };
    int expected_meatadata_1[] = { 1 };
    t_check(actual_l[0].start_chunk_loc == 0);
    t_check(actual_l[1].start_chunk_loc == 2);
    t_check(actual_l[2].start_chunk_loc == 4);
    t_check(actual_l[3].start_chunk_loc == 0);
    t_check(!memcmp(
      expected_poses_1, actual_l[1].chunk_poses, sizeof(expected_poses_1)));
    t_check(!memcmp(expected_meatadata_1,
                    actual_l[1].metadatas,
                    sizeof(expected_meatadata_1)));
    return true;
}

TEST(get_chunk_info){
    uint8_t data[] = { 0x41, 0x21, 0x01, 0x89, 0x72, 0x49, 0x29 };
    int metadata1 = 42;
    int metadata2 = 91;
    Chunk c1{ .data = data, .size = 3, .metadata = &metadata1 };
    Chunk c2{ .data = data+1, .size = 6, .metadata = &metadata2 };

    Blockifier* gen =
      create_blockifier(4, 3, sizeof(metadata1));

    add_chunk(gen, &c1);
    add_chunk(gen, &c2);
    {
    Block b = pop_block(gen);
    t_check(*(int*)(get_chunk_info(gen, &b, 2).metadata) == metadata1);
    t_check(*(int*)(get_chunk_info(gen, &b, 3).metadata) == metadata2);
    t_check(get_chunk_info(gen, &b, 4).metadata == nullptr);
    t_check(get_chunk_info(gen, &b, 2).dataidx == 2);
    t_check(get_chunk_info(gen, &b, 3).dataidx == 0);
    }
    {
    Block b = pop_block(gen);
    t_check(*(int*)(get_chunk_info(gen, &b, 1).metadata) == metadata2);
    t_check(get_chunk_info(gen, &b, b.end).metadata == nullptr);
    t_check(get_chunk_info(gen, &b, 1).dataidx == 2);
    }
    return true;
}
