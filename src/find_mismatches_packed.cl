#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable

#define uint64_t unsigned long
#define uint32_t unsigned int
#define size_t unsigned long
#define uint16_t unsigned short
#define blocks_avail (sizeof(block_ty) * 2)
#define blocks_per_pattern ((pattern_size - 1 + blocks_avail) / blocks_avail)

struct s_match
{
    uint64_t loc;
    uint32_t pattern_idx;
    uint32_t mismatches;
};
typedef struct s_match match;

__kernel void find_matches(__global block_ty* genome,
                           __global block_ty* pattern_blocks,
                           int max_mismatches,
                           int max_matches,
                           __global match* match_buffer,
                           __global int* entrycount)
{
    size_t genome_idx = get_global_id(0);
    size_t pattern_block_idx = get_global_id(1);
    block_ty shifted_blocks[blocks_per_pattern + 1];
    for (size_t i = 0; i < blocks_per_pattern + 1; i++) {
        shifted_blocks[i] = genome[genome_idx + i];
    }
    // genome is expected to be at least BIGGER than the genome_size
    for (size_t k = 0; k < blocks_avail; k++) {
        uint32_t count = 0;
        for (size_t l = 0; l < blocks_per_pattern; l++) {
            block_ty cur = shifted_blocks[l];
            count += popcount(
              cur & pattern_blocks[pattern_block_idx * blocks_per_pattern + l]);
        }
        for (size_t l = 0; l < blocks_per_pattern; l++) {
            shifted_blocks[l] =
              (shifted_blocks[l] >> 4) |
              (shifted_blocks[l + 1] << ((blocks_avail - 1) * 4));
        }
        shifted_blocks[blocks_per_pattern] >>= 4;
        int mismatches = max_matches - count;
        if (mismatches <= max_mismatches) {
            int next_idx = atomic_inc(entrycount);
            match next_item = {
                .loc = genome_idx * blocks_avail + k,
                .pattern_idx = pattern_block_idx,
                .mismatches = mismatches,
            };
            match_buffer[next_idx] = next_item;
        }
    }
}
