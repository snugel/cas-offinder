#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable


#define uint64_t unsigned long 
#define uint32_t unsigned int
#define size_t unsigned long
#define uint16_t unsigned short
#define block_ty uint32_t
#define blocks_avail (sizeof(block_ty) * 2)

struct s_match{
    uint64_t loc;
    uint32_t mismatches;
    uint32_t pattern_idx;
};
typedef struct s_match match;

__kernel void find_matches_packed_helper(
    __global block_ty *genome,
    __global block_ty *pattern_blocks,
    size_t num_patterns,
    size_t blocks_per_pattern,
    int max_mismatches,
    int max_matches,
    __global match * match_buffer,
    __global int * entrycount
){
    size_t genome_idx = get_global_id(0);
    size_t pattern_block_idx = get_global_id(1);
    // genome is expected to be at least BIGGER than the genome_size
    for (size_t k = 0; k < blocks_avail; k++)
    {
        uint32_t count = 0;
        for (size_t l = 0; l < blocks_per_pattern; l++)
        {
            block_ty prev = genome[genome_idx + l];
            block_ty next = genome[genome_idx + l + 1];
            block_ty cur = k == 0 ? prev : (prev << (k * 4)) | (next >> ((blocks_avail - k) * 4));
            count += popcount(cur & pattern_blocks[pattern_block_idx*blocks_per_pattern + l]);
        }
        int mismatches = max_matches - count;
        if (mismatches <= max_mismatches)
        {   
            int next_idx = atomic_inc(entrycount);
            match next_item = {
                .loc = genome_idx * blocks_avail + k,
                .mismatches = mismatches,
                .pattern_idx = pattern_block_idx,
            };
            match_buffer[next_idx] = next_item;
        }
    }
}
