#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable


#define uint64_t unsigned long 
#define uint32_t unsigned int
#define size_t unsigned long
#define uint16_t unsigned short
#define LOCAL_BLOCK_SIZE 4
#define blocks_avail 16

struct s_match{
    uint64_t loc;
    uint32_t mismatches;
    uint32_t pattern_idx;
};
typedef struct s_match match;

__kernel void find_matches_packed_helper(
    __global uint64_t *genome,
    __global uint64_t *pattern_blocks,
    size_t num_patterns,
    size_t blocks_per_pattern,
    int max_mismatches,
    int max_matches,
    __global match * match_buffer,
    __global int * entrycount
){
    size_t genome_idx = get_global_id(0);
    size_t pattern_block_idx = get_global_id(1) * LOCAL_BLOCK_SIZE;
    // genome is expected to be at least BIGGER than the genome_size
    __private uint16_t pattern_counts[LOCAL_BLOCK_SIZE];
    for (size_t k = 0; k < blocks_avail; k++)
    {
        size_t local_pattern_size = min(num_patterns - pattern_block_idx, (uint64_t)LOCAL_BLOCK_SIZE);
        for (size_t x = 0; x < local_pattern_size; x++)
        {
            pattern_counts[x] = 0;
        }
        for (size_t l = 0; l < blocks_per_pattern; l++)
        {
            uint64_t prev = genome[genome_idx + l];
            uint64_t next = genome[genome_idx + l + 1];
            uint64_t cur = k == 0 ? prev : (prev << (k * 4)) | (next >> ((blocks_avail - k) * 4));
            for (size_t x = 0; x < local_pattern_size; x++)
            {
                pattern_counts[x] += __builtin_popcountl(cur & pattern_blocks[(pattern_block_idx + x)*blocks_per_pattern + l]);
            }
        }
        for (size_t x = 0; x < local_pattern_size; x++)
        {
            int mismatches = max_matches - pattern_counts[x];
            if (mismatches <= max_mismatches)
            {   
        		int next_idx = atomic_inc(entrycount);
                match next_item = {
                    .loc = genome_idx * blocks_avail + k,
                    .mismatches = mismatches,
                    .pattern_idx = pattern_block_idx + x,
                };
                match_buffer[next_idx] = next_item;
            }
        }
    }
}
