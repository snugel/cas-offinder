#pragma OPENCL EXTENSION cl_khr_global_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_base_atomics : enable
#pragma OPENCL EXTENSION cl_khr_global_int32_extended_atomics : enable
#pragma OPENCL EXTENSION cl_khr_local_int32_extended_atomics : enable


#define uint64_t unsigned long 
#define uint32_t unsigned int
#define size_t unsigned long
#define uint16_t unsigned short
#define block_ty uint32_t
#define LOCAL_BLOCK_SIZE 4
#define blocks_avail (sizeof(block_ty) * 2)

struct s_match{
    uint64_t loc;
    uint32_t mismatches;
    uint32_t pattern_idx;
};
typedef struct s_match match;

block_ty popcount_bit_twiddle64(block_ty n){
    n = (n & 0x5555555555555555ul) + ((n >> 1) & 0x5555555555555555ul);
    n = (n & 0x3333333333333333ul) + ((n >> 2) & 0x3333333333333333ul);
    n = (n & 0x0f0f0f0f0f0f0f0ful) + ((n >> 4) & 0x0f0f0f0f0f0f0f0ful);
    n = (n & 0x00ff00ff00ff00fful) + ((n >> 8) & 0x00ff00ff00ff00fful);
    n = (n & 0x0000ffff0000fffful) + ((n >>16) & 0x0000ffff0000fffful);
    n = (n & 0x00000000fffffffful) + ((n >>32) & 0x00000000fffffffful);
    return n;
}

block_ty popcount_bit_twiddle32(block_ty n){
    n = (n & 0x55555555u) + ((n >> 1) & 0x55555555u);
    n = (n & 0x33333333u) + ((n >> 2) & 0x33333333u);
    n = (n & 0x0f0f0f0fu) + ((n >> 4) & 0x0f0f0f0fu);
    n = (n & 0x00ff00ffu) + ((n >> 8) & 0x00ff00ffu);
    n = (n & 0x0000ffffu) + ((n >>16) & 0x0000ffffu);
    return n;
}

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
            block_ty prev = genome[genome_idx + l];
            block_ty next = genome[genome_idx + l + 1];
            block_ty cur = k == 0 ? prev : (prev << (k * 4)) | (next >> ((blocks_avail - k) * 4));
            for (size_t x = 0; x < local_pattern_size; x++)
            {
                pattern_counts[x] += popcount(cur & pattern_blocks[(pattern_block_idx + x)*blocks_per_pattern + l]);
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
