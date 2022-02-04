#include "search.h"
#include "RangeIterator.h"
#include "ceildiv.h"
#include "oclfunctions.h"
#include "oclkernels.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

struct SearchFactory
{
    vector<cl_platform_id> platforms;
    vector<cl_device_id> devices;
    vector<cl_platform_id> dev_plat;
};
struct Searcher
{
    cl_platform_id platform;
    cl_device_id device;
    cl_context context;
    cl_program program;
    uint32_t block_size;
    uint32_t pattern_size;
    uint32_t num_patterns;
    uint64_t max_genome_size;
    uint64_t max_matches;
    cl_mem genome_buf;
    cl_mem pattern_buf;
    cl_mem match_buf;
    cl_mem count_buf;
    cl_kernel search;
    cl_command_queue queue;
};
cl_device_type to_dev_ty(DeviceType ty)
{
    switch (ty) {
        case GPU:
            return CL_DEVICE_TYPE_GPU;
        case CPU:
            return CL_DEVICE_TYPE_CPU;
        case ACCEL:
            return CL_DEVICE_TYPE_ACCELERATOR;
    }
    assert(false);
    return CL_DEVICE_TYPE_CPU;
}

SearchFactory* create_search_factory(DeviceType device_ty)
{
    constexpr size_t MAX_PLATFORM_NUM = 50;
    constexpr size_t MAX_DEVICE_NUM = 1000;

    cl_platform_id platforms[MAX_PLATFORM_NUM];
    cl_uint platform_cnt = 0;
    oclGetPlatformIDs(MAX_PLATFORM_NUM, platforms, &platform_cnt);
    if (platform_cnt == 0) {
        cerr << "No OpenCL platforms found. Check OpenCL installation!" << endl;
        exit(1);
    }
    SearchFactory* fact = new SearchFactory();
    fact->platforms.assign(platforms, platforms + platform_cnt);

    cl_device_type devtype = to_dev_ty(device_ty);
    cl_device_id devices[MAX_DEVICE_NUM];
    cl_uint device_cnt;
    for (cl_platform_id plat : fact->platforms) {
        oclGetDeviceIDs(plat, devtype, MAX_DEVICE_NUM, devices, &device_cnt);
        fact->devices.insert(fact->devices.end(), devices, devices + device_cnt);
        for (int _ : range(device_cnt)) {
            fact->dev_plat.push_back(plat);
        }
    }

    if (fact->devices.size() == 0) {
        cerr << "No OpenCL devices found." << endl;
        exit(1);
    }
    return fact;
}
int num_searchers_avaliable(SearchFactory* fact)
{
    return fact->devices.size();
}
static uint64_t num_blocks(Searcher* searcher, uint64_t size)
{
    return cdiv(size, searcher->block_size * 2);
}
static uint64_t num_bytes(Searcher* searcher, uint64_t size)
{
    return num_blocks(searcher, size) * searcher->block_size;
}
Searcher* create_searcher(SearchFactory* fact,
                          int searcher_idx,
                          uint64_t max_matches,
                          uint64_t max_genome_size,
                          uint8_t* bit4patterns,
                          uint64_t num_patterns,
                          uint64_t pattern_size)
{
    Searcher* searcher = new Searcher();
    searcher->pattern_size = pattern_size;
    searcher->num_patterns = num_patterns;

    const size_t src_len = strlen(program_src);
    searcher->device = fact->devices.at(searcher_idx);
    searcher->platform = fact->dev_plat.at(searcher_idx);
    const cl_context_properties contextProperties[] = {
        CL_CONTEXT_PLATFORM, reinterpret_cast<cl_context_properties>(searcher->platform), 0, 0
    };
    searcher->context = oclCreateContext(contextProperties, 1, &searcher->device, nullptr, nullptr);
    searcher->program = oclCreateProgramWithSource(searcher->context, 1, &program_src, &src_len);
    searcher->block_size = 4;
    string defs = "-Dpattern_size=" + to_string(pattern_size) + " -Dblock_ty=unsigned";
    oclBuildProgram(searcher->program, 1, &searcher->device, defs.c_str(), nullptr, nullptr);

    int old_pattern_size = cdiv(pattern_size, 2);
    int new_pattern_size = num_bytes(searcher, pattern_size);
    vector<uint8_t> padded_patterns(num_patterns * new_pattern_size);
    for (size_t j : range(num_patterns)) {
        for (size_t i : range(old_pattern_size)) {
            padded_patterns.at(new_pattern_size * j + i) = bit4patterns[old_pattern_size * j + i];
        }
    }
    searcher->genome_buf = oclCreateBuffer(
      searcher->context, CL_MEM_READ_WRITE, num_bytes(searcher, max_genome_size), nullptr);
    searcher->pattern_buf =
      oclCreateBuffer(searcher->context, CL_MEM_READ_WRITE, padded_patterns.size(), nullptr);
    searcher->match_buf =
      oclCreateBuffer(searcher->context, CL_MEM_READ_WRITE, max_matches * sizeof(Match), nullptr);
    searcher->count_buf =
      oclCreateBuffer(searcher->context, CL_MEM_READ_WRITE, sizeof(uint32_t), nullptr);
    searcher->search = oclCreateKernel(searcher->program, "find_matches");
    searcher->queue = oclCreateCommandQueue(searcher->context, searcher->device, 0);

    oclEnqueueWriteBuffer(searcher->queue,
                          searcher->pattern_buf,
                          CL_TRUE,
                          0,
                          padded_patterns.size(),
                          padded_patterns.data(),
                          0,
                          0,
                          0);

    oclSetKernelArg(searcher->search, 0, sizeof(cl_mem), &searcher->genome_buf);
    oclSetKernelArg(searcher->search, 1, sizeof(cl_mem), &searcher->pattern_buf);
    oclSetKernelArg(searcher->search, 3, sizeof(cl_uint), &pattern_size);
    oclSetKernelArg(searcher->search, 4, sizeof(cl_mem), &searcher->match_buf);
    oclSetKernelArg(searcher->search, 5, sizeof(cl_mem), &searcher->count_buf);
    return searcher;
}
void search(Searcher* searcher,
            uint8_t* bit4genome,
            uint64_t genome_size,
            uint32_t max_mismatches,
            Match** match_result,
            uint64_t* num_matches)
{
    oclSetKernelArg(searcher->search, 2, sizeof(cl_uint), &max_mismatches);
    oclEnqueueWriteBuffer(
      searcher->queue, searcher->genome_buf, CL_TRUE, 0, cdiv(genome_size, 2), bit4genome, 0, 0, 0);
    int zero = 0;
    oclEnqueueWriteBuffer(
      searcher->queue, searcher->count_buf, CL_TRUE, 0, sizeof(zero), &zero, 0, 0, 0);
    size_t num_pattern_blocks = num_blocks(searcher, searcher->pattern_size);
    size_t num_genome_blocks = num_blocks(searcher, genome_size) + 1 - num_pattern_blocks;
    size_t num_patterns = searcher->num_patterns;
    size_t work_sizes[] = { num_genome_blocks, num_patterns };
    oclEnqueueNDRangeKernel(searcher->queue, searcher->search, 2, 0, work_sizes, 0, 0, 0, 0);
    oclFinish(searcher->queue);
    cl_uint count;
    oclEnqueueReadBuffer(
      searcher->queue, searcher->count_buf, CL_TRUE, 0, sizeof(count), &count, 0, 0, 0);
    if (count > 0) {
        *num_matches = count;
        *match_result = (Match*)malloc(count * sizeof(Match));
        oclEnqueueReadBuffer(searcher->queue,
                             searcher->match_buf,
                             CL_TRUE,
                             0,
                             sizeof(Match) * count,
                             *match_result,
                             0,
                             0,
                             0);

    } else {
        *num_matches = 0;
        *match_result = nullptr;
    }
}
void free_searcher(Searcher**) {}
void free_search_factory(SearchFactory**) {}
