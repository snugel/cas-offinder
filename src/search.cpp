#include "search.h"
#include "RangeIterator.h"
#include "ceildiv.h"
#include "oclkernels.h"
#include "opencl_executor.h"

#include <algorithm>
#include <cassert>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <vector>

using namespace std;

struct SearchFactory
{
    std::vector<ExecutorTemplate> executor_templates;
};
struct Searcher
{
    OpenCLExecutor executor;
    uint32_t block_size;
    uint32_t pattern_size;
    uint32_t num_patterns;
    uint64_t max_genome_size;
    uint64_t max_matches;
    CLBuffer<uint8_t> genome_buf;
    CLBuffer<uint8_t> pattern_buf;
    CLBuffer<Match> match_buf;
    CLBuffer<uint32_t> count_buf;
    CLKernel search;
};
cl_device_type to_dev_ty(DeviceType ty)
{
    switch (ty) {
        case GPU: return CL_DEVICE_TYPE_GPU;
        case CPU: return CL_DEVICE_TYPE_CPU;
        case ACCEL: return CL_DEVICE_TYPE_ACCELERATOR;
    }
    throw "err";
    return CL_DEVICE_TYPE_CPU;
}

SearchFactory* create_search_factory(DeviceType device_ty)
{

    SearchFactory* fact = new SearchFactory();
    fact->executor_templates = get_executor_templates(to_dev_ty(device_ty));
    cerr << "Using devices:\n";
    for (ExecutorTemplate temp : fact->executor_templates) {
        cerr << get_template_info(temp) << "\n";
    }
    return fact;
}
int num_searchers_avaliable(SearchFactory* fact)
{
    return fact->executor_templates.size();
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
    ExecutorTemplate temp = fact->executor_templates.at(searcher_idx);
    Searcher* searcher = new Searcher();
    searcher->pattern_size = pattern_size;
    searcher->num_patterns = num_patterns;

    const size_t src_len = strlen(program_src);
    string defs = "-Dpattern_size=" + to_string(pattern_size);
    bool is_cpu = true;
    if (is_cpu) {
        searcher->block_size = 4;
        defs += " -Dblock_ty=uint32_t";
    } else {
        searcher->block_size = 8;
        defs += " -Dblock_ty=uint64_t";
    }
    searcher->executor =
      OpenCLExecutor(program_src, temp.plat, temp.device, defs);

    int old_pattern_size = cdiv(pattern_size, 2);
    int new_pattern_size = num_bytes(searcher, pattern_size);
    vector<uint8_t> padded_patterns(num_patterns * new_pattern_size);
    for (size_t j : range(num_patterns)) {
        for (size_t i : range(old_pattern_size)) {
            padded_patterns.at(new_pattern_size * j + i) =
              bit4patterns[old_pattern_size * j + i];
        }
    }

    searcher->genome_buf = searcher->executor.new_clbuffer<uint8_t>(
      num_bytes(searcher, max_genome_size) + searcher->block_size * 2);
    searcher->pattern_buf =
      searcher->executor.new_clbuffer<uint8_t>(padded_patterns.size());
    searcher->match_buf = searcher->executor.new_clbuffer<Match>(max_matches);
    searcher->count_buf = searcher->executor.new_clbuffer<uint32_t>(1);
    searcher->search = searcher->executor.new_clkernel("find_matches");

    searcher->pattern_buf.write_buffer(padded_patterns);

    return searcher;
}
void search(Searcher* searcher,
            uint8_t* bit4genome,
            uint64_t genome_size,
            uint32_t max_mismatches,
            Match** match_result,
            uint64_t* num_matches)
{
    searcher->search.set_args({
      searcher->genome_buf.k_arg(),
      searcher->pattern_buf.k_arg(),
      CLKernelArg(cl_int(max_mismatches)),
      CLKernelArg(cl_int(searcher->pattern_size)),
      searcher->match_buf.k_arg(),
      searcher->count_buf.k_arg(),
    });
    // searcher->genome_buf.clear_buffer();
    searcher->genome_buf.write_buffer(bit4genome, cdiv(genome_size, 2));
    searcher->count_buf.clear_buffer();
    size_t num_pattern_blocks = num_blocks(searcher, searcher->pattern_size);
    size_t num_genome_blocks =
      num_blocks(searcher, genome_size) + 1 - num_pattern_blocks;
    size_t num_patterns = searcher->num_patterns;
    size_t work_sizes[] = { num_genome_blocks, num_patterns };
    searcher->search.run(
      CL_NDRange(num_genome_blocks, num_patterns), CL_NDRange(), CL_NDRange());
    cl_uint count;
    searcher->count_buf.read_buffer(&count, 1);
    if (count > 0) {
        *num_matches = count;
        *match_result = (Match*)malloc(count * sizeof(Match));
        searcher->match_buf.read_buffer(*match_result, count);
    } else {
        *num_matches = 0;
        *match_result = nullptr;
    }
}
void free_searcher(Searcher**) {}
void free_search_factory(SearchFactory**) {}
