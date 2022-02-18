#include "RangeIterator.h"
#include "async_search.h"
#include "bit4ops.h"
#include "blockify.h"
#include "bulge_logic.h"
#include "ceildiv.h"
#include "parse_input.h"
#include "search.h"
#include "bidirect_search.h"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace std;

static DeviceType get_dev_ty(char c)
{
    switch (c) {
        case 'C': return CPU;
        case 'G': return GPU;
        case 'A': return ACCEL;
    }
    throw runtime_error(
      "not a valid device name, must be one of 'G', 'C', 'A'");
}
static ostream* out_str_ptr;
static InFileInfo input;
static string all_compares;
static string pattern_str;

static void async_callback(const GenomeMatch* gm)
{
    if (!fits_pattern(gm, pattern_str.data(), input.pattern_size)) {
        return;
    }
    char dir = get_dir(gm);
    auto pair = get_dna_rna_match(all_compares.data(),input.pattern_size, gm);
    std::string dna = pair.first;
    std::string rna = pair.second;
    (*out_str_ptr) << rna << '\t' << gm->chrom_name << '\t' << gm->chrom_loc
                   << '\t' << dna << '\t' << dir << '\t' << gm->mismatches;
    if (input.ids) {
        const char* idstr = input.ids[gm->pattern_idx / 2];
        (*out_str_ptr) << '\t' << idstr;
    }
    (*out_str_ptr) << '\n';
}

int main(int argc, char** argv)
{
    if (argc != 4) {
        cerr << "requires 3 CLI arguments, in_file, device, out_file\n";
        return 1;
    }
    const char* in_fname = argv[1];
    if (strlen(argv[2]) != 1) {
        cerr << "device should be 'C','G', or 'A'\n";
        return 1;
    }
    char device_chr = argv[2][0];
    const char* out_fname = argv[3];

    auto start = std::chrono::system_clock::now();

    input = read_file(in_fname);
    DeviceType device_ty = get_dev_ty(device_chr);
    ofstream out_file;

    all_compares = mirror_pattern(input.compares,input.num_patterns, input.pattern_size);
    pattern_str = mirror_pattern(input.pattern, 1, input.pattern_size);
    if (strlen(out_fname) == 1 && out_fname[0] == '-') {
        out_str_ptr = &cout;
    } else {
        out_file.open(out_fname);
        out_str_ptr = &out_file;
    }
    //    *out_str_ptr <<
    //    "RNA\tChromosome\tLocation\tDNA\tDirection\tMismatches\tID\n";

    async_search(input.genome_path,
                 device_ty,
                 all_compares.data(),
                 input.pattern_size,
                 input.num_patterns * 2,
                 input.mismatches,
                 nullptr,
                 async_callback);

    auto end = std::chrono::system_clock::now();
    auto duration = (end - start);
    auto micros =
      std::chrono::duration_cast<std::chrono::microseconds>(duration);
    cerr << "successfully finished in: " << micros.count() / 1000000.
         << " seconds\n";
    return 0;
}
