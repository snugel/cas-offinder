#include "RangeIterator.h"
#include "async_search.h"
#include "bidirect_search.h"
#include "bit4ops.h"
#include "blockify.h"
#include "bulge_logic.h"
#include "ceildiv.h"
#include "parse_input.h"
#include "search.h"
#include <chrono>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
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
struct OutputDataV2
{
    ostream* out_str_ptr;
    InFileInfo input;
    string all_compares;
    string pattern_str;
};

static void async_callback(const GenomeMatch* gm, void* user_data)
{
    OutputDataV2* data = (OutputDataV2*)(user_data);
    ostream& outs = *data->out_str_ptr;
    if (!fits_pattern(gm, data->pattern_str.data(), data->input.pattern_size)) {
        return;
    }
    const char* compares = data->all_compares.data();
    size_t pattern_size = data->input.pattern_size;
    char dir = gm->pattern_idx % 2 ? '-' : '+';
    std::string rna =
      std::string(compares + gm->pattern_idx * pattern_size,
                  compares + (gm->pattern_idx + 1) * pattern_size);
    std::string dna(gm->dna_match);

    if (gm->pattern_idx % 2 == 1) {
        i_reverse_compliment(rna.data(), rna.size());
        i_reverse_compliment(dna.data(), dna.size());
    }
    indicate_mismatches_dna(dna, rna);
    outs << rna << '\t' << gm->chrom_name << '\t' << gm->chrom_loc << '\t'
         << dna << '\t' << dir << '\t' << gm->mismatches;
    if (data->input.ids) {
        const char* idstr = data->input.ids[gm->pattern_idx / 2];
        outs << '\t' << idstr;
    }
    outs << '\n';
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

    InFileInfo input = read_file(in_fname);
    if (input.dna_bulges || input.rna_bulges) {
        cerr << "dna bulges and rna bulges must be 0 for v2\n";
        return 1;
    }
    DeviceType device_ty = get_dev_ty(device_chr);
    ofstream out_file;
    ostream* out_str_ptr;

    string all_compares =
      mirror_pattern(input.compares, input.num_patterns, input.pattern_size);
    string pattern_str = mirror_pattern(input.pattern, 1, input.pattern_size);
    if (strlen(out_fname) == 1 && out_fname[0] == '-') {
        out_str_ptr = &cout;
    } else {
        out_file.open(out_fname);
        out_str_ptr = &out_file;
    }
    //    *out_str_ptr <<
    //    "RNA\tChromosome\tLocation\tDNA\tDirection\tMismatches\tID\n";
    OutputDataV2 out_data{
        .out_str_ptr = out_str_ptr,
        .input = input,
        .all_compares = all_compares,
        .pattern_str = pattern_str,
    };
    async_search(input.genome_path,
                 device_ty,
                 all_compares.data(),
                 input.pattern_size,
                 input.num_patterns * 2,
                 input.mismatches,
                 nullptr,
                 &out_data,
                 async_callback);

    auto end = std::chrono::system_clock::now();
    auto duration = (end - start);
    auto micros =
      std::chrono::duration_cast<std::chrono::microseconds>(duration);
    cerr << "successfully finished in: " << micros.count() / 1000000.
         << " seconds\n";
    return 0;
}
