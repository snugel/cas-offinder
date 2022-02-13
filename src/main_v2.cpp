#include "RangeIterator.h"
#include "async_search.h"
#include "bit4ops.h"
#include "blockify.h"
#include "bulge_logic.h"
#include "ceildiv.h"
#include "parse_input.h"
#include "read_genome.h"
#include "search.h"
#include <cassert>
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
static ostream* out_str_ptr;
static InFileInfo input;
static string all_compares;
static string pattern_str;
static string rev_pattern_str;

static bool matches(const char* dna, const char* rna, size_t size)
{
    for (size_t i = 0; i < size; i++) {
        if (!is_match(dna[i], rna[i])) {
            return false;
        }
    }
    return true;
}
bool fits_pattern(const GenomeMatch* gm)
{
    return gm->pattern_idx % 2 == 0
             ? matches(gm->dna_match, pattern_str.data(), input.pattern_size)
             : matches(
                 gm->dna_match, rev_pattern_str.data(), input.pattern_size);
}
static void async_callback(const GenomeMatch* gm)
{
    if (!fits_pattern(gm)) {
        return;
    }
    std::string rna = all_compares.substr(gm->pattern_idx * input.pattern_size,
                                          input.pattern_size);
    std::string dna(gm->dna_match);
    char dir = '+';
    if (gm->pattern_idx % 2 == 1) {
        dna = reverse_compliment(dna);
        rna = reverse_compliment(rna);
        dir = '-';
    }
    indicate_mismatches_dna(dna, rna);
    (*out_str_ptr) << rna << '\t' << gm->chrom_name << '\t' << gm->chrom_loc
                   << '\t' << dna << '\t' << dir << '\t' << gm->mismatches
                   << '\n';
}

int main(int argc, char** argv)
{
    assert(argc == 4 && "requires 3 CLI arguments, in_file, device, out_file");
    //    assert(argc =)
    const char* in_fname = argv[1];
    assert(strlen(argv[2]) == 1 && "device should be 'C','G', or 'A'");
    char device_chr = argv[2][0];
    const char* out_fname = argv[3];

    input = read_file(in_fname);
    DeviceType device_ty = get_dev_ty(device_chr);
    ofstream out_file;

    for (size_t i : range(input.num_patterns)) {
        string forward_cmp(input.compares + i * input.pattern_size,
                           input.compares + (i + 1) * input.pattern_size);
        string backward_cmp = reverse_compliment(forward_cmp);
        all_compares += forward_cmp;
        all_compares += backward_cmp;
    }
    pattern_str = string(input.pattern, input.pattern + input.pattern_size);
    rev_pattern_str = reverse_compliment(pattern_str);
    if (strlen(out_fname) == 1 && out_fname[0] == '-') {
        out_str_ptr = &cout;
    } else {
        out_file.open(out_fname);
        out_str_ptr = &out_file;
    }
    //    *out_str_ptr <<
    //    "RNA\tChromosome\tLocation\tDNA\tDirection\tMismatches\n";

    async_search(input.genome_path,
                 device_ty,
                 all_compares.data(),
                 input.pattern_size,
                 input.num_patterns * 2,
                 input.mismatches,
                 async_callback);

    cerr << "successfully finished\n";
    return 0;
}
