#include "bidirect_search.h"
#include "bit4ops.h"
#include "bulge_logic.h"
#include "RangeIterator.h"

using namespace std;

std::string mirror_pattern(const char * compares, size_t num_compares, size_t pattern_size){
    std::string all_compares;
    for (size_t i : range(num_compares)) {
        string forward_cmp(compares + i * pattern_size,
                           compares + (i + 1) * pattern_size);
        string backward_cmp = reverse_compliment(forward_cmp);
        all_compares += forward_cmp;
        all_compares += backward_cmp;
    }
    return all_compares;
}
static bool matches(const char* dna, const char* rna, size_t size)
{
    for (size_t i = 0; i < size; i++) {
        if (!is_match(dna[i], rna[i])) {
            return false;
        }
    }
    return true;
}
bool fits_pattern(const GenomeMatch* gm, const char * doubled_pattern, size_t pattern_size){
    return gm->pattern_idx % 2 == 0
             ? matches(gm->dna_match, doubled_pattern, pattern_size)
             : matches(
                 gm->dna_match, doubled_pattern + pattern_size, pattern_size);
}
std::pair<std::string, std::string> get_dna_rna_match(const char * compares, size_t pattern_size, const GenomeMatch* gm){
    std::string rna = std::string(compares + gm->pattern_idx * pattern_size,
                                  compares + (gm->pattern_idx+1) * pattern_size);
    std::string dna(gm->dna_match);

    if (gm->pattern_idx % 2 == 1) {
        rna = reverse_compliment(rna);
        dna = reverse_compliment(dna);
    }
    indicate_mismatches_dna(dna, rna);
    return std::make_pair(dna, rna);
}

char get_dir(const GenomeMatch* gm){
    return gm->pattern_idx % 2 ? '-' : '+';
}
