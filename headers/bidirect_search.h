#pragma once
#include <string>
#include "async_search.h"


std::string mirror_pattern(const char * compares, size_t num_compares, size_t pattern_size);
bool matches(const char* dna, const char* rna, size_t size);
bool fits_pattern(const GenomeMatch* gm, const char * doubled_pattern, size_t pattern_size);
std::pair<std::string, std::string> get_dna_rna_match(const char * compares, size_t pattern_size, const GenomeMatch* gm);
char get_dir(const GenomeMatch* gm);
