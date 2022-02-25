#pragma once
#include <stddef.h>
#include <stdint.h>

struct InFileInfo
{
    char* genome_path;
    char* pattern;
    char* compares;
    size_t pattern_size;
    size_t num_patterns;
    int dna_bulges;
    int rna_bulges;
    char** ids;
    uint32_t mismatches;
};

InFileInfo read_file(const char* in_fname);
