#include "bulge_logic.h"
#include "RangeIterator.h"
#include "bit4ops.h"
#include "test/test_framework.h"
#include <algorithm>
#include <cassert>
#include <iostream>

std::string get_bulge_type_name(BulgeType type)
{
    switch (type) {
        case BULGE_DNA: return "DNA";
        case BULGE_RNA: return "RNA";
        case BULGE_NONE: return "X";
        default: throw std::runtime_error("bad bulge type");
    }
}

int num_leading_ns(bulge_augment& augment, int dna_bulges)
{
    switch (augment.bulge_type) {
        case BULGE_DNA: return dna_bulges - augment.bulge_size;
        case BULGE_RNA: return dna_bulges + augment.bulge_size;
        case BULGE_NONE: return dna_bulges;
    }
    throw std::runtime_error("invalid bulge type");
}
std::string get_rna_match(std::string base_rna_match, bulge_augment augment)
{
    switch (augment.bulge_type) {
        case BULGE_DNA:
            return base_rna_match.substr(0, augment.bulge_pos) +
                   std::string(augment.bulge_size, '-') +
                   base_rna_match.substr(augment.bulge_pos);
        case BULGE_RNA: return base_rna_match;
        case BULGE_NONE: return base_rna_match;
    }
    throw std::runtime_error("invalid bulge type");
}

std::string get_dna_match(std::string base_dna_match, bulge_augment augment)
{
    switch (augment.bulge_type) {
        case BULGE_DNA: return base_dna_match;
        case BULGE_RNA:
            return base_dna_match.substr(0, augment.bulge_pos) +
                   std::string(augment.bulge_size, '-') +
                   base_dna_match.substr(augment.bulge_pos);
        case BULGE_NONE: return base_dna_match;
    }
    throw std::runtime_error("invalid bulge type");
}

void indicate_mismatches_dna(std::string& dna_match, std::string& rna_match)
{
    assert(dna_match.size() == rna_match.size());
    for (size_t i : range(dna_match.size())) {
        if (dna_match[i] != '-' && rna_match[i] != '-' &&
            !is_match(dna_match[i], rna_match[i])) {
            dna_match[i] |= 0x20;
        }
    }
}
std::string reverse_compliment(std::string seq)
{
    for (size_t i : range(seq.size())) {
        if (seq[i] == 'A')
            seq[i] = 'T';
        else if (seq[i] == 'T')
            seq[i] = 'A';
        else if (seq[i] == 'G')
            seq[i] = 'C';
        else if (seq[i] == 'C')
            seq[i] = 'G';
        else if (seq[i] == 'R')
            seq[i] = 'Y';
        else if (seq[i] == 'Y')
            seq[i] = 'R';
        else if (seq[i] == 'M')
            seq[i] = 'K';
        else if (seq[i] == 'K')
            seq[i] = 'M';
        else if (seq[i] == 'H')
            seq[i] = 'D';
        else if (seq[i] == 'D')
            seq[i] = 'H';
        else if (seq[i] == 'B')
            seq[i] = 'V';
        else if (seq[i] == 'V')
            seq[i] = 'B';
    }
    std::reverse(seq.begin(), seq.end());
    return seq;
}

bulge_info get_bulge_info(std::string base_dna_match,
                          std::string base_rna_match,
                          bulge_augment augment,
                          int orig_loc,
                          int dna_bulges,
                          int rna_bulges)
{
    int leading_ns = num_leading_ns(augment, dna_bulges);
    base_dna_match = base_dna_match.substr(leading_ns);
    base_dna_match = get_dna_match(base_dna_match, augment);
    base_rna_match = get_rna_match(base_rna_match, augment);
    indicate_mismatches_dna(base_dna_match, base_rna_match);

    return bulge_info{ .dna = base_dna_match,
                       .rna = base_rna_match,
                       .loc = orig_loc + leading_ns };
}

std::vector<bulge_pair> augment_patterns_with_bulges(
  std::vector<std::string> patterns,
  int dna_bulges,
  int rna_bulges)
{
    std::vector<bulge_pair> augmented;
    int pad_size = dna_bulges;
    for (std::string orig : patterns) {
        augmented.emplace_back(std::string(pad_size, 'N') + orig,
                               bulge_augment{ .bulge_pos = 0,
                                              .bulge_size = 0,
                                              .bulge_type = BULGE_NONE });
        for (int i : range(1, dna_bulges + 1)) {
            for (int j : range(1, orig.size())) {
                augmented.emplace_back(
                  std::string(pad_size - i, 'N') +
                    std::string(orig.begin(), orig.begin() + j) +
                    std::string(i, 'N') +
                    std::string(orig.begin() + j, orig.end()),
                  bulge_augment{
                    .bulge_pos = j, .bulge_size = i, .bulge_type = BULGE_DNA });
            }
        }
        for (int i : range(1, rna_bulges + 1)) {
            for (int j : range(i, orig.size() + 1)) {
                augmented.emplace_back(
                  std::string(pad_size + i, 'N') +
                    std::string(orig.begin(), orig.begin() + j - i) +
                    std::string(orig.begin() + j, orig.end()),
                  bulge_augment{ .bulge_pos = j - i,
                                 .bulge_size = i,
                                 .bulge_type = BULGE_RNA });
            }
        }
    }
    return augmented;
}
