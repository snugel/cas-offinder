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

static int num_leading_ns(bulge_augment& augment, int dna_bulges)
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
                   base_rna_match.substr(augment.bulge_pos +
                                         augment.bulge_size);
        case BULGE_RNA:
            return base_rna_match.substr(0, augment.bulge_pos) +
                   augment.removed_part +
                   base_rna_match.substr(augment.bulge_pos);
        case BULGE_NONE: return base_rna_match;
    }
    throw std::runtime_error("invalid bulge type");
}

static std::string get_dna_match(std::string base_dna_match,
                                 bulge_augment augment,
                                 int dna_bulges)
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
        char c = to_upper(seq[i]);
        if (c == 'A')
            seq[i] = 'T';
        else if (c == 'T')
            seq[i] = 'A';
        else if (c == 'G')
            seq[i] = 'C';
        else if (c == 'C')
            seq[i] = 'G';
        else if (c == 'R')
            seq[i] = 'Y';
        else if (c == 'Y')
            seq[i] = 'R';
        else if (c == 'M')
            seq[i] = 'K';
        else if (c == 'K')
            seq[i] = 'M';
        else if (c == 'H')
            seq[i] = 'D';
        else if (c == 'D')
            seq[i] = 'H';
        else if (c == 'B')
            seq[i] = 'V';
        else if (c == 'V')
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
    base_rna_match = base_rna_match.substr(leading_ns);
    base_dna_match = get_dna_match(base_dna_match, augment, dna_bulges);
    base_rna_match = get_rna_match(base_rna_match, augment);

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
                                              .bulge_type = BULGE_NONE,
                                              .removed_part = "" });
        for (int i : range(1, dna_bulges + 1)) {
            for (int j : range(1, orig.size())) {
                augmented.emplace_back(
                  std::string(pad_size - i, 'N') +
                    std::string(orig.begin(), orig.begin() + j) +
                    std::string(i, 'N') +
                    std::string(orig.begin() + j, orig.end()),
                  bulge_augment{ .bulge_pos = j,
                                 .bulge_size = i,
                                 .bulge_type = BULGE_DNA,
                                 .removed_part = "" });
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
                                 .bulge_type = BULGE_RNA,
                                 .removed_part = orig.substr(j - i, i) });
            }
        }
    }
    return augmented;
}

size_t get_num_augments(int pattern_size, int dna_bulges, int rna_bulges)
{
    std::vector<std::string> pats{ std::string(pattern_size, 'N') };
    assert(pats.size() == 1);
    std::vector<bulge_pair> augs =
      augment_patterns_with_bulges(pats, dna_bulges, rna_bulges);
    return augs.size();
}
bool is_reversed_pam(std::string pattern)
{
    // As of today (2021.07.29), there is no reversed PAM starts with 'N'
    return pattern.at(0) != 'N';
}

std::string space_out_pattern(std::string pattern, int dna_bulge_size)
{
    if (is_reversed_pam(pattern)) {
        return pattern + std::string(dna_bulge_size, 'N');
    } else {
        return std::string(dna_bulge_size, 'N') + pattern;
    }
}
