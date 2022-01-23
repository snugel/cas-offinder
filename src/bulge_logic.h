#pragma once
#include<vector>
#include<string>
#include"find_mismatches.h"

enum BulgeType{BULGE_DNA, BULGE_RNA, BULGE_NONE};
struct bulge_info{
    std::string dna;
    std::string rna;
    int loc;
};
struct bulge_augment{
    int bulge_pos;
    int bulge_size;
    BulgeType bulge_type;
};
using bulge_pair = std::pair<std::string,bulge_augment>;
std::string get_bulge_type_name(BulgeType type);
bulge_info get_bulge_info(std::string base_dna_match, std::string base_rna_match, bulge_augment augment, int orig_loc, int dna_bulges, int rna_bulges);
std::vector<bulge_pair> augment_patterns_with_bulges(std::vector<std::string> patterns, int dna_bulges, int rna_bulges);
std::string reverse_compliment(std::string s);
