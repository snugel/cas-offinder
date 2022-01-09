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
    int orig_idx;
    int bulge_pos;
    BulgeType bulge_type;
};
void get_bulge_info(bulge_info & result, std::string & genome, bulge_augment & augment, int dna_bulges, int rna_bulges);
std::vector<std::pair<std::string,bulge_augment>> augment_patterns_with_bulges(std::vector<std::string> & patterns, int dna_bulges, int rna_bulges);
