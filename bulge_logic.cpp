#include "bulge_logic.h"
#include "RangeIterator.h"
#include "test/test_framework.h"
#include <iostream>


std::string get_bulge_type_name(BulgeType type){
    switch(type){
        case BULGE_DNA: return "DNA";
        case BULGE_RNA: return "RNA";
        case BULGE_NONE: return "X";
        default: throw std::runtime_error("bad bulge type");
    }
}

bulge_info get_bulge_info(std::string & genome, bulge_augment & augment, int orig_loc, int dna_bulges, int rna_bulges){
    return bulge_info{.dna="placeholder",.rna="placebolder",.loc=-1};
}

std::vector<bulge_pair> augment_patterns_with_bulges(std::vector<std::string> patterns, int dna_bulges, int rna_bulges){
    std::vector<bulge_pair> augmented;
    int pad_size = dna_bulges;
    for(std::string orig : patterns){
        augmented.emplace_back(
            std::string(pad_size, 'N') + orig,
            bulge_augment{.bulge_pos=0,.bulge_size=0,.bulge_type=BULGE_NONE}
        );
        for(int i : range(1, dna_bulges + 1)){
            for(int j : range(1, orig.size())){
                augmented.emplace_back(
                    std::string(pad_size-i, 'N') + std::string(orig.begin(), orig.begin() + j) + std::string(i, 'N') + std::string(orig.begin() + j, orig.end()),
                    bulge_augment{.bulge_pos=j,.bulge_size=i,.bulge_type=BULGE_DNA}
                );
            }
        }
        for(int i : range(1, rna_bulges + 1)){
            for(int j : range(i, orig.size()+1)){
                augmented.emplace_back(
                    std::string(pad_size+i, 'N') + std::string(orig.begin(), orig.begin() + j - i) + std::string(orig.begin() + j, orig.end()),
                    bulge_augment{.bulge_pos=j,.bulge_size=i,.bulge_type=BULGE_RNA}
                );
            }
        }
    }
    return augmented;
}

TEST(test_bulge_logic){
    std::string pattern = "ACT";
    std::vector<bulge_pair> augmented = augment_patterns_with_bulges({pattern}, 2, 2);
    std::vector<bulge_pair> expected = {
        {"NNACT", bulge_augment{.bulge_pos=0, .bulge_size=0, .bulge_type=BULGE_NONE}},
        {"NANCT", bulge_augment{.bulge_pos=1, .bulge_size=1, .bulge_type=BULGE_DNA}},
        {"NACNT", bulge_augment{.bulge_pos=2, .bulge_size=1, .bulge_type=BULGE_DNA}},
        {"ANNCT", bulge_augment{.bulge_pos=1, .bulge_size=2, .bulge_type=BULGE_DNA}},
        {"ACNNT", bulge_augment{.bulge_pos=2, .bulge_size=2, .bulge_type=BULGE_DNA}},
        {"NNNCT", bulge_augment{.bulge_pos=1, .bulge_size=1, .bulge_type=BULGE_RNA}},
        {"NNNAT", bulge_augment{.bulge_pos=2, .bulge_size=1, .bulge_type=BULGE_RNA}},
        {"NNNAC", bulge_augment{.bulge_pos=3, .bulge_size=1, .bulge_type=BULGE_RNA}},
        {"NNNNT", bulge_augment{.bulge_pos=2, .bulge_size=2, .bulge_type=BULGE_RNA}},
        {"NNNNA", bulge_augment{.bulge_pos=3, .bulge_size=2, .bulge_type=BULGE_RNA}},
    };
    for(auto p : augmented){
        std::cout << p.first << "\t" << p.second.bulge_pos << "\t" << p.second.bulge_size << "\t" << get_bulge_type_name(p.second.bulge_type) << "\n";
    }
    return augmented.size() == expected.size() &&
           std::equal(augmented.begin(), augmented.end(), expected.begin(), [](bulge_pair a, bulge_pair b){
               return 
                a.first == b.first &&
                a.second.bulge_pos == b.second.bulge_pos &&
                a.second.bulge_type == b.second.bulge_type &&
                a.second.bulge_size == b.second.bulge_size;
           });
}