#include "bulge_logic.h"
#include "test/test_framework.h"
#include <iostream>


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