#include "bulge_logic.h"
#include "test/test_framework.h"
#include "RangeIterator.h"
#include <iostream>

//NNACT   0       0               X
//NANCT   1       1               DNA
//NACNT   2       1               DNA
//ANNCT   1       2               DNA
//ACNNT   2       2               DNA
//NNNCT   0       1       A       RNA
//NNNAT   1       1       C       RNA
//NNNAC   2       1       T       RNA
//NNNNT   0       2       AC      RNA
//NNNNA   1       2       CT      RNA

static bool equal(bulge_info b1, bulge_info b2)
{
    std::cout << b1.dna << "\t" << b1.rna << "\t" << b1.loc << "\t\t";
    std::cout << b2.dna << "\t" << b2.rna << "\t" << b2.loc << "\n";
    return b1.loc == b2.loc && b1.dna == b2.dna && b1.rna == b2.rna;
}

constexpr int orig_off = 3;
TEST(test_get_bulge_info1)
{
    return equal(bulge_info{ .dna = "TCA", .rna = "ACT", .loc = 2 + orig_off }, get_bulge_info("AGTCA", "NNACT", bulge_augment{ .bulge_pos = 0, .bulge_size = 0, .bulge_type = BULGE_NONE }, orig_off, 2));
}
TEST(test_get_bulge_info2)
{
    return equal(bulge_info{ .dna = "GTCA", .rna = "A-CT", .loc = 1 + orig_off }, get_bulge_info("AGTCA", "NANCT", bulge_augment{ .bulge_pos = 1, .bulge_size = 1, .bulge_type = BULGE_DNA }, orig_off, 2));
}
TEST(test_get_bulge_info3)
{
    return equal(bulge_info{ .dna = "AGTCA", .rna = "AC--T", .loc = 0 + orig_off }, get_bulge_info("AGTCA", "ACNNT", bulge_augment{ .bulge_pos = 2, .bulge_size = 2, .bulge_type = BULGE_DNA }, orig_off, 2));
}
TEST(test_get_bulge_info4)
{
    return equal(bulge_info{ .dna = "-CA", .rna = "ACT", .loc = 3 + orig_off }, get_bulge_info("AGTCA", "NNNCT", bulge_augment{ .bulge_pos = 0, .bulge_size = 1, .bulge_type = BULGE_RNA, .removed_part="A" }, orig_off, 2));
}
TEST(test_get_bulge_info5)
{
    return equal(bulge_info{ .dna = "C-A", .rna = "ACT", .loc = 3 + orig_off }, get_bulge_info("AGTCA", "NNNAT", bulge_augment{ .bulge_pos = 1, .bulge_size = 1, .bulge_type = BULGE_RNA, .removed_part="C" }, orig_off, 2));
}
TEST(test_get_bulge_info6)
{
    return equal(bulge_info{ .dna = "--A", .rna = "ACT", .loc = 4 + orig_off }, get_bulge_info("AGTCA", "NNNNT", bulge_augment{ .bulge_pos = 0, .bulge_size = 2, .bulge_type = BULGE_RNA, .removed_part="AC" }, orig_off, 2));
}
TEST(test_get_bulge_info7)
{
    return equal(bulge_info{ .dna = "A--", .rna = "ACT", .loc = 4 + orig_off }, get_bulge_info("AGTCA", "NNNNA", bulge_augment{ .bulge_pos = 1, .bulge_size = 2, .bulge_type = BULGE_RNA, .removed_part="CT" }, orig_off, 2));
}


TEST(test_augment_bulge_logic)
{
    std::string pattern = "ACT";
    std::vector<bulge_pair> augmented = augment_patterns_with_bulges({ pattern }, 2, 2);
    std::vector<bulge_pair> expected = {
        { "NNACT", bulge_augment{ .bulge_pos = 0, .bulge_size = 0, .bulge_type = BULGE_NONE } }, { "NANCT", bulge_augment{ .bulge_pos = 1, .bulge_size = 1, .bulge_type = BULGE_DNA } },
        { "NACNT", bulge_augment{ .bulge_pos = 2, .bulge_size = 1, .bulge_type = BULGE_DNA } },  { "ANNCT", bulge_augment{ .bulge_pos = 1, .bulge_size = 2, .bulge_type = BULGE_DNA } },
        { "ACNNT", bulge_augment{ .bulge_pos = 2, .bulge_size = 2, .bulge_type = BULGE_DNA } },  { "NNNCT", bulge_augment{ .bulge_pos = 0, .bulge_size = 1, .bulge_type = BULGE_RNA, .removed_part="A" } },
        { "NNNAT", bulge_augment{ .bulge_pos = 1, .bulge_size = 1, .bulge_type = BULGE_RNA, .removed_part="C" } },  { "NNNAC", bulge_augment{ .bulge_pos = 2, .bulge_size = 1, .bulge_type = BULGE_RNA, .removed_part="T" } },
        { "NNNNT", bulge_augment{ .bulge_pos = 0, .bulge_size = 2, .bulge_type = BULGE_RNA, .removed_part="AC" } },  { "NNNNA", bulge_augment{ .bulge_pos = 1, .bulge_size = 2, .bulge_type = BULGE_RNA, .removed_part="CT" } },
    };
    t_assert(augmented.size() == expected.size());
    for(size_t i : range(augmented.size())){
        {
            auto p = augmented[i];
            std::cout << p.first << "\t" << p.second.bulge_pos << "\t" << p.second.bulge_size  << "\t" << p.second.removed_part << "\t" << get_bulge_type_name(p.second.bulge_type) << "\n";
        }
//        {
//            auto p = expected[i];
//            std::cout << p.first << "\t" << p.second.bulge_pos << "\t" << p.second.bulge_size << "\t" << p.second.removed_part << "\t" << get_bulge_type_name(p.second.bulge_type) << "\n";
//        }
        bulge_pair a = augmented[i];
        bulge_pair b = expected[i];
        t_check(a.first == b.first);
        t_check(a.second.bulge_pos == b.second.bulge_pos);
        t_check(a.second.bulge_type == b.second.bulge_type);
        t_check(a.second.removed_part == b.second.removed_part);
        t_check(a.second.bulge_size == b.second.bulge_size);
    }
    return true;
}
TEST(test_reverse_compliment)
{
    return reverse_compliment("AGCTN") == "NAGCT";
}
