#include "bulge_logic.h"
#include "test/test_framework.h"
#include <iostream>

static bool equal(bulge_info b1, bulge_info b2)
{
    std::cout << b1.dna << "\t" << b1.rna << "\t" << b1.loc << "\t\t";
    std::cout << b2.dna << "\t" << b2.rna << "\t" << b2.loc << "\n";
    return b1.loc == b2.loc && b1.dna == b2.dna && b1.rna == b2.rna;
}

constexpr int orig_off = 3;
TEST(test_get_bulge_info1)
{
    return equal(bulge_info{ .dna = "tCa", .rna = "ACT", .loc = 2 + orig_off },
                 get_bulge_info("AGTCA",
                                "ACT",
                                bulge_augment{ .bulge_pos = 0,
                                               .bulge_size = 0,
                                               .bulge_type = BULGE_NONE },
                                orig_off,
                                2,
                                2));
}
TEST(test_get_bulge_info2)
{
    return equal(
      bulge_info{ .dna = "gTCa", .rna = "A-CT", .loc = 1 + orig_off },
      get_bulge_info("AGTCA",
                     "ACT",
                     bulge_augment{ .bulge_pos = 1,
                                    .bulge_size = 1,
                                    .bulge_type = BULGE_DNA },
                     orig_off,
                     2,
                     2));
}
TEST(test_get_bulge_info3)
{
    return equal(
      bulge_info{ .dna = "AgTCa", .rna = "AC--T", .loc = 0 + orig_off },
      get_bulge_info("AGTCA",
                     "ACT",
                     bulge_augment{ .bulge_pos = 2,
                                    .bulge_size = 2,
                                    .bulge_type = BULGE_DNA },
                     orig_off,
                     2,
                     2));
}
TEST(test_get_bulge_info4)
{
    return equal(bulge_info{ .dna = "-Ca", .rna = "ACT", .loc = 3 + orig_off },
                 get_bulge_info("AGTCA",
                                "ACT",
                                bulge_augment{ .bulge_pos = 0,
                                               .bulge_size = 1,
                                               .bulge_type = BULGE_RNA },
                                orig_off,
                                2,
                                2));
}
TEST(test_get_bulge_info5)
{
    return equal(bulge_info{ .dna = "c-a", .rna = "ACT", .loc = 3 + orig_off },
                 get_bulge_info("AGTCA",
                                "ACT",
                                bulge_augment{ .bulge_pos = 1,
                                               .bulge_size = 1,
                                               .bulge_type = BULGE_RNA },
                                orig_off,
                                2,
                                2));
}
TEST(test_get_bulge_info6)
{
    return equal(bulge_info{ .dna = "--a", .rna = "ACT", .loc = 4 + orig_off },
                 get_bulge_info("AGTCA",
                                "ACT",
                                bulge_augment{ .bulge_pos = 0,
                                               .bulge_size = 2,
                                               .bulge_type = BULGE_RNA },
                                orig_off,
                                2,
                                2));
}
TEST(test_get_bulge_info7)
{
    return equal(bulge_info{ .dna = "A--", .rna = "ACT", .loc = 4 + orig_off },
                 get_bulge_info("AGTCA",
                                "ACT",
                                bulge_augment{ .bulge_pos = 1,
                                               .bulge_size = 2,
                                               .bulge_type = BULGE_RNA },
                                orig_off,
                                2,
                                2));
}

TEST(test_augment_bulge_logic)
{
    std::string pattern = "ACT";
    std::vector<bulge_pair> augmented =
      augment_patterns_with_bulges({ pattern }, 2, 2);
    std::vector<bulge_pair> expected = {
        { "NNACT",
          bulge_augment{
            .bulge_pos = 0, .bulge_size = 0, .bulge_type = BULGE_NONE } },
        { "NANCT",
          bulge_augment{
            .bulge_pos = 1, .bulge_size = 1, .bulge_type = BULGE_DNA } },
        { "NACNT",
          bulge_augment{
            .bulge_pos = 2, .bulge_size = 1, .bulge_type = BULGE_DNA } },
        { "ANNCT",
          bulge_augment{
            .bulge_pos = 1, .bulge_size = 2, .bulge_type = BULGE_DNA } },
        { "ACNNT",
          bulge_augment{
            .bulge_pos = 2, .bulge_size = 2, .bulge_type = BULGE_DNA } },
        { "NNNCT",
          bulge_augment{
            .bulge_pos = 0, .bulge_size = 1, .bulge_type = BULGE_RNA } },
        { "NNNAT",
          bulge_augment{
            .bulge_pos = 1, .bulge_size = 1, .bulge_type = BULGE_RNA } },
        { "NNNAC",
          bulge_augment{
            .bulge_pos = 2, .bulge_size = 1, .bulge_type = BULGE_RNA } },
        { "NNNNT",
          bulge_augment{
            .bulge_pos = 0, .bulge_size = 2, .bulge_type = BULGE_RNA } },
        { "NNNNA",
          bulge_augment{
            .bulge_pos = 1, .bulge_size = 2, .bulge_type = BULGE_RNA } },
    };
    for (auto p : augmented) {
        std::cout << p.first << "\t" << p.second.bulge_pos << "\t"
                  << p.second.bulge_size << "\t"
                  << get_bulge_type_name(p.second.bulge_type) << "\n";
    }
    return augmented.size() == expected.size() &&
           std::equal(augmented.begin(),
                      augmented.end(),
                      expected.begin(),
                      [](bulge_pair a, bulge_pair b) {
                          return a.first == b.first &&
                                 a.second.bulge_pos == b.second.bulge_pos &&
                                 a.second.bulge_type == b.second.bulge_type &&
                                 a.second.bulge_size == b.second.bulge_size;
                      });
}
TEST(test_reverse_compliment)
{
    return reverse_compliment("AGCTN") == "NAGCT";
}
