#include "parse_input.h"
#include "test/test_framework.h"
#include <cstring>

TEST(test_parse_input1)
{
    //    test/t_data
    //    NNNNNNNNNNNNNNNNNNNNNRG
    //    GGCCGACCTGTCGCTGACGCNNN 5
    //    CGCCGACCTGTCGCTGACGCNNN 5
    const char* in_fname = "../cas-offinder/test/t_data/example.in";
    auto out = read_file(in_fname);
    return !strcmp(out.genome_path, "test/t_data") && out.mismatches == 5 &&
           out.num_patterns == 2 && out.pattern_size == 23 &&
           !strcmp(out.pattern, "NNNNNNNNNNNNNNNNNNNNNRG") &&
           !strcmp(out.compares,
                   "GGCCGACCTGTCGCTGACGCNNNCGCCGACCTGTCGCTGACGCNNN") &&
           out.ids == nullptr;
}

TEST(test_parse_input2)
{
    //    test/t_data
    //    NNNNNNNNNNNNNNNNNNNNNRG
    //    GGCCGACCTGTCGCTGACGCNNN 5 IDS
    //    GGCCGACCTGTCGCTGACGCNNN 5 IDS2
    const char* in_fname = "../cas-offinder/test/t_data/example2.in";
    auto out = read_file(in_fname);
    return out.ids && !strcmp(out.ids[0], "IDS") && !strcmp(out.ids[1], "IDS2");
}
