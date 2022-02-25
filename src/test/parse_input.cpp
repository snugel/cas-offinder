#include "parse_input.h"
#include "test/test_framework.h"
#include <cstring>

TEST(test_parse_input1)
{
    //    test/t_data
    //    NNNNNNNNNNNNNNNNNNNNNRG
    //    GGCCGACCTGTCGCTGACGCNNN 5
    //    CGCCGACCTGTCGCTGACGCNNN 5
    const char* in_fname = "../cas-offinder/src/test/t_data/example.in";
    auto out = read_file(in_fname);
    t_check(!strcmp(out.genome_path, "../cas-offinder/src/test/t_data"));
    t_check(out.mismatches == 7);
    t_check(out.num_patterns == 2 && out.pattern_size == 23);
    t_check(!strcmp(out.pattern, "NNNNNNNNNNNNNNNNNNNNNRG"));
    t_check(
      !strcmp(out.compares, "GGCCGACCTGTCGCTGACGCNNNCGCCGACCTGTCGCTGACGCNNN"));
    t_check(out.ids == nullptr);
    t_check(out.dna_bulges == 0);
    t_check(out.rna_bulges == 0);
    return true;
}

TEST(test_parse_input2)
{
    //    test/t_data
    //    NNNNNNNNNNNNNNNNNNNNNRG
    //    GGCCGACCTGTCGCTGACGCNNN 5 IDS
    //    GGCCGACCTGTCGCTGACGCNNN 5 IDS2
    const char* in_fname = "../cas-offinder/src/test/t_data/example2.in";
    auto out = read_file(in_fname);
    return out.ids && !strcmp(out.ids[0], "IDS") && !strcmp(out.ids[1], "IDS2");
}

TEST(test_parse_inputv3)
{
    //    test/t_data
    //    NNNNNNNNNNNNNNNNNNNNNRG 1 2
    //    GGCCGACCTGTCGCTGACGCNNN 5
    //    CGCCGACCTGTCGCTGACGCNNN 5
    const char* in_fname = "../cas-offinder/src/test/t_data/examplev3.in";
    auto out = read_file(in_fname);
    t_check(!strcmp(out.genome_path, "../cas-offinder/src/test/t_data"));
    t_check(!strcmp(out.pattern, "NNNNNNNNNNNNNNNNNNNNNRG"));
    t_check(out.ids == nullptr);
    t_check(out.dna_bulges == 1);
    t_check(out.rna_bulges == 2);
    return true;
}
