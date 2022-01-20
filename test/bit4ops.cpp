#include "bit4ops.h"
#include "test/test_framework.h"
#include <iostream>

TEST(test_byte_transform){
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint32_t> genomeb4 = make4bitpackedint32(genome);
    std::vector<uint64_t> actual = bit64tobit32(genomeb4);
    std::vector<uint64_t> expected = make4bitpackedint64(genome);
    for(uint64_t i : expected){
        std::cout << std::hex << i << "\t";
    }
    std::cout << "\n";
    for(uint64_t i : actual){
        std::cout << std::hex << i << "\t";
    }
    std::cout << "\n";
    return expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin());
}
TEST(test_make4bitpackedint64)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint64_t> expected = {0x4282814842841248, 0x1284128148214812, 0x1841800000000000};
    std::vector<uint64_t> actual = make4bitpackedint64(genome);
    return expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin());
}

TEST(test_make4bitpackedint32)
{
    std::string genome = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint32_t> expected = {0x42828148, 0x42841248, 0x12841281, 0x48214812, 0x18418000};
    std::vector<uint32_t> actual = make4bitpackedint32(genome);
    return expected.size() == actual.size() && std::equal(expected.begin(), expected.end(), actual.begin());
}


TEST(test_bit4tostr)
{
    //data = "ACGCGTAGACGATCAGTCGATCGTAGCTAGTCTGATG";
    std::vector<uint32_t> data = {0x42828148, 0x42841248, 0x12841281, 0x48214812, 0x18418000};
    size_t start = 5;
    size_t end = 30;
    std::string expected = "TAGACGATCAGTCGATCGTAGCTAG";
    std::string actual = bit4tostr(data, start, end);
    std::cout << expected << "\n" << actual << "\n";
    return expected == actual;
}
