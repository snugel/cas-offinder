#include "bit4ops.h"
#include "test/test_framework.h"

inline std::array<char, 256> make4bitmap(){
    constexpr char T = 0x1;
    constexpr char C = 0x2;
    constexpr char A = 0x4;
    constexpr char G = 0x8;
    std::array<char, 256> arr;
    arr['G'] = G;
    arr['C'] = C;
    arr['A'] = A;
    arr['T'] = T;
    arr['R'] = A | G;
    arr['Y'] = C | T;
    arr['S'] = G | C;
    arr['W'] = A | T;
    arr['K'] = G | T;
    arr['M'] = A | C;
    arr['B'] = C | G | T;
    arr['D'] = A | G | T;
    arr['H'] = A | C | T;
    arr['V'] = A | C | G;
    arr['N'] = A | C | G | T;
    return arr;
}
std::array<char, 256> to4bitmap = make4bitmap();


template<typename int_ty>
std::vector<int_ty> make4bitpackedint_generic(const std::string & genome)
{
    constexpr size_t bit4_c = sizeof(int_ty) * 8 / 4;
    std::vector<int_ty> result((genome.size() + bit4_c - 1) / bit4_c);
    
    for (size_t i = 0; i < genome.size() / bit4_c; i++)
    {
        for(size_t j = 0; j < bit4_c; j++){
            int shift = bit4_c - 1 - j;
            result[i] |= int_ty(to4bit(genome[i*bit4_c+j])) << (shift * 4);
        }
    }
    size_t i = genome.size() / bit4_c;
    for(size_t j = 0; j < genome.size() - i*bit4_c; j++){
        int shift = bit4_c - 1 - j;
        result[i] |= int_ty(to4bit(genome[i*bit4_c+j])) << (shift * 4);
    }
    return result;
}
std::vector<uint64_t> make4bitpackedint64(const std::string & genome){
    return make4bitpackedint_generic<uint64_t>(genome);
}
std::vector<uint32_t> make4bitpackedint32(const std::string & genome){
    return make4bitpackedint_generic<uint32_t>(genome);
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
