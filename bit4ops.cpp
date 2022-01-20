#include "bit4ops.h"
#include <algorithm>

inline std::array<char, 256> make4bitmap(){
    constexpr char T = 0x1;
    constexpr char C = 0x2;
    constexpr char A = 0x4;
    constexpr char G = 0x8;
    std::array<char, 256> arr;
    std::fill(arr.begin(), arr.end(), 0);
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
void clean_bogus(std::string & genome){
    for(char & c : genome){
        if(c == 'N' || c == ';'){
            c = 0;
        }
    }
}

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
std::vector<uint64_t> bit64tobit32(const std::vector<uint32_t> & bit32){   
    std::vector<uint64_t> bit64((bit32.size()+1)/2);
    uint32_t * newb64 = reinterpret_cast<uint32_t*>(&bit64[0]);
    for(size_t i = 0; i < bit64.size() - bit32.size()%2; i++){
        newb64[i*2] = bit32[i*2+1];
    
        newb64[i*2+1] = bit32[i*2];

    } 
    if(bit32.size()%2){
        bit64.back() = uint64_t(bit32.back()) << 32;
    }
    return bit64; 
}
