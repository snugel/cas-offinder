#include "bit4ops.h"
#include <array>
#include <stddef.h>

constexpr uint8_t T = 0x1;
constexpr uint8_t C = 0x2;
constexpr uint8_t A = 0x4;
constexpr uint8_t G = 0x8;
static std::array<uint8_t, 256> makebit4patternmap(){
    std::array<uint8_t, 256> arr;
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
static std::array<uint8_t, 256> makebit4map(){
    std::array<uint8_t, 256> arr;
    // all other values, including 'N' mapped to 0
    std::fill(arr.begin(), arr.end(), 0);
    arr['G'] = G;
    arr['C'] = C;
    arr['A'] = A;
    arr['T'] = T;
    return arr;
}
std::array<uint8_t, 256> tobit4patternmap = makebit4patternmap();
std::array<uint8_t, 256> tobit4map = makebit4map();

void str2bit4_impl(uint8_t * dest, const char * src, int64_t write_offset, uint64_t n_chrs, uint8_t * arrdata){
    dest += write_offset / 2;
    write_offset %= 2;
    if(write_offset && n_chrs > 0){
        dest[0] |= arrdata[uint8_t(src[0])] << 4; 
        dest += 1;
        src += 1;
        n_chrs -= 1;
    }
    for(size_t i = 0; i < n_chrs/2; i++){
        dest[i] = arrdata[uint8_t(src[i*2])] | arrdata[uint8_t(src[i*2+1])] << 4;
    }
    if(n_chrs%2){
        dest[n_chrs/2] |= arrdata[uint8_t(src[n_chrs-1])];
    }
}

void str2bit4_pattern(uint8_t * dest, const char * src, int64_t write_offset, uint64_t n_chrs){
    str2bit4_impl(dest, src, write_offset, n_chrs, &tobit4patternmap[0]);
}
void str2bit4(uint8_t * dest, const char * src, int64_t write_offset, uint64_t n_chrs){
    str2bit4_impl(dest, src, write_offset, n_chrs, &tobit4map[0]);
}
char bit42chr(uint8_t v){
    // currently no support for patterns, as it is not needed
    switch(v){
        case C: return 'C';
        case G: return 'G';
        case T: return 'T';
        case A: return 'A';
        default: return 'N';// to capture weird values
    }
}
void bit42str(char * dest, const uint8_t * src, int64_t read_offset, uint64_t n_chrs){
    src += read_offset / 2;
    read_offset %= 2;
    if(read_offset && n_chrs > 0){
        dest[0] = bit42chr(src[0] >> 4);
        dest += 1;
        src += 1;
        n_chrs -= 1;
    }
    for(size_t i = 0; i < n_chrs/2; i++){
        dest[i*2] = bit42chr(src[i] & 0xf);
        dest[i*2+1] = bit42chr(src[i] >> 4);
    }
    if(n_chrs % 2){
        dest[n_chrs-1] = bit42chr(src[n_chrs/2] & 0xf);
    } 
}
