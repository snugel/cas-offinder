#include "bit4ops.h"
#include <array>
#include <cassert>
#include <cstring>
#include <stddef.h>

constexpr uint8_t T = 0x1;
constexpr uint8_t C = 0x2;
constexpr uint8_t A = 0x4;
constexpr uint8_t G = 0x8;
static std::array<uint8_t, 256> makebit4patternmap()
{
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
static std::array<uint8_t, 256> makebit4map()
{
    std::array<uint8_t, 256> arr;
    // all other values, including 'N' mapped to 0
    std::fill(arr.begin(), arr.end(), 0);
    arr['G'] = G;
    arr['C'] = C;
    arr['A'] = A;
    arr['T'] = T;
    return arr;
}
static std::array<uint8_t, 256> tobit4patternmap = makebit4patternmap();
static std::array<uint8_t, 256> tobit4map = makebit4map();

void str2bit4_impl(uint8_t* dest,
                   const char* src,
                   int64_t write_offset,
                   uint64_t n_chrs,
                   uint8_t* arrdata)
{
    dest += write_offset / 2;
    write_offset %= 2;
    if (write_offset && n_chrs > 0) {
        dest[0] |= arrdata[uint8_t(src[0])] << 4;
        dest += 1;
        src += 1;
        n_chrs -= 1;
    }
    for (size_t i = 0; i < n_chrs / 2; i++) {
        dest[i] =
          arrdata[uint8_t(src[i * 2])] | arrdata[uint8_t(src[i * 2 + 1])] << 4;
    }
    if (n_chrs % 2) {
        dest[n_chrs / 2] |= arrdata[uint8_t(src[n_chrs - 1])];
    }
}

void str2bit4pattern(uint8_t* dest,
                     const char* src,
                     int64_t write_offset,
                     uint64_t n_chrs)
{
    str2bit4_impl(dest, src, write_offset, n_chrs, &tobit4patternmap[0]);
}
void str2bit4(uint8_t* dest,
              const char* src,
              int64_t write_offset,
              uint64_t n_chrs)
{
    str2bit4_impl(dest, src, write_offset, n_chrs, &tobit4map[0]);
}
char bit42chr(uint8_t v)
{
    // currently no support for patterns, as it is not needed
    switch (v) {
        case C:
            return 'C';
        case G:
            return 'G';
        case T:
            return 'T';
        case A:
            return 'A';
        default:
            return 'N'; // to capture weird values
    }
}
void bit42str(char* dest,
              const uint8_t* src,
              int64_t read_offset,
              uint64_t n_chrs)
{
    src += read_offset / 2;
    read_offset %= 2;
    if (read_offset && n_chrs > 0) {
        dest[0] = bit42chr(src[0] >> 4);
        dest += 1;
        src += 1;
        n_chrs -= 1;
    }
    for (size_t i = 0; i < n_chrs / 2; i++) {
        dest[i * 2] = bit42chr(src[i] & 0xf);
        dest[i * 2 + 1] = bit42chr(src[i] >> 4);
    }
    if (n_chrs % 2) {
        dest[n_chrs - 1] = bit42chr(src[n_chrs / 2] & 0xf);
    }
}

void memsetbit4(uint8_t* dest, uint8_t bit4val, uint64_t start, uint64_t end)
{
    assert(bit4val <= 0xf);
    if (start < end && start % 2) {
        dest[start / 2] &= 0xf;
        dest[start / 2] |= bit4val << 4;
        start += 1;
    }
    if (start < end && end % 2) {
        dest[end / 2] &= 0xf0;
        dest[end / 2] |= bit4val;
        end -= 1;
    }
    if (start < end) {
        size_t bstart = start / 2;
        size_t bend = end / 2;
        uint8_t bval = bit4val | (bit4val << 4);
        memset(dest + bstart, bval, bend - bstart);
    }
}
char base_mapping[] = { 'T', 'C', 'A', 'G' };
static std::array<uint16_t, 256> create_block_map()
{
    std::array<uint16_t, 256> arr;
    for (size_t i = 0; i < 256; i++) {
        uint16_t val = 0;
        for (size_t j = 0; j < 4; j++) {
            val |= tobit4map[base_mapping[((i >> ((4 - j - 1) * 2)) & 0x3)]]
                   << (j * 4);
        }
        arr[i] = val;
    }
    return arr;
}
static std::array<uint16_t, 256> adv_mapping = create_block_map();

void twobit2bit4(uint8_t* dest, const uint8_t* src, uint64_t n_chrs)
{
    assert(n_chrs % 4 == 0 && "twobit2bit4 only support block of 4 writes");
    assert(reinterpret_cast<size_t>(dest) % 2 == 0 &&
           "dest must be alligned to 2 byte boundaries");
    uint16_t* blkdest = reinterpret_cast<uint16_t*>(dest);
    uint64_t n_blks = n_chrs / 4;
    for (size_t i = 0; i < n_blks; i++) {
        blkdest[i] = adv_mapping[src[i]];
    }
}
bool is_mixedbase(const char* src, uint64_t n_chrs)
{
    for (size_t i = 0; i < n_chrs; i++) {
        if (!tobit4patternmap[src[i]]) {
            return false;
        }
    }
    return true;
}
bool is_match(char nucl, char pattern){
    return 0 != (tobit4map[uint8_t(nucl)] & tobit4patternmap[uint8_t(pattern)]);
}
