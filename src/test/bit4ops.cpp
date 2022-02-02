#include "bit4ops.h"
#include "test/test_framework.h"
#include <string.h>

TEST(test_str2bit4){
    char inpt_arr[] = "ACTGC";
    uint8_t expected_out[] = {0x24,0x81,0x02};
    uint8_t actual_out[sizeof(expected_out)] = {0};
    str2bit4(actual_out, inpt_arr, 0, sizeof(inpt_arr));
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_str2bit4_offset_1){
    int offset = 1;
    char inpt_arr[] = "ACTGC";
    uint8_t expected_out[] = {0x40,0x12,0x28};
    uint8_t actual_out[sizeof(expected_out)] = {0};
    str2bit4(actual_out, inpt_arr, offset, sizeof(inpt_arr)-offset);
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_str2bit4_offset_3){
    int offset = 3;
    char inpt_arr[] = "ACTGC";
    uint8_t expected_out[] = {0x00,0x40,0x12, 0x28,0x00};
    uint8_t actual_out[sizeof(expected_out)] = {0};
    str2bit4(actual_out, inpt_arr, offset, sizeof(inpt_arr));
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_bit42str){
    uint8_t input[] = {0x24,0x81,0x02};
    char expected_out[] = "ACTGC";
    constexpr int out_size = sizeof(expected_out)-1;
    char actual_out[out_size+1] = {0};
    bit42str(actual_out, input, 0, out_size);

    return !memcmp(expected_out, actual_out, out_size);
}

TEST(test_bit42str_offset1){
    uint8_t input[] = {0x40,0x12,0x28};
    char expected_out[] = "ACTGC";
    constexpr int out_size = sizeof(expected_out)-1;
    char actual_out[out_size] = {0};
    bit42str(actual_out, input, 1, out_size);
    return !memcmp(expected_out, actual_out, out_size);
}

TEST(test_bit42str_offset3){
    uint8_t input[] = {0x00, 0x40,0x12,0x28};
    char expected_out[] = "ACTGC";
    constexpr int out_size = sizeof(expected_out)-1;
    char actual_out[out_size] = {0};
    bit42str(actual_out, input, 3, out_size);
    return !memcmp(expected_out, actual_out, out_size);
}

