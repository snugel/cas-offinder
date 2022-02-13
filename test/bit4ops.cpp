#include "bit4ops.h"
#include "test_framework.h"
#include <string.h>

TEST(test_str2bit4)
{
    char inpt_arr[] = "ACTGC";
    uint8_t expected_out[] = { 0x24, 0x81, 0x02 };
    uint8_t actual_out[sizeof(expected_out)] = { 0 };
    str2bit4(actual_out, inpt_arr, 0, sizeof(inpt_arr));
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_str2bit4_offset_1)
{
    int offset = 1;
    char inpt_arr[] = "ACTGC";
    uint8_t expected_out[] = { 0x40, 0x12, 0x28 };
    uint8_t actual_out[sizeof(expected_out)] = { 0 };
    str2bit4(actual_out, inpt_arr, offset, sizeof(inpt_arr) - offset);
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_str2bit4_offset_3)
{
    int offset = 3;
    char inpt_arr[] = "ACTGC";
    uint8_t expected_out[] = { 0x00, 0x40, 0x12, 0x28, 0x00 };
    uint8_t actual_out[sizeof(expected_out)] = { 0 };
    str2bit4(actual_out, inpt_arr, offset, sizeof(inpt_arr));
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_str2bit4_large)
{
    char inpt_arr[] = "ACTGCAACTGCA";
    uint8_t expected_out[] = { 0x24, 0x81, 0x42, 0x24, 0x81, 0x42 };
    uint8_t actual_out[sizeof(expected_out)] = { 0 };
    str2bit4(actual_out, inpt_arr, 0, sizeof(inpt_arr));
    return !memcmp(expected_out, actual_out, sizeof(expected_out));
}

TEST(test_bit42str)
{
    uint8_t input[] = { 0x24, 0x81, 0x02 };
    char expected_out[] = "ACTGC";
    constexpr int out_size = sizeof(expected_out) - 1;
    char actual_out[out_size + 1] = { 0 };
    bit42str(actual_out, input, 0, out_size);

    return !memcmp(expected_out, actual_out, out_size);
}

TEST(test_bit42str_offset1)
{
    uint8_t input[] = { 0x40, 0x12, 0x28 };
    char expected_out[] = "ACTGC";
    constexpr int out_size = sizeof(expected_out) - 1;
    char actual_out[out_size] = { 0 };
    bit42str(actual_out, input, 1, out_size);
    return !memcmp(expected_out, actual_out, out_size);
}

TEST(test_bit42str_offset3)
{
    uint8_t input[] = { 0x00, 0x40, 0x12, 0x28 };
    char expected_out[] = "ACTGC";
    constexpr int out_size = sizeof(expected_out) - 1;
    char actual_out[out_size] = { 0 };
    bit42str(actual_out, input, 3, out_size);
    return !memcmp(expected_out, actual_out, out_size);
}

TEST(test_memsetbit4_middle)
{
    uint8_t input[] = { 0x24, 0x81, 0x42, 0x02 };
    uint8_t expected[] = { 0x04, 0x00, 0x40, 0x02 };
    memsetbit4(input, 0, 1, 5);
    return !memcmp(input, expected, sizeof(expected));
}

TEST(test_memsetbit4_edge)
{
    uint8_t input[] = { 0x24, 0x81, 0x42, 0x02 };
    uint8_t expected[] = { 0x04, 0x80, 0x42, 0x02 };
    memsetbit4(input, 0, 1, 3);
    return !memcmp(input, expected, sizeof(expected));
}

TEST(test_memsetbit4_chnk)
{
    uint8_t input[] = { 0x24, 0x81, 0x42, 0x02 };
    uint8_t expected[] = { 0x24, 0xff, 0x42, 0x02 };
    memsetbit4(input, 0xf, 2, 4);
    return !memcmp(input, expected, sizeof(expected));
}

TEST(test_twobit2bit4)
{
    // ACTGTGAC
    uint8_t expected[] = { 0x24, 0x81, 0x81, 0x24 };
    uint8_t input[] = { (2 << 6) | (1 << 4) | (0 << 2) | (3 << 0),
                        (0 << 6) | (3 << 4) | (2 << 2) | (1 << 0) };
    uint8_t result[sizeof(expected)] = { 0 };
    twobit2bit4(result, input, sizeof(input) * 4);
    return !memcmp(result, expected, sizeof(expected));
}

TEST(test_is_mixedbase_true)
{
    char data[] = "ACTNRVG";
    return is_mixedbase(data, sizeof(data) - 1);
}

TEST(test_is_mixedbase_false)
{
    char data[] = "ACTNRVQG";
    return !is_mixedbase(data, sizeof(data) - 1);
}

TEST(test_to_upper)
{
    char data[] = "acao\0TRd";
    char expected[] = "ACAO\0TRD";
    for (size_t i = 0; i < sizeof(data); i++) {
        data[i] = to_upper(data[i]);
    }
    return !memcmp(data, expected, sizeof(data));
}
