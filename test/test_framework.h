#pragma once

typedef bool (*test_fn)();
bool add_test(test_fn test, const char* func_name);
void run_all_tests();
bool notify_check_failure(const char* message);
bool _throw_assert(const char* message);

// make sure asserts are kept in code
#undef NDEBUG

#define TOKENPASTE(x, y) x##y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)
#define STRINGIFY(x) #x

#define STRINGIZE_DETAIL(x) #x
#define STRINGIZE(x) STRINGIZE_DETAIL(x)

#define TEST_HELP(test_name)                                                   \
    static bool test_name();                                                   \
    bool TOKENPASTE2(test_name, __LINE__) =                                    \
      add_test(test_name, STRINGIFY(test_name));                               \
    static bool test_name()
#define TEST(test_name) TEST_HELP(test_name)
// #define TESTTEST
//#ifdef RUN_TESTS
#define MAX_NUM_TESTS 10000
#define t_assert(expr)                                                         \
    !!(expr) || _throw_assert("At line: " STRINGIZE(__LINE__) " test assertion: '" #expr "' failed.")
#define t_check(expr)                                                          \
    !!(expr) || notify_check_failure("At line: " STRINGIZE(__LINE__) " test check: '" #expr "' failed.")
