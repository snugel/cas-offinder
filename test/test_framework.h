#pragma once

typedef bool (*test_fn)();
bool add_test(test_fn test, const char* func_name);
void run_all_tests(const char* name_filter = nullptr);
void notify_check_failure(const char* message);

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
    if (!(expr))                                                               \
    throw("At line: " STRINGIZE(__LINE__) " test assertion: '" #expr           \
                                          "' failed.")
#define t_check(expr)                                                          \
    if (!(expr))                                                               \
    notify_check_failure(                                                      \
      "At line: " STRINGIZE(__LINE__) " test check: '" #expr "' failed.")
