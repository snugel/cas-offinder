#pragma once
#include <vector>
#include <string>
#include <utility>
#include <chrono>

typedef bool (*test_fn)();

//make sure asserts are kept in code
#undef NDEBUG

#define TOKENPASTE(x, y) x ## y
#define TOKENPASTE2(x, y) TOKENPASTE(x, y)
#define STRINGIFY(x) #x

#define TEST_HELP(test_name) bool test_name(); bool TOKENPASTE2(test_name,__LINE__) = all_tests.add_test(test_name, STRINGIFY(test_name)); bool test_name()
#define TEST(test_name) TEST_HELP(test_name)
// #define TESTTEST
//#ifdef RUN_TESTS


class TestObj{
public:
    using test_ty = std::pair<test_fn,std::string>;
    std::vector<test_ty> tests;
    bool add_test(test_fn test, std::string func_name);
    void run_all();
};

extern TestObj all_tests;
