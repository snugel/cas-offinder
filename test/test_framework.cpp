#include "test_framework.h"
#include <array>
#include <exception>
#include <iostream>
#include <vector>

using namespace std;

struct TestPair
{
    const char* name;
    test_fn fn;
};
// no constructor needs to be called to use these variables--they are
// initialized statically at compile time, allowing them to be used
// any time during startup.
std::array<TestPair, MAX_NUM_TESTS> tests;
size_t cur_n_tests; // initted to 0
vector<const char*> cur_t_check_messages;
struct TestResult
{
    const char* test_name;
    bool retval;
    vector<const char*> t_check_messages;
    std::string exception_msg;
    bool asserted;
};

void run_all_tests(const char* name_filter)
{
    vector<TestResult> test_failures;
    int num_passed = 0;
    for (size_t i = 0; i < cur_n_tests; i++) {
        TestPair test = tests[i];
        cout.clear();
        try {
            cur_t_check_messages.clear();
            bool ret_val = test.fn();
            if (ret_val && cur_t_check_messages.size() == 0) {
                num_passed++;
            } else {
                test_failures.push_back(
                  TestResult{ .test_name = test.name,
                              .retval = ret_val,
                              .t_check_messages = cur_t_check_messages,
                              .asserted = false });
            }
        } catch (std::exception& except) {
            test_failures.push_back(
              TestResult{ .test_name = test.name,
                          .retval = false,
                          .t_check_messages = cur_t_check_messages,
                          .exception_msg = string(except.what()),
                          .asserted = false });
        } catch (const char* message) {
            test_failures.push_back(
              TestResult{ .test_name = test.name,
                          .retval = false,
                          .t_check_messages = cur_t_check_messages,
                          .exception_msg = string(message),
                          .asserted = true });
        } catch (...) {
            test_failures.push_back(
              TestResult{ .test_name = test.name,
                          .retval = false,
                          .t_check_messages = cur_t_check_messages,
                          .exception_msg = string("unknown exception"),
                          .asserted = false });
        }
    }
    for (TestResult& result : test_failures) {
        cerr << "test '" << result.test_name << "'\n";
        if (result.t_check_messages.size()) {
            for (const char* check : result.t_check_messages) {
                cerr << check << "\n";
            }
        }
        if (result.asserted) {
            cerr << result.exception_msg << "\n";
        } else if (result.exception_msg.size()) {
            cerr << "Exception raised: '" << result.exception_msg << "'\n";
        } else if (!result.retval) {
            cerr << "returned false\n";
        }
        cerr << "\n";
    }
    if (test_failures.size()) {
        cerr << "\n " << test_failures.size() << " TESTS FAILED \n "
             << num_passed << " TESTS PASSED" << endl;
    } else {
        cerr << "\n\nALL " << num_passed << " TESTS PASSED.\n" << endl;
    }
}

bool add_test(test_fn test, const char* func_name)
{
    tests.at(cur_n_tests) = TestPair{ .name = func_name, .fn = test };
    cur_n_tests++;
    return true;
}

void notify_check_failure(const char* message)
{
    cur_t_check_messages.push_back(message);
}

//#define TESTTEST

#ifdef TESTTEST

TEST(testtestpass)
{
    return true;
}

TEST(testtestfail)
{
    return false;
}

TEST(testtesterror)
{
    throw runtime_error("testtesterror raised error as expected");
}

static void t_check_test_helper()
{
    t_check(false && "helper check failure");
}

TEST(test_t_check)
{
    t_check(false && "first check failure");
    t_check_test_helper();
    return true;
}

TEST(test_t_assert)
{
    t_assert(false && "first check failure");
    cout << "THIS SHOULD NOT PRINT\n";
    return true;
}

#endif
