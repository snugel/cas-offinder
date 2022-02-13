#include "test_framework.h"
#include <array>
#include <cstdio>
#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

#ifdef WIN32
const char* STDOUT_REFRESH = "CON";
#else
const char* STDOUT_REFRESH = "/dev/tty";
#endif

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
    string stdout_s;
    string stderr_s;
};
string read_file(const char* fname)
{
    ifstream file(fname);
    stringstream ss;
    ss << file.rdbuf();
    return ss.str();
}

void run_all_tests()
{
    int num_passed = 0;
    for (size_t i = 0; i < cur_n_tests; i++) {
        TestPair test = tests[i];
        cout.clear();
        bool ret_val;
        cur_t_check_messages.clear();
        const char* stdout_red = "_tmp_redirect_stdout.txt";
        const char* stderr_red = "_tmp_redirect_stderr.txt";
        freopen(stdout_red, "w+", stdout);
        freopen(stderr_red, "w+", stderr);
        string exception_msg = "";
        bool asserted = false;
        bool passed = false;
        try {
            ret_val = test.fn();
            if (ret_val && cur_t_check_messages.size() == 0) {
                num_passed++;
                passed = true;
            }
        } catch (std::exception& except) {
            exception_msg = except.what();
        } catch (const char* message) {
            exception_msg = string(message);
            asserted = true;
        } catch (...) {
            exception_msg = "unknown exception";
        }
        fflush(stdout);
        fflush(stderr);
        string stdout_str = read_file(stdout_red);
        string stderr_str = read_file(stderr_red);
        fclose(stdout);
        fclose(stderr);
        freopen(STDOUT_REFRESH, "w", stderr);
        freopen(STDOUT_REFRESH, "w", stdout);
        if (!passed) {
            cerr << "========test '" << test.name << "'========\n";
            if (cur_t_check_messages.size()) {
                for (const char* check : cur_t_check_messages) {
                    cerr << check << "\n";
                }
            }
            if (asserted) {
                cerr << exception_msg << "\n";
            } else if (exception_msg.size()) {
                cerr << "Exception raised: '" << exception_msg << "'\n";
            } else if (!ret_val) {
                cerr << "returned false\n";
            }
            if (stderr_str.size()) {
                cerr << "\nCAPTURED STDERR\n";
                cerr << "---------------\n";
                cerr << stderr_str;
            }
            if (stdout_str.size()) {
                cerr << "\nCAPTURED STDOUT\n";
                cerr << "---------------\n";
                cerr << stdout_str;
            }
            cerr << "\n";
        }
    }
    cerr << "========Summary========\n";
    int n_failures = cur_n_tests - num_passed;
    if (n_failures) {
        cerr << "\n " << n_failures << " TESTS FAILED \n " << num_passed
             << " TESTS PASSED" << endl;
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

bool notify_check_failure(const char* message)
{
    cur_t_check_messages.push_back(message);
    // return value shouldn't really matter...
    return true;
}
bool _throw_assert(const char* message)
{
    throw message;
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

TEST(test_t_stderr)
{
    cerr << "printing to stderr\n";
    return false;
}

TEST(test_t_stderr_noprint)
{
    cerr << "THIS SHOULD NOT PRINT\n";
    return true;
}

TEST(test_t_stdout)
{
    cout << "printing to stdout\n";
    return false;
}

TEST(test_t_stdout_noprint)
{
    cout << "THIS SHOULD NOT PRINT\n";
    return true;
}
#endif
