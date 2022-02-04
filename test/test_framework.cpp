#include "test_framework.h"
#include <exception>
#include <iostream>

using namespace std;

void TestObj::run_all()
{
    int num_failed = 0;
    int num_passed = 0;

    for (test_ty& test : tests) {
        cout.clear();
        try {
            cout << test.second << ": " << endl;
            if (test.first()) {
                cout << "PASSED" << endl;
                num_passed++;
            } else {
                cout << "FAILED" << endl;
                num_failed++;
            }
        } catch (std::exception& except) {
            cout << "EXCEPTION RAISED\n" << except.what() << endl;
            num_failed++;
        } catch (...) {
            cout << "threw some weird exception" << endl;
            num_failed++;
        }
    }
    if (num_failed) {
        cout << "\n " << num_failed << " TESTS FAILED \n " << num_passed << " TESTS PASSED" << endl;
    } else {
        cout << "\n\nALL " << num_passed << " TESTS PASSED.\n" << endl;
    }
    // cout.close();
    // system("cat < argname.txt");
    // system("rm argname.txt");
}

bool TestObj::add_test(test_fn test, std::string func_name)
{
    tests.emplace_back(test, func_name);
    return true;
}

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

#endif