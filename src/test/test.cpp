#include "test_framework.h"
#include <cstdlib>

int main(int argc, char ** argv)
{
    int capture_stdout = 1;
    if(argc > 1){
        capture_stdout = !atoi(argv[1]);
    }

    return run_all_tests(capture_stdout);
}
