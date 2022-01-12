cmake_minimum_required (VERSION 3.1)

project (Cas-OFFinder
         VERSION 4.0.0
         HOMEPAGE_URL http://github.com/snugel/cas-offinder
)

set (CMAKE_CXX_STANDARD 11)

set (SOURCES
    test/test.cpp
	find_mismatches.cpp
    test/test_framework.cpp
    bulge_logic.cpp
    find_mismatches_packed.cpp
    find_mismatches_opencl.cpp
)

set (HEADERS
	find_mismatches.h
    RangeIterator.h
    bulge_logic.h
    test/test_framework.h
)

add_executable (cas-offinder ${SOURCES} ${HEADERS})
target_compile_options(cas-offinder PRIVATE -Wall -Wextra)

link_directories(/usr/lib/x86_64-linux-gnu/)
target_link_libraries(cas-offinder OpenCL)
