cmake_minimum_required (VERSION 3.1)

project (Cas-OFFinder
         VERSION 4.0.0
         HOMEPAGE_URL http://github.com/snugel/cas-offinder
)

set (CMAKE_CXX_STANDARD 11)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")
set (SOURCES
    bit4ops.cpp
    bulge_logic.cpp
    find_mismatches_cpu.cpp
    find_mismatches_opencl.cpp
    postprocess.cpp
)

set (HEADERS
	find_mismatches.h
    RangeIterator.h
    bulge_logic.h
    test/test_framework.h
    bit4ops.h
    postprocess.h
)

add_executable (test test/test.cpp test/test_framework.cpp ${SOURCES} ${HEADERS})
add_executable (cas-offinder search_main.cpp ${SOURCES} ${HEADERS})
add_executable (dump_file_data 
    dump_file_data.cpp
    read_twobit.cpp
    bit4ops.cpp
    read_twobit.h
)
add_executable (infinite_buffer large_pipe_buffer.cpp)
find_package(Threads REQUIRED)
target_link_libraries(infinite_buffer PRIVATE Threads::Threads)

target_compile_options(test PRIVATE -Wall -Wextra -Wno-deprecated-declarations)

target_compile_definitions(test PRIVATE ENABLE_TESTS)

link_directories(/usr/lib/x86_64-linux-gnu/)
target_link_libraries(test PUBLIC OpenCL)
target_link_libraries(cas-offinder PUBLIC OpenCL)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(cas-offinder PUBLIC OpenMP::OpenMP_CXX)
endif()
