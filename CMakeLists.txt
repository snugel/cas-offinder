cmake_minimum_required (VERSION 3.1)

project (Cas-OFFinder
         VERSION 4.0.0
         HOMEPAGE_URL http://github.com/snugel/cas-offinder
)

find_package(OpenMP)
find_package(Threads REQUIRED)
include_directories(.)

set (CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "-O3")

set (HEADERS
	find_mismatches.h
    RangeIterator.h
    bulge_logic.h
    test/test_framework.h
    bit4ops.h
    postprocess.h
    read_twobit.h
    read_fasta.h
    read_folder.h
)

add_executable (test 
    test/test.cpp 
    test/test_framework.cpp 
    bit4ops.cpp
    bulge_logic.cpp
    find_mismatches_cpu.cpp
    find_mismatches_opencl.cpp
    postprocess.cpp
    test/bit4ops.cpp
    test/bulge_logic.cpp
    test/find_mismatches.cpp
    test/postprocess.cpp
    ${SOURCES} ${HEADERS}
)

add_executable (cas-offinder 
    search_main.cpp 
    bit4ops.cpp
    find_mismatches_cpu.cpp
    find_mismatches_opencl.cpp
    ${HEADERS}
)

add_executable (dump_file_data 
    dump_file_data.cpp
    bit4ops.cpp
    read_twobit.cpp
    read_fasta.cpp
    read_folder.cpp
    ${HEADERS}
)
add_executable (infinite_buffer 
    large_pipe_buffer.cpp
)

target_link_libraries(infinite_buffer PRIVATE Threads::Threads)

target_compile_options(test PRIVATE -Wall -Wextra -Wno-deprecated-declarations)

target_link_libraries(test PRIVATE Threads::Threads)
target_link_libraries(test PUBLIC OpenCL)
target_link_libraries(cas-offinder PUBLIC OpenCL)

link_directories(/usr/lib/x86_64-linux-gnu/)

if(OpenMP_CXX_FOUND)
    target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(cas-offinder PUBLIC OpenMP::OpenMP_CXX)
endif()
