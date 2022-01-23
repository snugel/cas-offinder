cmake_minimum_required (VERSION 3.1)

project (Cas-OFFinder
         VERSION 4.0.0
         HOMEPAGE_URL http://github.com/snugel/cas-offinder
)
include_directories(
    ${PROJECT_SOURCE_DIR}/src
)


find_package(OpenMP)
find_package(Threads REQUIRED)
include_directories(.)

set (CMAKE_CXX_STANDARD 14)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")

set (HEADERS
	src/find_mismatches.h
    src/RangeIterator.h
    src/bulge_logic.h
    src/test/test_framework.h
    src/bit4ops.h
    src/chromloc.h
    src/postprocess.h
    src/read_twobit.h
    src/read_fasta.h
    src/read_folder.h
    src/search_genome.h
)

add_executable (test 
    src/bit4ops.cpp
    src/chromloc.cpp
    src/bulge_logic.cpp
    src/find_mismatches_cpu.cpp
    src/find_mismatches_opencl.cpp
    src/find_mismatches_utils.cpp
    src/postprocess.cpp
    src/test/test.cpp 
    src/test/test_framework.cpp 
    src/test/postprocess.cpp
    src/test/find_mismatches.cpp
    src/test/chromloc.cpp
    src/test/bit4ops.cpp
    src/test/bulge_logic.cpp
    ${SOURCES} ${HEADERS}
)

add_executable (cas-offinder 
    src/search_main.cpp 
    src/bit4ops.cpp
    src/chromloc.cpp
    src/find_mismatches_opencl.cpp
    src/search_genome.cpp
    src/find_mismatches_utils.cpp
    src/postprocess.cpp
    ${HEADERS}
)

add_executable (dump_file_data 
    src/dump_file_data.cpp
    src/bit4ops.cpp
    src/chromloc.cpp
    src/read_twobit.cpp
    src/read_fasta.cpp
    src/read_folder.cpp
    ${HEADERS}
)
add_executable (infinite_buffer 
    src/large_pipe_buffer.cpp
)

target_link_libraries(infinite_buffer PRIVATE Threads::Threads)

target_compile_options(test PRIVATE -Wall -Wextra)

target_link_libraries(test PRIVATE Threads::Threads)
target_link_libraries(test PUBLIC OpenCL)
target_link_libraries(cas-offinder PUBLIC OpenCL)

link_directories(/usr/lib/x86_64-linux-gnu/)

if(OpenMP_CXX_FOUND)
    target_link_libraries(test PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(cas-offinder PUBLIC OpenMP::OpenMP_CXX)
endif()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")
execute_process(COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")
add_custom_command (TARGET cas-offinder PRE_BUILD COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")
