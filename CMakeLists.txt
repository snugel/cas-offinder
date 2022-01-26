cmake_minimum_required (VERSION 3.9)

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

#set(CMAKE_INTERPROCEDURAL_OPTIMIZATION TRUE)
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -s")

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")
#execute_process(COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")
add_custom_command (OUTPUT ${PROJECT_SOURCE_DIR}/oclkernels.h DEPENDS ${PROJECT_SOURCE_DIR}/src/find_mismatches_packed.cl COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")

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
    src/read_genome.h
    src/chunkify_data.h
    ${PROJECT_SOURCE_DIR}/oclkernels.h
)

add_library(genome_search_lib STATIC
    src/bit4ops.cpp
    src/chromloc.cpp
    src/bulge_logic.cpp
    src/find_mismatches_opencl.cpp
    src/find_mismatches_utils.cpp
    src/postprocess.cpp
    src/read_genome.cpp
    src/search_genome.cpp
    src/find_mismatches_cpu.cpp
    src/chunkify_data.cpp
    ${HEADERS}
)

add_executable (genome_search
    src/search_main.cpp
    ${HEADERS}
)
add_executable (cas-offinder-2
    src/cas-offinder-2.cpp
    ${HEADERS}
)
add_executable (cas-offinder-3
    src/cas-offinder-3.cpp
    ${HEADERS}
)


add_executable (test 
    src/test/test.cpp 
    src/test/test_framework.cpp 
    src/test/postprocess.cpp
    src/test/find_mismatches.cpp
    src/test/chromloc.cpp
    src/test/bit4ops.cpp
    src/test/bulge_logic.cpp
    src/test/chunkify_data.cpp
    ${SOURCES} ${HEADERS}
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

target_compile_options(test PRIVATE -Wall -Wextra)

target_link_libraries(genome_search_lib PUBLIC Threads::Threads)
target_link_libraries(genome_search_lib PUBLIC OpenCL)

target_link_libraries(test PUBLIC genome_search_lib)
target_link_libraries(genome_search PUBLIC genome_search_lib)
target_link_libraries(cas-offinder-2 PUBLIC genome_search_lib)
target_link_libraries(cas-offinder-3 PUBLIC genome_search_lib)


link_directories(/usr/lib/x86_64-linux-gnu/)
