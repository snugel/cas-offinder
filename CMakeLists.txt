cmake_minimum_required (VERSION 3.9)

project (Cas-OFFinder
         VERSION 4.0.0
         HOMEPAGE_URL http://github.com/snugel/cas-offinder
)
include_directories(
    ${PROJECT_SOURCE_DIR}/src
)


find_package(Threads REQUIRED)
include_directories(.)
include_directories(headers/)

set (CMAKE_CXX_STANDARD 17)

if(WIN32)
    link_directories(${PROJECT_SOURCE_DIR}/bin)
else()
link_directories(/usr/lib/x86_64-linux-gnu/)
endif()

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")
#execute_process(COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")
add_custom_command (OUTPUT ${PROJECT_SOURCE_DIR}/oclkernels.h DEPENDS ${PROJECT_SOURCE_DIR}/src/find_mismatches_packed.cl COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")


set (HEADERS
    headers/bit4ops.h
    headers/search.h
    headers/RangeIterator.h
    headers/ceildiv.h
    headers/oclfunctions.h
    headers/read_genome.h
    headers/parse_input.h
    headers/blockify.h
    headers/async_search.h
    headers/bidirect_search.h
    headers/bulge_logic.h
    headers/opencl_executor.h
    headers/channel.h
    oclkernels.h
)
set(SOURCES
    src/bit4ops.cpp
    src/search.cpp
    src/oclfunctions.cpp
    src/read_genome.cpp
    src/read_fasta.cpp
    src/read_2bit.cpp
    src/parse_input.cpp
    src/blockify.cpp
    src/async_search.cpp
    src/bidirect_search.cpp
    src/bulge_logic.cpp
)

set(NODEP_SOURCES
    src/bit4ops.cpp
    src/search_nodep.cpp
    src/read_genome.cpp
    src/read_fasta.cpp
    src/read_2bit.cpp
    src/parse_input.cpp
    src/blockify.cpp
    src/async_search.cpp
    src/bidirect_search.cpp
    src/bulge_logic.cpp
)

add_executable (test
    test/test.cpp
    test/test_framework.cpp
    test/test_framework.h
    test/bit4ops.cpp
    test/search.cpp
    test/read_file.cpp
    test/parse_input.cpp
    test/blockify.cpp
    test/bulge_logic.cpp
    test/async_search.cpp
    test/bidirect_search.cpp
    ${NODEP_SOURCES} ${HEADERS}
)
add_executable (cas-offinder-2
    src/main_v2.cpp
    ${SOURCES} ${HEADERS}
)
add_executable (cas-offinder-2-nocl
    src/main_v2.cpp
    ${NODEP_SOURCES} ${HEADERS}
)

target_link_libraries(cas-offinder-2 PUBLIC OpenCL)
target_link_libraries(test PUBLIC OpenCL)
target_link_libraries(cas-offinder-2 PUBLIC Threads::Threads)
target_link_libraries(test PUBLIC Threads::Threads)
target_link_libraries(cas-offinder-2-nocl PUBLIC Threads::Threads)

target_compile_options(test PRIVATE -Wall -Wextra)

