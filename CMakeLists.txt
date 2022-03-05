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
add_library(genome_base_lib STATIC
    src/bit4ops.cpp
    src/read_genome.cpp
    src/bidirect_search.cpp
    src/read_fasta.cpp
    src/read_2bit.cpp
    src/parse_input.cpp
    src/blockify.cpp
    src/bidirect_search.cpp
    src/bulge_logic.cpp
    ${HEADERS}
    )
target_link_libraries(genome_base_lib PUBLIC Threads::Threads)

add_library(genome_search_lib STATIC
    src/async_search.cpp
    src/search.cpp
    src/oclfunctions.cpp
    ${HEADERS}
    )
target_link_libraries(genome_search_lib PUBLIC OpenCL)
target_link_libraries(genome_search_lib PUBLIC genome_base_lib)

add_library(genome_nocl_search_lib STATIC
    src/async_search.cpp
    src/search_nodep.cpp
    ${HEADERS}
    )
target_link_libraries(genome_nocl_search_lib PUBLIC genome_base_lib)

add_executable (test
    src/test/test.cpp
    src/test/test_framework.cpp
    src/test/test_framework.h
    src/test/bit4ops.cpp
    src/test/search.cpp
    src/test/read_file.cpp
    src/test/parse_input.cpp
    src/test/blockify.cpp
    src/test/bulge_logic.cpp
    src/test/async_search.cpp
    src/test/bidirect_search.cpp
    ${HEADERS}
)
target_link_libraries(test PUBLIC genome_nocl_search_lib )

add_executable (cas-offinder-2
    src/main_v2.cpp
    ${HEADERS}
)
target_link_libraries(cas-offinder-2 PUBLIC genome_search_lib )

add_executable (cas-offinder-2-nocl
    src/main_v2.cpp
    ${HEADERS}
)
target_link_libraries(cas-offinder-2-nocl PUBLIC genome_nocl_search_lib )
add_executable (cas-offinder-3-nocl
 src/main_v3.cpp
   ${HEADERS}
   )
target_link_libraries(cas-offinder-3-nocl PUBLIC genome_nocl_search_lib )
add_executable (cas-offinder-3
 src/main_v3.cpp
   ${HEADERS}
   )
target_link_libraries(cas-offinder-3 PUBLIC genome_search_lib )

if(WIN32)
else()
    target_compile_options(genome_base_lib PUBLIC -Wall -Wextra)
endif()
