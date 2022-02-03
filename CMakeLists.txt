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

set (CMAKE_CXX_STANDARD 17)



set (HEADERS
    headers/bit4ops.h
    headers/search.h
    headers/RangeIterator.h
    headers/ceildiv.h
)
set(SOURCES
    src/bit4ops.cpp
    src/search.cpp
)

add_executable (test 
    src/test/test.cpp 
    src/test/test_framework.cpp 
    src/test/bit4ops.cpp
    src/test/search.cpp
    ${SOURCES} ${HEADERS}
)

target_compile_options(test PRIVATE -Wall -Wextra)

