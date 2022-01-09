cmake_minimum_required (VERSION 3.1)

project (Cas-OFFinder
         VERSION 3.0.0
         HOMEPAGE_URL http://github.com/snugel/cas-offinder
)

set (CMAKE_CXX_STANDARD 11)

set (SOURCES
	main.cpp
	cas-offinder.cpp
	read_fasta.cpp
	read_twobit.cpp
	cas-offinder.cl
	oclfunctions.cpp
)

set (HEADERS
	config.h
	cas-offinder.h
	read_fasta.h
	read_twobit.h
	oclfunctions.h
	oclkernels.h
)

add_executable (cas-offinder ${SOURCES} ${HEADERS})

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-ignored-attributes")
endif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")

if (WIN32)
	message ("Windows system detected, assume no proper dirent.h is available.")
	message ("As a workaround, a custom dirent.h for Windows (by Toni Ronkko) will be used.")
	message ("For more information, please visit https://github.com/tronkko/dirent")
endif (WIN32)

find_package (OpenCL REQUIRED)
target_link_libraries(cas-offinder OpenCL::OpenCL)

find_package (OpenMP)
if (OPENMP_FOUND)
	target_link_libraries(cas-offinder OpenMP::OpenMP_CXX)
endif (OPENMP_FOUND)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")
execute_process(COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")
add_custom_command (TARGET cas-offinder PRE_BUILD COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_SOURCE_DIR}/config.h"
)
