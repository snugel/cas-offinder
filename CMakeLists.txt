cmake_minimum_required (VERSION 3.1)
project (Cas-OFFinder)

set (Cas-OFFinder_VERSION_MAJOR 2)
set (Cas-OFFinder_VERSION_MINOR 4)

include (CheckIncludeFiles)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")

find_package (OpenCL REQUIRED)
include_directories ("${OpenCL_INCLUDE_DIRS}")
link_directories (${OpenCL_LIBRARY})

include_directories ("${PROJECT_SOURCE_DIR}")
if (WIN32)
	message ("Windows system detected, assume no proper dirent.h is available.")
	message ("As a workaround, a custom dirent.h for Windows (by Toni Ronkko) will be used.")
	message ("For more information, please visit https://github.com/tronkko/dirent")
endif (WIN32)

find_package(OpenMP)
if (OPENMP_FOUND)
	set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_FOUND)

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_SOURCE_DIR}/config.h"
)

execute_process(COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")

set (SOURCES
	main.cpp
	cas-offinder.cpp
	read_fasta.cpp
	read_twobit.cpp
	cas-offinder.cl
)

set (HEADERS
	config.h
	cas-offinder.h
	read_fasta.h
	read_twobit.h
	oclfunctions.h
	oclkernels.h
)

add_executable(cas-offinder ${SOURCES} ${HEADERS})
if (APPLE)
	target_link_libraries(cas-offinder "-framework OpenCL")
else (APPLE)
	target_link_libraries(cas-offinder ${OpenCL_LIBRARY})
endif (APPLE)
add_custom_command(TARGET cas-offinder PRE_BUILD COMMAND ${CMAKE_COMMAND} "-DPROJECT_SOURCE_DIR=${PROJECT_SOURCE_DIR}" -P "${CMAKE_MODULE_PATH}/CopyKernels.cmake")
