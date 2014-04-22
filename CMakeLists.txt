cmake_minimum_required (VERSION 2.6)
project (CAS-OFFinder)

set (CAS-OFFinderr_VERSION_MAJOR 2)
set (CAS-OFFinderr_VERSION_MINOR 1)

set(CMAKE_BUILD_TYPE Release)

include (CheckIncludeFiles)
check_include_files (dirent.h HAVE_DIRENT_H)
check_include_files (CL/cl.hpp HAVE_CL_HPP)

set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMakeModules")
if (APPLE)
	find_library(OPENCL_LIBRARY OpenCL)
else (APPLE)
	find_package (OpenCL REQUIRED)
	include_directories ("${OPENCL_INCLUDE_DIR}")
endif (APPLE)
include_directories ("${PROJECT_BINARY_DIR}")
if (NOT HAVE_DIRENT_H)
	message ("Could not find dirent.h on your system.")
	message ("For workaround, using a custom dirent.h written for MS system (by Toni Ronkko).")
	message ("More information, visit http://softagalleria.net/dirent.php")
	include_directories ("${PROJECT_BINARY_DIR}/ext")
endif (NOT HAVE_DIRENT_H)

configure_file (
	"${PROJECT_SOURCE_DIR}/config.h.in"
	"${PROJECT_BINARY_DIR}/config.h"
)

add_executable(cas-offinder main.cpp)
target_link_libraries(cas-offinder ${OPENCL_LIBRARY})
