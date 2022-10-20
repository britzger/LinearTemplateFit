# Works with 3.11 and tested through 3.19
cmake_minimum_required(VERSION 3.11...3.19)

# Project name and a few useful settings. Other commands can pick up the results
project(
  LTF
  VERSION 1.0
  DESCRIPTION "Linear Template Fit"
  LANGUAGES C CXX)

find_package(ROOT 6.26 CONFIG REQUIRED)

# Sets up global settings
#include("${ROOT_USE_FILE}")

# This is required on if there is more than one flag (like on macOS)
separate_arguments(ROOT_EXE_LINKER_FLAGS)

SET(CMAKE_CXX_STANDARD 17)

SET(CMAKE_INSTALL_INCLUDEDIR ${CMAKE_BINARY_DIR}/include)
SET(CMAKE_INSTALL_LIBDIR ${CMAKE_BINARY_DIR}/lib)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

set(HEADER_OUTPUT_PATH ${CMAKE_BINARY_DIR}/include)

include_directories(${ROOT_INCLUDE_DIRS})
include_directories("${CMAKE_SOURCE_DIR}/LTF/inc")

add_subdirectory(LTF)

add_subdirectory(examples)

get_property(__allHeaders GLOBAL PROPERTY ROOT_HEADER_TARGETS)
get_property(__allBuiltins GLOBAL PROPERTY ROOT_BUILTIN_TARGETS)
add_custom_target(move_headers ALL DEPENDS ${__allHeaders} ${__allBuiltins})
