# start configuration
cmake_minimum_required(VERSION 2.8)
project(tprf CXX)
set(CMAKE_BUILD_TYPE Release)
enable_testing()

# We use shared libraries
option(BUILD_SHARED_LIBS "Build shared libraries" ON)

# Load TRIQS and Cpp2Py
find_package(TRIQS 2.0 REQUIRED)
find_package(Cpp2Py 1.5 REQUIRED)

# Load TRIQS, including all predefined variables from TRIQS installation
# First append triqs installed files to the cmake load path
list(APPEND CMAKE_MODULE_PATH ${TRIQS_PATH}/share/triqs/cmake)

add_definitions(${TRIQS_CXX_DEFINITIONS})
include_directories(${TRIQS_INCLUDE_ALL})

# Default Install directory to TRIQS_ROOT if not given. Checks an absolute name is given.
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT OR (NOT IS_ABSOLUTE ${CMAKE_INSTALL_PREFIX}))
  message(STATUS " No install prefix given (or invalid). Defaulting to TRIQS_ROOT")
  set(CMAKE_INSTALL_PREFIX ${TRIQS_ROOT} CACHE PATH "default install path" FORCE)
endif()
message(STATUS "-------- CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX} -------------")

# All PRIVATE common options.
# The std for all targets
add_compile_options(-std=c++17)

# Common to all targets (ctint, python, tests).
include_directories(${CMAKE_SOURCE_DIR})

# Get the Git Hash
find_package(Git)

execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
 WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}"
 OUTPUT_VARIABLE   TPRF_GIT_HASH
 OUTPUT_STRIP_TRAILING_WHITESPACE
)
message(STATUS "TPRF git hash : ${TPRF_GIT_HASH}")

# To be fixed!!
triqs_set_rpath_for_all_targets()

# Compile C++ code
add_subdirectory(c++)

# Python interface (if is compiled with Python support).
# Remove this for a pure C++ code.
if (${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)

 # Build Sphinx based docs
 option(BUILD_DOC "Build documentation" OFF)
 if(${BUILD_DOC})
   if(NOT ${TRIQS_WITH_DOCUMENTATION})
     message(FATAL_ERROR "Error: TRIQS library has not been compiled with its documentation")
   endif()
   add_subdirectory(doc)
 endif(${BUILD_DOC})
endif()

# Add tests
enable_testing()
add_subdirectory(test)
