#
# We compile all the cpp files in this directory into a shared library
#

# Find all the .cpp files in the current directory
FILE(GLOB cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

# Set up the compilation of those files into a shared library
add_library(tprf_c ${cpp_sources})

# This library should be linked to all libraries against which TRIQS is linked
target_link_libraries(tprf_c ${TRIQS_LIBRARY_ALL})

# This macro sets up the rpath, to avoid any LD_LIBRARY_PATH issue
triqs_set_rpath_for_target(tprf_c)

# Installation path
install(TARGETS tprf_c DESTINATION lib)
