#
# We compile all the cpp files in this directory into a shared library
#

# Find all the .cpp files in the current directory
FILE(GLOB cpp_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
FILE(GLOB cpp_sources_lattice RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} lattice/*.cpp)

# Set up the compilation of those files into a shared library
add_library(tprf_c ${cpp_sources} ${cpp_sources_lattice})

# This library should be linked to all libraries against which TRIQS is linked
target_link_libraries(tprf_c PUBLIC triqs)

# Installation path
install(TARGETS tprf_c DESTINATION lib)
