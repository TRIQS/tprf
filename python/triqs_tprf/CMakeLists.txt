# Build the python modules
add_cpp2py_module(lattice)
target_link_libraries(lattice triqs_tprf_c)

add_cpp2py_module(linalg)
target_link_libraries(linalg triqs_tprf_c)

add_cpp2py_module(chi_from_gg2)
target_link_libraries(chi_from_gg2 triqs_tprf_c)

add_cpp2py_module(freq_conv)
target_link_libraries(freq_conv triqs_tprf_c)

# Build the C++ extension module
add_cpp2py_module(hubbard_atom)
target_link_libraries(hubbard_atom triqs_tprf_c)

# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
list(FILTER python_sources EXCLUDE REGEX "_desc.py")
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python module to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/triqs_tprf)
install(TARGETS lattice linalg chi_from_gg2 freq_conv hubbard_atom DESTINATION ${PYTHON_LIB_DEST})
install(FILES ${python_sources} ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
