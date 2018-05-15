# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
 
# List all tests
set(all_tests chi4_iw_from_tau gk_to_from_gr chi00_square_lattice)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath : put the build dir first (in case there is an installed version). 
set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
