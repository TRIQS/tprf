# Copy h5 files to binary dir

FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

FILE(GLOB all_tgz_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.tar.gz)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_tgz_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

# List all tests
set(all_tests chi4_iw_from_tau g_wk_to_from_g_wr chi00_square_lattice bse_and_rpa_loc_vs_latt mean_field mean_field_kanamori hartree_response)

foreach(t ${all_tests})
  add_test(NAME ${t} COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${t}.py)
endforeach()

# Set the PythonPath : put the build dir first (in case there is an installed version). 
set_property(TEST ${all_tests} PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} )
