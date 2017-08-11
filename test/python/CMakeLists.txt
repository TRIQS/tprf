
# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
 
triqs_add_python_test(g_iw_from_tau)
triqs_add_python_test(chi4_iw_from_tau)

