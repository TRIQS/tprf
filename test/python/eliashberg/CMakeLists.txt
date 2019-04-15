FILE(GLOB all_tgz_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.tar.gz)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_tgz_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set(PREFIX eliashberg-)

add_python_test(product_summation_vs_fft ${PREFIX})
add_python_test(eigenvalue_solver ${PREFIX})
add_python_test(previous_implementation ${PREFIX})
