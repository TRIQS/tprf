# Copy h5 files to binary dir
FILE(GLOB all_h5_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.h5)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${all_h5_files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

include_directories(${CMAKE_SOURCE_DIR}/c++)

#triqs_add_all_cpp_tests_in_current_dir()

# Linalg Test
#add_executable(linalg_ops ${CMAKE_CURRENT_SOURCE_DIR}/linalg_ops.cpp)
#target_link_libraries (linalg_ops ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST} tprf_c)
#add_test(linalg_ops ${CMAKE_CURRENT_BINARY_DIR}/linalg_ops)

set(SIMPLE_TESTS linalg_ops bubble gf_inverse gk_to_from_gr)

foreach(t ${SIMPLE_TESTS})
  add_executable(${t} ${CMAKE_CURRENT_SOURCE_DIR}/${t}.cpp)
  target_link_libraries (${t} ${TRIQS_LIBRARY_ALL} ${TRIQS_LIBRARY_GTEST} tprf_c)
  add_test(${t} ${CMAKE_CURRENT_BINARY_DIR}/${t})
endforeach(t)
