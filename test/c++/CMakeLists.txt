# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_ref_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5)
foreach(file ${all_h5_ref_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests linalg_ops bubble gf_inverse g_wk_to_from_g_wr chi0_wk_to_from_chi0_wr hubbard_atom_test mpi_and_openmp)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  add_executable(${test_name} ${test})
  target_link_libraries(${test_name} triqs_tprf_c gtest)
  add_test(${test_name} ${test_name})
  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    set_target_properties(${test_name} PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  endif()
  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    add_custom_command(
      TARGET ${test_name}
      COMMAND ${CPPCHECK_EXECUTABLE}
      --enable=warning,style,performance,portability
      --std=c++14
      --template=gcc
      --verbose
      --quiet
      ${CMAKE_CURRENT_SOURCE_DIR}/${test}.cpp
    )
  endif()
endforeach()
