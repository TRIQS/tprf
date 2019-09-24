# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_tgz_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5 *.tar.gz)
foreach(file ${all_h5_tgz_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests chi4_iw_from_tau g_wk_to_from_g_wr_py chi00_square_lattice bse_and_rpa_loc_vs_latt mean_field mean_field_kanamori hartree_response 1d_hubbard_hf_rpa 1d_hubbard_hf_spin_rot_inv 1d_hubbard_hf_rpa_2site_AFM eliashberg compare_general_rpa_to_matrix_rpa interaction_tensor_charge_spin_factorization)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endforeach()
