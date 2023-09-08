# Copy h5 files to binary dir
file(GLOB_RECURSE all_h5_tgz_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.ref.h5 *.tar.gz)
foreach(file ${all_h5_tgz_files})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# List of all tests
set(all_tests
  chi4_iw_from_tau
  g_wk_to_from_g_wr_py
  lattice_utility
  V_orb_to_band
  gf
  gw
  gw_hubbard_dimer
  gw_hubbard_dimer_dlr
  gw_hubbard_dimer_matrix
  gw_hubbard_dimer_scGW_G0W0
  gw_hubbard_atom_hf
  gw_hubbard_dimer_hf
  gw_singlekpoint_exactresult
  gw_singlekpoint_twoband
  gw_compare_spectralRep_and_direct
  gw_fk
  g0w_separate_kpoints
  fitdlr_hubbard_atom
  chi00_square_lattice
  chi00_square_lattice_fk
  chi00_square_lattice_dlr
  dynamical_screened_interaction_w
  dynamical_screened_interaction_w_from_generalized_susceptibility
  chi0_save_memory
  chi0_nr_at_specific_w
  chi0_nk_at_specific_w
  solve_lattice_bse_at_specific_w
  bse_and_rpa_loc_vs_latt
  mean_field
  mean_field_kanamori
  hartree_response
  1d_hubbard_hf_rpa
  1d_hubbard_hf_spin_rot_inv
  1d_hubbard_hf_rpa_2site_AFM
  compare_general_rpa_to_matrix_rpa
  interaction_tensor_charge_spin_factorization
  symmetrize_gf
  add_fake_bosonic_mesh  
  eliashberg/preprocessing_gamma
  eliashberg/gamma_creation
  eliashberg/solve_eliashberg_functionality
  eliashberg/semi_random_initial_delta
  eliashberg/product_summation_vs_fft
  eliashberg/eigenvalue_solver
  eliashberg/regression_test_one_band
  eliashberg/regression_test_two_band
  eliashberg/fft_product_constant_vs_full
  eliashberg/symmetrize_delta
  eliashberg/compare_dlr_and_direct
  eliashberg/dlr_eliashberg_solver
)

foreach(test ${all_tests})
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME Py_${test} COMMAND ${TRIQS_PYTHON_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${PROJECT_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endforeach()
