# where will the python end up in triqs python tree ?
set(python_destination triqs_tprf)

include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})

# Build the C++ extension module
#triqs_python_extension(tprf ${python_destination})
#target_link_libraries(tprf tprf_c)
#include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
#triqs_set_rpath_for_target(tprf)
#install(TARGETS tprf DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
add_cpp2py_module(fourier)
target_link_libraries(fourier tprf_c)
triqs_set_rpath_for_target(fourier)
install(TARGETS fourier DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
add_cpp2py_module(lattice)
target_link_libraries(lattice tprf_c)
triqs_set_rpath_for_target(lattice)
install(TARGETS lattice DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
add_cpp2py_module(linalg)
target_link_libraries(linalg tprf_c)
triqs_set_rpath_for_target(linalg)
install(TARGETS linalg DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
add_cpp2py_module(chi_from_gg2)
target_link_libraries(chi_from_gg2 tprf_c)
triqs_set_rpath_for_target(chi_from_gg2)
install(TARGETS chi_from_gg2 DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
add_cpp2py_module(freq_conv)
target_link_libraries(freq_conv tprf_c)
triqs_set_rpath_for_target(freq_conv)
install(TARGETS freq_conv DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
add_cpp2py_module(hubbard_atom)
target_link_libraries(hubbard_atom tprf_c)
triqs_set_rpath_for_target(hubbard_atom)
install(TARGETS hubbard_atom DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Install utilities
install(FILES logo.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES tight_binding.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES wannier90.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES lattice_utils.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES ParameterCollection.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES analytic_hubbard_atom.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES Dummy.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

install(FILES __init__.py DESTINATION ${CPP2PY_PYTHON_LIB_DEST_ROOT}/${python_destination})

