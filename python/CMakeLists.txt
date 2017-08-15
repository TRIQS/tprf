# where will the python end up in triqs python tree ?
set(python_destination pytriqs/applications/tprf)

# make a local pytriqs copy to be able to run in the current dir without installation
triqs_prepare_local_pytriqs(${python_destination})

include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})

# Build the C++ extension module
#triqs_python_extension(tprf ${python_destination})
#target_link_libraries(tprf tprf_c)
#include_directories(${TRIQS_INCLUDE_ALL} ${CMAKE_CURRENT_SOURCE_DIR})
#triqs_set_rpath_for_target(tprf)
#install(TARGETS tprf DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
triqs_python_extension(fourier ${python_destination})
target_link_libraries(fourier tprf_c)
triqs_set_rpath_for_target(fourier)
install(TARGETS fourier DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
triqs_python_extension(lattice ${python_destination})
target_link_libraries(lattice tprf_c)
triqs_set_rpath_for_target(lattice)
install(TARGETS lattice DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
triqs_python_extension(linalg ${python_destination})
target_link_libraries(linalg tprf_c)
triqs_set_rpath_for_target(linalg)
install(TARGETS linalg DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Build the C++ extension module
triqs_python_extension(chi_from_gg2 ${python_destination})
target_link_libraries(chi_from_gg2 tprf_c)
triqs_set_rpath_for_target(chi_from_gg2)
install(TARGETS chi_from_gg2 DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})

# Install utilities
install(FILES Dummy.py DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT}/${python_destination})
