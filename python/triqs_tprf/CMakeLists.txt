# Configure the version
configure_file(version.py.in version.py)

# All Python files. Copy them in the build dir to have a complete package for the tests.
file(GLOB_RECURSE python_sources RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.py)
foreach(file ${python_sources})
  configure_file(${file} ${file} COPYONLY)
endforeach()

# Install python files to proper location
set(PYTHON_LIB_DEST ${TRIQS_PYTHON_LIB_DEST_ROOT}/${PROJECT_NAME})
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/version.py DESTINATION ${PYTHON_LIB_DEST})
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION ${TRIQS_PYTHON_LIB_DEST_ROOT} FILES_MATCHING PATTERN "*.py")

# Add an interface target for all Python modules
add_library(${PROJECT_NAME}_python_modules INTERFACE)

# Add Python modules
foreach(module_name chi_from_gg2 freq_conv hubbard_atom lattice linalg)
  c2py_add_module(${module_name} LINK_LIBRARIES ${PROJECT_NAME}_c ${PROJECT_NAME}_warnings)
  add_library(${PROJECT_NAME}::${module_name} ALIAS ${module_name})
  add_dependencies(${PROJECT_NAME}_python_modules ${module_name})
  install(TARGETS ${module_name} DESTINATION ${PYTHON_LIB_DEST})
endforeach()
