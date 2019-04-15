# runs a python test
# Example: add_python_test(my_script)
#   where my_script.py is the script
macro(add_python_test test)
  get_filename_component(test_name ${test} NAME_WE)
  get_filename_component(test_dir ${test} DIRECTORY)
  add_test(NAME Py_${test_name} COMMAND ${TRIQS_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${test_dir}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${test_dir})
  set_property(TEST Py_${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endmacro(add_python_test)

add_subdirectory(c++)

if(${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
endif()
