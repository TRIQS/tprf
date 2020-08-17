# runs a python test
# Example: add_python_test(my_script)
#   where my_script.py is the script
macro(add_python_test test)
  get_filename_component(test_name ${test} NAME_WE)
  add_test(NAME Py_${ARGN}${test_name} COMMAND ${TRIQS_PYTHON_INTERPRETER} ${CMAKE_CURRENT_SOURCE_DIR}/${test_name}.py WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
  set_property(TEST Py_${ARGN}${test_name} APPEND PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH} ${SANITIZER_RT_PRELOAD})
endmacro(add_python_test)

add_subdirectory(c++)

if(PythonSupport)
 add_subdirectory(python)
endif()
