# runs a python test
# Example: add_python_test(my_script)
#   where my_script.py is the script
macro(add_python_test testname)

    add_test(NAME "${ARGN}${testname}" COMMAND python ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.py)
    set_property(TEST "${ARGN}${testname}" PROPERTY ENVIRONMENT PYTHONPATH=${CMAKE_BINARY_DIR}/python:$ENV{PYTHONPATH})

endmacro(add_python_test)

add_subdirectory(c++)

if(${TRIQS_WITH_PYTHON_SUPPORT})
 add_subdirectory(python)
endif()
