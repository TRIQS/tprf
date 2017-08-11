# Path to the Sphinx inventory file of the library
set(TRIQS_INVENTORY_DIR "${TRIQS_PATH}/share/doc/triqs" CACHE STRING
    "Path or URL to a directory containing a Sphinx inventory file of the TRIQS library")
mark_as_advanced(TRIQS_INVENTORY_DIR)

# Path to the Sphinx extensions provided by TRIQS
set(TRIQS_SPHINXEXT_PATH "${TRIQS_PATH}/share/triqs/sphinxext")

# generate the conf.py
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in ${CMAKE_CURRENT_BINARY_DIR}/conf.py)

# all rst files of the documentation
file(GLOB_RECURSE doc_sources *.rst)

# create documentation target
set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${doc_sources} py_copy fourier tprf_c
                   COMMAND ${CMAKE_BINARY_DIR}/build_pytriqs ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b html ${CMAKE_CURRENT_SOURCE_DIR} html)
add_custom_target(doc_sphinx ALL DEPENDS ${sphinx_top})

# use c++2doc to automatically build the documentation of C++
option(DocWithCpp2doc "Use c++2doc to generate C++ documentation of tprf classes and functions (requires libclang)" ON)
if (DocWithCpp2doc)

FUNCTION(PREPEND_LIST var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND_LIST)

PREPEND_LIST(CPP2DOC_INCLUDE_FLAGS "--includes=" ${TRIQS_INCLUDE_ALL})

 FILE(GLOB_RECURSE all_header_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ../c++/*.hpp ../c++/*.mhpp)

 message(STATUS "all_header_files")
 message(STATUS ${all_header_files})
 
add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/cpp2py.log DEPENDS ${all_header_files} 
                    COMMAND export PYTHONPATH=$$PYTHONPATH:${TRIQS_PATH}/share/triqs/cpp2doc/mako &&
		    ${PYTHON_INTERPRETER} 
                    ${TRIQS_PATH}/bin/c++2doc.py 
                    -N c++
                    ${CMAKE_CURRENT_SOURCE_DIR}/reference/doc_root.hpp 
                    --output_directory=${CMAKE_CURRENT_BINARY_DIR}/cpp2doc_generated
                    --includes=${CMAKE_SOURCE_DIR}/ 
		    ${CPP2DOC_INCLUDE_FLAGS}
                    --mako_location=${TRIQS_PATH}/share/triqs
		    ${all_header_files}
                    #2>&1 > cpp2py.log
                    ) 
add_custom_target(doc_cpp2doc DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/cpp2py.log) 
add_dependencies(doc_cpp2doc tprf_c) #
add_dependencies(doc_sphinx doc_cpp2doc) # we must have installed the sources before building
endif()

# install
install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html/ COMPONENT documentation DESTINATION share/doc/tprf
 FILES_MATCHING 
 PATTERN "*.html" 
 PATTERN "*.png"
 PATTERN "*.js"
 PATTERN "_*"
 PATTERN "*.png" 
 PATTERN "*.jpg" 
 PATTERN "*.gif" 
 PATTERN "*.xsl" 
 PATTERN "*.css"
 PATTERN "*.pdf"
 PATTERN "*.py"
 PATTERN "*.txt"
 PATTERN "*.inv"
 PATTERN "*.bib"
 )
