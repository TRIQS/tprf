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
set(DOC_SOURCE ${CMAKE_CURRENT_BINARY_DIR})

# create documentation target
set(sphinx_top ${CMAKE_CURRENT_BINARY_DIR}/html/index.html)
add_custom_command(OUTPUT ${sphinx_top} DEPENDS ${doc_sources} fourier tprf_c
                   COMMAND ${TRIQS_SPHINXBUILD_EXECUTABLE} -c . -b html ${DOC_SOURCE} html)
add_custom_target(doc_sphinx ALL DEPENDS ${sphinx_top})

# Create a temporary copy of the doc that will be completed by the automatically generated one.
# NB! the conf.py.in -> conf.py is generated at the configure step
# the conf.py in the repo is for readthedocs and we must take care to
# not overwrite the generated conf.py with the one in the repo
FILE(GLOB_RECURSE all_rst_files RELATIVE ${CMAKE_SOURCE_DIR}/doc *.rst _* *.py *.png *.cpp)
SET(rst_copy_tar  ${CMAKE_CURRENT_BINARY_DIR}/rst_copy.tar)
add_custom_command (OUTPUT ${rst_copy_tar} DEPENDS ${all_rst_files}
 COMMAND cd ${CMAKE_SOURCE_DIR}/doc && tar cf ${rst_copy_tar} ${all_rst_files} && cd ${CMAKE_BINARY_DIR}/doc && cp conf.py conf.py.bak && tar xf ${rst_copy_tar} && mv conf.py.bak conf.py )
add_custom_target(rst_copy DEPENDS ${rst_copy_tar})
add_dependencies(doc_sphinx rst_copy) # we must first have copied the sources

# use c++2doc to automatically build the documentation of C++
option(DocWithC++2rst "Use c++2rst to generate C++ documentation of tprf classes and functions (requires libclang)" ON)

if (DocWithC++2rst)

  # get include paths for the triqs target
  get_property(TRIQS_INCLUDE_DIRS TARGET triqs PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  foreach(I ${TRIQS_INCLUDE_DIRS})
   set (INC  ${INC} --includes ${I})
  endforeach()

  FILE(GLOB_RECURSE all_header_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} ../c++/*.hpp ../c++/lattice/*.hpp)

  add_custom_command (OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/c++2rst.log
                      DEPENDS ${all_header_files} 
		      COMMAND c++2rst
		      ${CMAKE_CURRENT_SOURCE_DIR}/reference/doc_root.hpp
                      --output_directory ${CMAKE_CURRENT_BINARY_DIR}/cpp2doc_generated
                      --namespace tprf
                      ${INC}		    
                      --includes ${TRIQS_ROOT}/include
                      --includes ${CMAKE_SOURCE_DIR}/c++/
                      --includes ${CMAKE_SOURCE_DIR}/c++/lattice
		      --cxxflags "-std=c++17 -DFOOBAR"
                      2>&1 > c++2rst.log
		      VERBATIM)

  add_custom_target(doc_cpp2doc DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/c++2rst.log) 
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
