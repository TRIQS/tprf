add_subdirectory(cmake)

if(NOT CMAKE_INSTALL_PREFIX STREQUAL TRIQS_ROOT)

  configure_file(triqs_tprf.modulefile.in triqs_tprf.modulefile @ONLY)
  configure_file(triqs_tprfvars.sh.in triqs_tprfvars.sh @ONLY)

  install(
    FILES
      ${CMAKE_CURRENT_BINARY_DIR}/triqs_tprf.modulefile
      ${CMAKE_CURRENT_BINARY_DIR}/triqs_tprfvars.sh
    DESTINATION share
  )
 
  message(STATUS "***************************************************************")
  message(STATUS "* Custom install Location. Use: 			         ")
  message(STATUS "*                                                              ")
  message(STATUS "*   source ${CMAKE_INSTALL_PREFIX}/share/triqs_tprfvars.sh      ")
  message(STATUS "*                                                              ") 
  message(STATUS "* to set up the environment variables                          ") 
  message(STATUS "***************************************************************")

endif()
