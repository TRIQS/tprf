configure_file(triqs_tprf-config.cmake.in triqs_tprf-config.cmake @ONLY)
configure_file(triqs_tprf-config-version.cmake.in triqs_tprf-config-version.cmake @ONLY)
install(
  FILES
    ${CMAKE_CURRENT_BINARY_DIR}/triqs_tprf-config.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/triqs_tprf-config-version.cmake
  DESTINATION lib/cmake/triqs_tprf
)

install(EXPORT triqs_tprf-targets NAMESPACE triqs_tprf:: DESTINATION lib/cmake/triqs_tprf)
