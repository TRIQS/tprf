file(GLOB_RECURSE sources *.cpp)
add_library(triqs_tprf_c ${sources})

# Link against triqs and enable warnings
target_link_libraries(triqs_tprf_c PUBLIC triqs PRIVATE project_warnings)

# Configure compilation
target_compile_options(triqs_tprf_c PUBLIC -fPIC -fopenmp)
target_link_libraries(triqs_tprf_c PUBLIC "-fopenmp") # We are using openmp
target_include_directories(triqs_tprf_c PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/c++>)
target_compile_definitions(triqs_tprf_c PUBLIC
  				TPRF_GIT_HASH=${PROJECT_GIT_HASH}
				TRIQS_GIT_HASH=${TRIQS_GIT_HASH}
				$<$<CONFIG:Debug>:TPRF_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_DEBUG>
				$<$<CONFIG:Debug>:TRIQS_ARRAYS_ENFORCE_BOUNDCHECK>
			  )

# Install library and headers
install(TARGETS triqs_tprf_c EXPORT triqs_tprf-targets DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} DESTINATION include FILES_MATCHING PATTERN "*.hpp" PATTERN "*.h")


# ========= Static Analyzer Checks ==========

option(ANALYZE_SOURCES OFF "Run static analyzer checks if found (clang-tidy, cppcheck)")
if(ANALYZE_SOURCES)

  # Locate static analyzer tools
  find_program(CPPCHECK_EXECUTABLE NAMES "cppcheck" PATHS ENV PATH)
  find_program(CLANG_TIDY_EXECUTABLE NAMES "clang-tidy" PATHS ENV PATH)

  # Run clang-tidy if found
  if(CLANG_TIDY_EXECUTABLE)
    message(STATUS "clang-tidy found: ${CLANG_TIDY_EXECUTABLE}")
    set_target_properties(triqs_tprf_c PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXECUTABLE}")
  else()
    message(STATUS "clang-tidy not found in $PATH. Please consider installing clang-tidy for additional checks!")
  endif()

  # Run cppcheck if found
  if(CPPCHECK_EXECUTABLE)
    message(STATUS "cppcheck found: ${CPPCHECK_EXECUTABLE}")
    add_custom_command(
	TARGET triqs_tprf_c
	COMMAND ${CPPCHECK_EXECUTABLE}
	  --enable=warning,style,performance,portability
	  --std=c++17
	  --template=gcc
	  --verbose
	  --quiet
          ${sources}
    )
  else()
    message(STATUS "cppcheck not found in $PATH. Please consider installing cppcheck for additional checks!")
  endif()

endif()


# ========= Dynamic Analyzer Checks ==========

option(ASAN OFF "Compile library and executables with LLVM Address Sanitizer")
option(UBSAN OFF "Compile library and executables with LLVM Undefined Behavior Sanitizer")

if(ASAN)
  find_package(sanitizer REQUIRED "asan")
  target_link_libraries(triqs_tprf_c PUBLIC asan)
  install(TARGETS asan EXPORT triqs_tprf-targets)
endif()
if(UBSAN)
  find_package(sanitizer REQUIRED "ubsan")
  target_link_libraries(triqs_tprf_c PUBLIC ubsan)
  install(TARGETS ubsan EXPORT triqs_tprf-targets)
endif()
