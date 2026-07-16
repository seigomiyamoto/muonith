# exec/CMakeLists.txt

function(setup_target TARGET_NAME)
  message("Setting up target: ${TARGET_NAME}")
  add_executable(${TARGET_NAME}.exe ${TARGET_NAME}.cpp)

  # Include directories
  target_include_directories(${TARGET_NAME}.exe
    PRIVATE
      ${PROJECT_SOURCE_DIR}/include
  )

  # Link libmuog (libmuog aggregates all library dependencies)
  target_link_libraries(${TARGET_NAME}.exe PRIVATE libmuog)

  # Require C++20 features
  target_compile_features(${TARGET_NAME}.exe PUBLIC cxx_std_20)

  # Compile-time constant definitions
  target_compile_definitions(${TARGET_NAME}.exe PUBLIC MY_ORIGINAL_MACRO_CONST=1)

  # Per-build-type compile options
  if(CMAKE_BUILD_TYPE MATCHES Debug)
    target_compile_options(${TARGET_NAME}.exe PUBLIC -O0 -Wall -Wfatal-errors -Wundef -Wreturn-type -Wuninitialized)
    message("Debug build with -O0")
    add_definitions(-DDEBUG)
  elseif(CMAKE_BUILD_TYPE MATCHES Release)
    target_compile_options(${TARGET_NAME}.exe PUBLIC -O3 -Wall -Wfatal-errors -Wundef -Wreturn-type -Wuninitialized)
    message("Release build with -O3")
  endif()

  # Register spdlog_pch
  if(USE_PCH)
    target_precompile_headers(${TARGET_NAME}.exe PRIVATE ${PROJECT_SOURCE_DIR}/include/spdlog_pch.hpp)
  endif()

endfunction()

# Executable targets to build
# Generates penetrating-flux and dF/dR muon flux tables consumed by the pipeline.
setup_target(make_peneflux_dFdR)
# Main tool: full 3D muon-density reconstruction and parameter sweep.
setup_target(muonith)
# Depth vs. spatial resolution sweep producing detection-limit heatmaps.
setup_target(depth_reso)
# Unit test: to_blas_int() int64 -> BLAS-int conversion with overflow guards.
setup_target(test_to_blas_int)
# Guards the Grid1d/Grid2d save-load pass-through contract (tf_half_shift removal).
setup_target(test_grid1d_roundtrip)
# Guards Grid2dXYZ::load_g2zbin against v2 files in both z precisions
# (float32/float64), matching the muonith-gsi-dem Python writer contract.
setup_target(test_g2zbin_read)
# User-facing tool: rotate vectors by yaw/pitch/roll (LOCAL/GLOBAL) and show the attitude.
setup_target(rotate_vector)
