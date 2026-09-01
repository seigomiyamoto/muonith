# ./CMakeLists.txt

cmake_minimum_required(VERSION 3.22.2)

project(muonith CXX)

option(ENABLE_GC_SECTIONS "Enable --gc-sections and report removed sections" ON)
option(ENABLE_UNUSED_WARNINGS "Enable -Wunused-function warnings" ON)
option(NODEBUG "Define NODEBUG" OFF)
option(DBGPRINT "Define DBGPRINT" OFF)

if(ENABLE_UNUSED_WARNINGS)
  add_compile_options(-Wunused-function)
endif()

if(APPLE)
  add_link_options(-Wl,-dead_strip)
  # Allow unresolved symbols in shared libraries (resolved at load time).
  # Required because shared libs have circular dependencies (e.g. nagamu_util <-> nagamu_detector).
  # Linux ld permits this by default; macOS ld does not.
  add_link_options(-Wl,-undefined,dynamic_lookup)
elseif(ENABLE_GC_SECTIONS)
  add_compile_options(-ffunction-sections -fdata-sections)
  add_link_options(-Wl,--gc-sections)
endif()

if(NODEBUG)
  add_compile_definitions(NODEBUG)
endif()
if(DBGPRINT)
  add_compile_definitions(DBGPRINT)
endif()

# Pin floating-point contraction to keep numerical results identical across
# compilers and targets: GCC defaults to -ffp-contract=fast and Clang to "on",
# so FMA fusion could silently change last digits on arm64 or with -march.
add_compile_options(-ffp-contract=off)

# Generate compile_commands.json (per-file compile information)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set RPATH automatically from the link-time paths
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(APPLE)
  set(CMAKE_MACOSX_RPATH ON)
endif()

add_compile_definitions(
  SPDLOG_SHARED_LIB
  SPDLOG_FMT_EXTERNAL
  EIGEN_MPL2_ONLY
)

# Add SPDLOG_ACTIVE_LEVEL to the compile definitions for Release builds
if(CMAKE_BUILD_TYPE MATCHES Release)
  add_compile_definitions(SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_TRACE)
  message(STATUS "Setting SPDLOG_ACTIVE_LEVEL for Release build")
endif()

# C++ standard settings
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
option(USE_PCH "Use precompiled headers" OFF)

# External library lookup and other settings
find_package(OpenMP REQUIRED)
# BLAS/LAPACK backend (platform auto-detection)
if(APPLE)
  find_library(ACCELERATE_FRAMEWORK Accelerate REQUIRED)
  set(BLAS_LINK_LIBRARIES ${ACCELERATE_FRAMEWORK})
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  find_library(LAPACKE_LIB lapacke REQUIRED)
  set(BLAS_LINK_LIBRARIES BLAS::BLAS LAPACK::LAPACK ${LAPACKE_LIB})
endif()

# --- FetchContent: vendored dependencies ---
# Ensure all static libraries are built with -fPIC (required when
# linking into shared libraries such as spdlog).
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
include(FetchContent)

# Eigen 3.4.0 (header-only)
FetchContent_Declare(Eigen3
  GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
  GIT_TAG 3.4.0
  GIT_SHALLOW TRUE)
set(EIGEN_BUILD_DOC OFF CACHE BOOL "" FORCE)
set(BUILD_TESTING OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(Eigen3)

# nlohmann/json 3.11.3 (header-only)
FetchContent_Declare(json
  GIT_REPOSITORY https://github.com/nlohmann/json.git
  GIT_TAG v3.11.3
  GIT_SHALLOW TRUE)
set(JSON_BuildTests OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(json)

# fmt 11.2.0
FetchContent_Declare(fmt
  GIT_REPOSITORY https://github.com/fmtlib/fmt.git
  GIT_TAG 11.2.0
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(fmt)

# spdlog 1.15.3 (uses external fmt)
set(SPDLOG_FMT_EXTERNAL ON CACHE BOOL "" FORCE)
set(SPDLOG_BUILD_SHARED ON CACHE BOOL "" FORCE)
FetchContent_Declare(spdlog
  GIT_REPOSITORY https://github.com/gabime/spdlog.git
  GIT_TAG v1.15.3
  GIT_SHALLOW TRUE)
FetchContent_MakeAvailable(spdlog)

# Sub-directories
add_subdirectory(src)
add_subdirectory(exec)
