#!/usr/bin/env bash

set -e

# Add ccache settings (for C++ caching)
export CCACHE_CPP2=yes
export CCACHE_NOHASHDIR=1

# Common settings
BUILD_TYPE=${1:-Debug} # Default is Debug
BUILD_DIR="build-${BUILD_TYPE,,}" # Convert to lowercase to generate the directory name

# Debug options
ENABLE_NODEBUG=${ENABLE_NODEBUG:-OFF} # Can be switched via environment variable or script
ENABLE_DBGPRINT=${ENABLE_DBGPRINT:-OFF}

# Log file names
CMAKE_LOG="cmake-${BUILD_TYPE,,}.log"
BUILD_LOG="build-${BUILD_TYPE,,}.log"

# Cross-platform logical CPU count
if command -v nproc > /dev/null 2>&1; then
  logical_cpus=$(nproc)
elif command -v sysctl > /dev/null 2>&1; then
  logical_cpus=$(sysctl -n hw.logicalcpu)
else
  logical_cpus=4
fi

# Build parallelism: BUILD_JOBS env var overrides everything.
# Otherwise: take stricter of (a) core-headroom and (b) memory cap (~2GB/job).
if [ -n "${BUILD_JOBS:-}" ]; then
  total_threads="${BUILD_JOBS}"
else
  # (a) Core headroom: reserve max(2, ncpu/8) cores for IDE/browser
  headroom=$(( logical_cpus / 8 ))
  [ "$headroom" -lt 2 ] && headroom=2
  core_jobs=$(( logical_cpus - headroom ))
  [ "$core_jobs" -lt 1 ] && core_jobs=1

  # (b) Memory cap: ~2 GB per C++ compile job (Eigen-heavy templates)
  if [ -r /proc/meminfo ]; then
    mem_kb=$(awk '/^MemTotal:/ {print $2}' /proc/meminfo)
    mem_gb=$(( mem_kb / 1024 / 1024 ))
  elif command -v sysctl > /dev/null 2>&1; then
    mem_bytes=$(sysctl -n hw.memsize 2>/dev/null || echo 8589934592)
    mem_gb=$(( mem_bytes / 1024 / 1024 / 1024 ))
  else
    mem_gb=8
  fi
  mem_jobs=$(( mem_gb / 2 ))
  [ "$mem_jobs" -lt 1 ] && mem_jobs=1

  if [ "$core_jobs" -lt "$mem_jobs" ]; then
    total_threads="$core_jobs"
  else
    total_threads="$mem_jobs"
  fi
fi
export CMAKE_BUILD_PARALLEL_LEVEL="$total_threads"

echo "logical_cpus: $logical_cpus"
echo "total_threads for process: $total_threads"

# Compiler selection: honor CC/CXX env vars.
# Nix shell: use cc/c++ (Nix wrapper, Clang on macOS, GCC on Linux).
# Non-Nix macOS: Homebrew GCC if available, otherwise Apple Clang.
# Non-Nix Linux: gcc/g++.
if [ -n "${IN_NIX_SHELL:-}" ]; then
  : "${CC:=cc}"
  : "${CXX:=c++}"
elif [ "$(uname -s)" = "Darwin" ]; then
  if command -v gcc-14 > /dev/null 2>&1; then
    : "${CC:=gcc-14}"
    : "${CXX:=g++-14}"
  elif command -v gcc-13 > /dev/null 2>&1; then
    : "${CC:=gcc-13}"
    : "${CXX:=g++-13}"
  else
    # Use Apple Clang by absolute path so a Homebrew LLVM clang earlier on
    # PATH is not picked up (it can fail to build the bundled fmt version).
    : "${CC:=/usr/bin/clang}"
    : "${CXX:=/usr/bin/clang++}"
  fi
else
  : "${CC:=gcc}"
  : "${CXX:=g++}"
fi

# Homebrew libomp hint for Apple Clang (not needed in a Nix environment)
EXTRA_CMAKE_ARGS=""
if [ "$(uname -s)" = "Darwin" ] && [ -z "${IN_NIX_SHELL:-}" ]; then
  LIBOMP_PREFIX="$(brew --prefix libomp 2>/dev/null)"
  if [ -d "$LIBOMP_PREFIX" ]; then
    EXTRA_CMAKE_ARGS="-DOpenMP_ROOT=$LIBOMP_PREFIX"
  fi
fi

# Remove logs
rm -f "$CMAKE_LOG" "$BUILD_LOG"

# CMake generation (specify Ninja as the generator)
cmake \
  -G Ninja \
  -S . \
  -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DNODEBUG="$ENABLE_NODEBUG" \
  -DDBGPRINT="$ENABLE_DBGPRINT" \
  $EXTRA_CMAKE_ARGS \
  > "$CMAKE_LOG" 2>&1 || {
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo " CMake configure failed. Check $CMAKE_LOG"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    exit 1
  }

# Build
cmake \
  --build "$BUILD_DIR" \
  -j "$total_threads" \
  --verbose \
  > "$BUILD_LOG" 2>&1 || {
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo " Build failed. Check $BUILD_LOG"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    exit 1
  }

echo "####################################################"
echo " OK! Build completed successfully for $BUILD_TYPE."
echo "####################################################"
