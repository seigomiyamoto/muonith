#!/usr/bin/env bash
# basan.sh -- AddressSanitizer debug build.
# Based on bdebug.sh; injects ASan flags via CMAKE_CXX_FLAGS / CMAKE_C_FLAGS
# and linker flags.  Does NOT modify CMakeLists.txt or build.sh.
#
# Usage:
#   bash basan.sh
#
# Running the resulting binary (recommended):
#   ASAN_OPTIONS=detect_leaks=0:halt_on_error=0 OMP_NUM_THREADS=1 ./your_binary
#
# Notes:
#   - OpenMP is kept enabled (CMakeLists requires it), but you should set
#     OMP_NUM_THREADS=1 at runtime to avoid false positives.
#   - Build type is Debug with -O1 (minimum optimisation for readable ASan
#     traces while keeping some inlining for header-only libs like Eigen).

set -euo pipefail

# ── Paths & logging ──
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
LOG_DIR="${LOG_DIR:-logs}"
mkdir -p "${LOG_DIR}"

ts="$(date +%Y%m%d_%H%M%S)"
log="${LOG_DIR}/basan_${ts}.log"

BUILD_TYPE="Debug"
BUILD_DIR="build-asan"

# ── ccache ──
export CCACHE_CPP2=yes
export CCACHE_NOHASHDIR=1

# ── CPU count (same logic as build.sh) ──
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
echo "total_threads for build: $total_threads"

# ── Compiler selection (same logic as build.sh) ──
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
    : "${CC:=clang}"
    : "${CXX:=clang++}"
  fi
else
  : "${CC:=gcc}"
  : "${CXX:=g++}"
fi

# ── Homebrew libomp hint (macOS only, non-Nix) ──
EXTRA_CMAKE_ARGS=""
if [ "$(uname -s)" = "Darwin" ] && [ -z "${IN_NIX_SHELL:-}" ]; then
  LIBOMP_PREFIX="$(brew --prefix libomp 2>/dev/null)"
  if [ -d "$LIBOMP_PREFIX" ]; then
    EXTRA_CMAKE_ARGS="-DOpenMP_ROOT=$LIBOMP_PREFIX"
  fi
fi

# ── ASan flags ──
ASAN_FLAGS="-fsanitize=address -fno-omit-frame-pointer -fno-optimize-sibling-calls"
# -O2 reproduces Release-like optimisation while keeping ASan instrumentation.
# (Use -O1 if stack traces are unreadable due to inlining.)
ASAN_OPT="-O2"

CMAKE_LOG="cmake-asan.log"
BUILD_LOG="build-asan.log"
rm -f "$CMAKE_LOG" "$BUILD_LOG"

echo "=== ASan debug build ==="
echo "Build dir : $BUILD_DIR"
echo "CC        : $CC"
echo "CXX       : $CXX"
echo "ASan flags: $ASAN_FLAGS $ASAN_OPT"

# ── CMake configure ──
cmake \
  -G Ninja \
  -S "$SCRIPT_DIR" \
  -B "$BUILD_DIR" \
  -DCMAKE_BUILD_TYPE="$BUILD_TYPE" \
  -DCMAKE_C_COMPILER="$CC" \
  -DCMAKE_CXX_COMPILER="$CXX" \
  -DCMAKE_CXX_STANDARD=20 \
  -DCMAKE_CXX_STANDARD_REQUIRED=ON \
  -DCMAKE_CXX_EXTENSIONS=OFF \
  -DCMAKE_C_COMPILER_LAUNCHER=ccache \
  -DCMAKE_CXX_COMPILER_LAUNCHER=ccache \
  -DCMAKE_C_FLAGS="${ASAN_FLAGS} ${ASAN_OPT}" \
  -DCMAKE_CXX_FLAGS="${ASAN_FLAGS} ${ASAN_OPT}" \
  -DCMAKE_EXE_LINKER_FLAGS="${ASAN_FLAGS}" \
  -DCMAKE_SHARED_LINKER_FLAGS="${ASAN_FLAGS}" \
  -DNODEBUG=OFF \
  -DDBGPRINT=OFF \
  $EXTRA_CMAKE_ARGS \
  > "$CMAKE_LOG" 2>&1 || {
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo " CMake configure failed. Check $CMAKE_LOG"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    cat "$CMAKE_LOG"
    exit 1
  }

# ── Build ──
cmake \
  --build "$BUILD_DIR" \
  -j "$total_threads" \
  --verbose \
  > "$BUILD_LOG" 2>&1 || {
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo " Build failed. Check $BUILD_LOG"
    echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    tail -80 "$BUILD_LOG"
    exit 1
  }

echo "####################################################"
echo " OK! ASan build completed successfully."
echo " Binary: $BUILD_DIR/exec/muonith (or similar)"
echo ""
echo " Run with:"
echo "   ASAN_OPTIONS=detect_leaks=0:halt_on_error=0 \\"
echo "   OMP_NUM_THREADS=1 \\"
echo "   ./$BUILD_DIR/exec/<binary>"
echo "####################################################"

# ── Copy logs to timestamped files ──
cp "$CMAKE_LOG" "${LOG_DIR}/cmake-asan_${ts}.log"
cp "$BUILD_LOG" "${LOG_DIR}/build-asan_${ts}.log"
echo "Log (cmake) : ${LOG_DIR}/cmake-asan_${ts}.log"
echo "Log (build) : ${LOG_DIR}/build-asan_${ts}.log"
