#!/usr/bin/env bash
set -euo pipefail

BUILD_TYPE="${1:-Debug}"
LOG_DIR="${LOG_DIR:-logs}"
mkdir -p "${LOG_DIR}"

ts="$(date +%Y%m%d_%H%M%S)"
log="${LOG_DIR}/bdebug_${BUILD_TYPE}_${ts}.log"

# Important: save stdout/stderr combined into one (no fabrication, ensure evidence)
{ time ENABLE_NODEBUG=OFF ENABLE_DBGPRINT=OFF ./build.sh "${BUILD_TYPE}"; } \
  > "${log}" 2>&1

echo "Log: ${log}"
