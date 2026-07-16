#!/usr/bin/env bash
set -e

echo "==== Cleaning ccache (Debug) ===="
ccache --clear

echo "==== Removing build-debug directory ===="
rm -rf build-debug

echo "==== Starting clean Debug build ===="
time ENABLE_NODEBUG=OFF ENABLE_DBGPRINT=ON bash ./build.sh Debug
