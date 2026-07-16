#!/usr/bin/env bash
set -e

echo "==== Cleaning ccache (Release) ===="
ccache --clear

echo "==== Removing build-release directory ===="
rm -rf build-release

echo "==== Starting clean Release build ===="
time ENABLE_NODEBUG=ON ENABLE_DBGPRINT=OFF ./build.sh Release
