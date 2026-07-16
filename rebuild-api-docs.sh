#!/usr/bin/env bash

rm -rf apidocs/

# Generate doxygen documentation
doxygen Doxyfile &

wait  # wait for doxygen to finish

