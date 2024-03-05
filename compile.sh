#!/bin/bash

rm -rf build

# Create a build directory
mkdir -p build
cd build

# Run CMake to generate the build files
cmake ..

# Build the project using the generated build files
cmake --build . -- -j4

# Run the compiled executable
#./main

