#!/bin/bash

# Script to compile the klft module

set -e  # Exit on error

# Define paths
KLFT_DIR="$(dirname $(realpath $0))/../extern/klft"
BUILD_DIR="$KLFT_DIR/build"
CMAKE_CONFIG_FILE="$(dirname $(realpath $0))/setup_home.cmake"

# Create the build directory if it doesn't exist
if [ ! -d "$BUILD_DIR" ]; then
    # Check if CMake configuration file exists
    if [ -f "$CMAKE_CONFIG_FILE" ]; then
        echo "Creating build directory at $BUILD_DIR"
        mkdir -p "$BUILD_DIR"

        # Change to the build directory
        cd "$BUILD_DIR"
        echo "CMake configuration file found. Configuring with CMake."
        cmake -C "$CMAKE_CONFIG_FILE" -S .. -B .
    else
        echo "CMake configuration file not found. Please configure manually before building."
        exit 1
    fi
fi

# Change to the build directory
cd "$BUILD_DIR"

# Run make with 6 parallel jobs
echo "Running make -j"
make -j

echo "Build completed successfully."