#!/bin/bash
#
# SLURM Job Script for Building klft
#
#SBATCH --job-name=build_klft
#SBATCH --output=output/slurm-build-%j.out
#SBATCH --error=output/slurm-build-%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

set -e  # Exit on error

# Define paths
KLFT_DIR="$/extern/klft"
BUILD_DIR="$KLFT_DIR/build"
CMAKE_CONFIG_FILE="scripts/setup_home.cmake"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    CMAKE_CONFIG_FILE="scripts/setup_cluster.cmake"
fi


# Create the build directory if it doesn't exist
if [ ! -d "$BUILD_DIR" ]; then
    # Check if CMake configuration file exists
    if [ -f "$CMAKE_CONFIG_FILE" ]; then
        echo "Creating build directory at $BUILD_DIR"
        mkdir -p "$BUILD_DIR"

        # Change to the build directory
        cd "$BUILD_DIR"
        echo "CMake configuration file found. Configuring with CMake."
        cmake -C "$CMAKE_CONFIG_FILE" -S ..
    else
        echo "CMake configuration file not found. Please configure manually before building."
        exit 1
    fi
fi

# Change to the build directory
cd "$BUILD_DIR"

# Run make with all available CPU cores
if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "Running make with srun using all available CPU cores"
    srun make -j$(nproc)
else
    echo "Running make without srun using all available CPU cores"
    make -j$(nproc)
fi

echo "Build completed successfully."