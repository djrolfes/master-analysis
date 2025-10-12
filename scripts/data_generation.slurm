#!/bin/bash
#
# SLURM Job Script for Data Generation
#
#SBATCH --job-name=wf_data_generation
#SBATCH --output=slurm-data-%j.out
#SBATCH --error=slurm-data-%j.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=32G
#SBATCH --gres:pascal:1

set -euo pipefail
module purge
module load CUDA
module load GCC/11.2.0
module load CMake
module load OpenMPI

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
export KLFT_TUNING=0
export KLFT_VERBOSITY=1

# --- Setup ---
INPUT_YAML=$1
OUTPUT_DIR=$2

if [ -z "$INPUT_YAML" ] || [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing arguments. Usage: sbatch data_generation.slurm <input_yaml> <output_dir>"
    exit 1
fi

mkdir -p "$OUTPUT_DIR"
echo "Output directory: $OUTPUT_DIR"

# Copy input yaml to output directory
cp "$INPUT_YAML" "$OUTPUT_DIR/input.yaml"

# --- Build and Data Generation ---
echo "Running build script..."
echo "Running build script..."
if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "Building with srun..."
    srun ./scripts/build_klft.sh
else
    echo "Building without srun..."
    ./scripts/build_klft.sh
fi

export KLFT_TUNING=0
export KLFT_VERBOSITY=1

# Use srun if within a SLURM job, otherwise run directly
if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "Running with srun..."
    srun ./extern/klft/build/tests/wilsonflow_sp_test -f "$OUTPUT_DIR/input.yaml" -o "$OUTPUT_DIR"
else
    echo "Running without srun..."
    ./extern/klft/build/tests/wilsonflow_sp_test -f "$OUTPUT_DIR/input.yaml" -o "$OUTPUT_DIR"
fi

echo "Data generation finished."