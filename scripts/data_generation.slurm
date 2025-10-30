#!/bin/bash
#
# SLURM Job Script for Data Generation
#
#SBATCH --job-name=wf_data_generation
#SBATCH --output=output/slurm-data-%j.out
#SBATCH --error=output/slurm-data-%j.err
#SBATCH --time=3-24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=32G
#SBATCH --gres=gpu:pascal:1
#SBATCH --nodelist=lnode[13-17]

set -euo pipefail
if [ -n "${SLURM_JOB_ID:-}" ]; then
    module purge
    module load CUDA
    module load GCC/11.2.0
    module load CMake
    module load OpenMPI
    export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
else
    export OMP_NUM_THREADS=8
fi


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

# Use srun if within a SLURM job, otherwise run directly
if [ -n "${SLURM_JOB_ID:-}" ]; then
    # Append the SLURM job ID to the output directory
    OUTPUT_DIR="${OUTPUT_DIR}_$SLURM_JOB_ID"
    mkdir -p "$OUTPUT_DIR"
    echo "Output directory: $OUTPUT_DIR"
    # Copy input yaml to output directory
    cp "$INPUT_YAML" "$OUTPUT_DIR/$(basename "$INPUT_YAML")"
    echo "Running with srun..."
    srun ./extern/klft/build/binaries/hmc -f "$OUTPUT_DIR/$(basename "$INPUT_YAML")" -o "$OUTPUT_DIR"
else
    # Copy input yaml to output directory
    mkdir -p "$OUTPUT_DIR"
    echo "Output directory: $OUTPUT_DIR"
    cp "$INPUT_YAML" "$OUTPUT_DIR/$(basename "$INPUT_YAML")"
    echo "Running without srun..."
    ./extern/klft/build/binaries/hmc -f "$OUTPUT_DIR/$(basename "$INPUT_YAML")" -o "$OUTPUT_DIR"
fi

echo "Data generation finished."