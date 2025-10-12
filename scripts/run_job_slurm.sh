#!/bin/bash
#
# SLURM Submission Script for klft
#
# Usage: ./run_job.sh <input_yaml>

INPUT_YAML=$1
if [ -z "$INPUT_YAML" ]; then
    echo "Error: Please provide an input yaml file."
    exit 1
fi

# Convert input path to full path
INPUT_YAML=$(realpath "$INPUT_YAML")

# Create a unique output directory
OUTPUT_DIR="$(realpath output)/$(basename ${INPUT_YAML%.*})_$(date +%Y%m%d_%H%M%S)"
echo "Output directory: $OUTPUT_DIR"

# Submit data generation job
sbatch scripts/data_generation.slurm "$INPUT_YAML" "$OUTPUT_DIR"

# Submit analysis job with dependency on data generation
sbatch --dependency=singleton scripts/analysis.slurm "$OUTPUT_DIR"