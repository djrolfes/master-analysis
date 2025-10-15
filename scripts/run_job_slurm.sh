#!/bin/bash
#
# SLURM Submission Script for klft
#
# Usage: ./run_job_slurm.sh <input_yaml1> <input_yaml2> ...

if [ "$#" -lt 1 ]; then
    echo "Error: Please provide at least one input yaml file."
    exit 1
fi

# Submit the build job
BUILD_JOB_ID=$(sbatch --parsable scripts/build_klft.slurm)
echo "Build job submitted with Job ID: $BUILD_JOB_ID"

# Iterate over all input YAML files
for INPUT_YAML in "$@"; do
    # Convert input path to full path
    INPUT_YAML=$(realpath "$INPUT_YAML")

    # Create a base output directory
    BASE_OUTPUT_DIR="$(realpath output)/$(basename ${INPUT_YAML%.*})_$(date +%Y%m%d_%H%M%S)"

    # Submit data generation job with dependency on build job
    DATA_JOB_ID=$(sbatch --dependency=afterok:$BUILD_JOB_ID --parsable scripts/data_generation.slurm "$INPUT_YAML" "$BASE_OUTPUT_DIR")
    echo "Data generation job submitted with Job ID: $DATA_JOB_ID"

    # Update the output directory to include the DATA_JOB_ID
    OUTPUT_DIR="${BASE_OUTPUT_DIR}_$DATA_JOB_ID"

    # Submit analysis job with dependency on successful data generation
    sbatch --dependency=afterok:$DATA_JOB_ID scripts/analysis.slurm "$OUTPUT_DIR"
    echo "Analysis job submitted with dependency on Job ID: $DATA_JOB_ID"

done