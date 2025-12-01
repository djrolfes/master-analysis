#!/bin/bash
#
# Shell Script to Dispatch Multiple Analysis Jobs
# Usage: ./scripts/analysis.sh [-skip <n>] <output_dir1> [<output_dir2> ...]
#

set -euo pipefail

# --- Parse optional -skip flag ---
SKIP_STEPS=250
if [ "$#" -ge 2 ] && [ "$1" = "-skip" ]; then
    SKIP_STEPS="$2"
    shift 2
fi

# --- Expect one or more output directories ---
if [ "$#" -lt 1 ]; then
    echo "Error: Missing argument. Usage: ./scripts/analysis.sh [-skip <n>] <output_dir1> [<output_dir2> ...]"
    exit 1
fi

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

# --- Submit Python Environment Setup Job ---
echo "Submitting Python environment setup job..."
VENV_JOB_OUTPUT=$(sbatch --parsable "${SCRIPT_DIR}/setup_venv.slurm")
VENV_JOB_ID=$(echo "$VENV_JOB_OUTPUT" | grep -oP '^\d+' || echo "$VENV_JOB_OUTPUT")

echo "Virtual environment setup job submitted: $VENV_JOB_ID"
echo "Waiting for environment setup to complete..."
echo ""

# Track submitted job IDs
JOB_IDS=()

echo "Submitting analysis jobs with skip=$SKIP_STEPS..."
echo ""

# Submit one job per directory
for OUTPUT_DIR in "$@"; do
    ABS_OUTPUT_DIR="$(realpath "$OUTPUT_DIR")"
    
    if [ ! -d "$ABS_OUTPUT_DIR" ]; then
        echo "Warning: Directory does not exist: $ABS_OUTPUT_DIR. Skipping."
        continue
    fi
    
    # Extract directory name for job name
    DIR_NAME="$(basename "$ABS_OUTPUT_DIR")"
    
    JOB_OUTPUT=$(sbatch \
        --dependency=afterok:${VENV_JOB_ID} \
        --job-name="analysis_${DIR_NAME}" \
        --output="output/slurm-analysis-${DIR_NAME}-%j.out" \
        --error="output/slurm-analysis-${DIR_NAME}-%j.err" \
        "${SCRIPT_DIR}/analysis_single_output.slurm" -skip "$SKIP_STEPS" "$ABS_OUTPUT_DIR")
    
    # Extract job ID from "Submitted batch job 12345"
    JOB_ID=$(echo "$JOB_OUTPUT" | grep -oP '\d+$' || echo "$JOB_OUTPUT")
    JOB_IDS+=("$JOB_ID")
    
    echo "Submitted job $JOB_ID for: $ABS_OUTPUT_DIR"
done
echo ""
echo "All jobs submitted!"
echo "Setup job ID: $VENV_JOB_ID"
echo "Analysis job IDs: ${JOB_IDS[*]}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all with: scancel $VENV_JOB_ID ${JOB_IDS[*]}"
echo "Cancel all with: scancel ${JOB_IDS[*]}"
