#!/bin/bash
#
# Shell Script to Dispatch Multiple Analysis Jobs
# Usage: ./scripts/analysis.sh [-setup] [-skip <n>] <output_dir1> [<output_dir2> ...]
#

set -euo pipefail

# --- Parse optional flags ---
SETUP_VENV=false
SKIP_STEPS=250

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -setup)
            SETUP_VENV=true
            shift
            ;;
        -skip)
            if [ "$#" -lt 2 ]; then
                echo "Error: -skip requires a value"
                exit 1
            fi
            SKIP_STEPS="$2"
            shift 2
            ;;
        *)
            break
            ;;
    esac
done

# --- Expect one or more output directories ---
if [ "$#" -lt 1 ]; then
    echo "Error: Missing argument. Usage: ./scripts/analysis.sh [-setup] [-skip <n>] <output_dir1> [<output_dir2> ...]"
    exit 1
fi

# Get the script directory and project root
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

cd "$PROJECT_ROOT"

# --- Conditionally Submit Python Environment Setup Job ---
VENV_JOB_ID=""
DEPENDENCY_FLAG=""

if [ "$SETUP_VENV" = true ]; then
    echo "Submitting Python environment setup job..."
    VENV_JOB_OUTPUT=$(sbatch --parsable "${SCRIPT_DIR}/setup_venv.slurm")
    VENV_JOB_ID=$(echo "$VENV_JOB_OUTPUT" | grep -oP '^\d+' || echo "$VENV_JOB_OUTPUT")
    
    echo "Virtual environment setup job submitted: $VENV_JOB_ID"
    echo "Analysis jobs will wait for environment setup to complete."
    DEPENDENCY_FLAG="--dependency=afterok:${VENV_JOB_ID}"
    echo ""
else
    echo "Using existing virtual environment (add -setup flag to recreate)"
    echo ""
fi

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
        ${DEPENDENCY_FLAG} \
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
if [ -n "$VENV_JOB_ID" ]; then
    echo "Setup job ID: $VENV_JOB_ID"
    echo "Analysis job IDs: ${JOB_IDS[*]}"
    echo ""
    echo "Monitor with: squeue -u \$USER"
    echo "Cancel all with: scancel $VENV_JOB_ID ${JOB_IDS[*]}"
else
    echo "Analysis job IDs: ${JOB_IDS[*]}"
    echo ""
    echo "Monitor with: squeue -u \$USER"
    echo "Cancel all with: scancel ${JOB_IDS[*]}"
fi
echo "Cancel all with: scancel ${JOB_IDS[*]}"
