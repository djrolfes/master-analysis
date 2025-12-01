#!/bin/bash
#
# Shell Script to Dispatch Multiple Analysis Jobs
# Usage: ./scripts/analysis.sh [-skip <n>] <output_dir1> [<output_dir2> ...]
#

set -euo pipefail

module purge
module load Python

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

# --- Setup Python Environment (once, before submitting jobs) ---
echo "Setting up Python virtual environment..."
cd "$PROJECT_ROOT"

if [ -d ".venv" ]; then
    echo "Removing existing .venv to avoid corruption..."
    rm -rf .venv
fi

echo "Creating fresh virtual environment..."
python3 -m venv .venv
source .venv/bin/activate

echo "Installing Python packages with compatible versions..."
pip install --upgrade pip
pip install "numpy>=1.16.5,<1.23.0"
pip install pandas PyYAML matplotlib seaborn

echo "Python environment ready."
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
echo "Job IDs: ${JOB_IDS[*]}"
echo ""
echo "Monitor with: squeue -u \$USER"
echo "Cancel all with: scancel ${JOB_IDS[*]}"
