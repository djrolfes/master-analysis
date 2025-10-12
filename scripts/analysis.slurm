#!/bin/bash
#
# SLURM Job Script for Analysis
#
#SBATCH --job-name=wf_analysis
#SBATCH --output=slurm-analysis-%j.out
#SBATCH --error=slurm-analysis-%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

set -euo pipefail
module purge
module load GCC/11.2.0
module load R
module load Python

cd "$(dirname "$0")/.."

# --- Setup ---
OUTPUT_DIR=$(realpath "$1")

if [ -z "$OUTPUT_DIR" ]; then
    echo "Error: Missing argument. Usage: sbatch analysis.slurm <output_dir>"
    exit 1
fi


# --- Check and Install Hadron ---
echo "Running check_and_install_hadron.sh..."
if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "Running with srun..."
    srun ./scripts/check_and_install_hadron.sh
else
    echo "Running without srun..."
    ./scripts/check_and_install_hadron.sh
fi

# --- Python Environment Setup ---
if [ ! -d ".venv" ]; then
    python3 -m venv .venv
fi

source .venv/bin/activate
pip install --upgrade pip
pip install pandas PyYAML
pip install matplotlib seaborn

# --- Run Analysis ---
cd analysis
echo "Starting analysis job for directory: $OUTPUT_DIR"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    echo "Running with srun..."
    srun python3 analysis.py "$OUTPUT_DIR"
else
    echo "Running without srun..."
    python3 analysis.py "$OUTPUT_DIR"
fi
echo "Analysis job finished."

mv *.pdf "$OUTPUT_DIR"
echo "Copied analysis PDFs to $OUTPUT_DIR"