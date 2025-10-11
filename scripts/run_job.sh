#!/bin/bash
#
# SLURM Submission Script for klft
#
# Usage: sbatch run_job.sh <input_yaml>

#SBATCH --job-name=klft_job
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

# --- Setup ---
INPUT_YAML=$1
if [ -z "$INPUT_YAML" ]; then
    echo "Error: Please provide an input yaml file."
    exit 1
fi

# Create a unique output directory
OUTPUT_DIR="output/$(basename ${INPUT_YAML%.*})_$(date +%Y%m%d_%H%M%S)"
mkdir -p $OUTPUT_DIR
echo "Output directory: $OUTPUT_DIR"

# Copy input yaml to output directory
cp $INPUT_YAML $OUTPUT_DIR/input.yaml

# --- Data Generation Job ---
echo "Starting data generation job..."
# This is where you would run the klft executable
# For example:
# ./extern/klft/klft $OUTPUT_DIR/input.yaml --output-dir $OUTPUT_DIR
# As a placeholder, we'll just create some dummy files
touch $OUTPUT_DIR/output_data.dat
touch $OUTPUT_DIR/simulation.log
echo "Data generation finished."

# --- Analysis Job ---
echo "Submitting analysis job..."
# This job will run after the data generation job is successfully completed.
sbatch --dependency=afterok:$SLURM_JOB_ID <<EOF
#!/bin/bash
#SBATCH --job-name=analysis_job
#SBATCH --output=${OUTPUT_DIR}/analysis-%j.out
#SBATCH --error=${OUTPUT_DIR}/analysis-%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

echo "Starting analysis job for directory: $OUTPUT_DIR"
# Activate your python environment if needed
# source /path/to/your/venv/bin/activate
python analysis/analysis.py $OUTPUT_DIR
echo "Analysis job finished."
EOF