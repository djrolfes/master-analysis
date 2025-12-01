#!/bin/bash
#
# SLURM Job to Setup Python Virtual Environment
#
#SBATCH --job-name=setup_venv
#SBATCH --output=output/slurm-setup-venv-%j.out
#SBATCH --error=output/slurm-setup-venv-%j.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G

set -euo pipefail

module purge
module load Python

echo "Setting up Python virtual environment..."

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

echo "Python environment setup completed successfully."
echo "Virtual environment is ready at: $(pwd)/.venv"
