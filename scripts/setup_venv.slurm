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

echo "=== Python Virtual Environment Setup ==="
echo "Working directory: $(pwd)"
echo "Python version: $(python3 --version)"
echo ""

if [ -d ".venv" ]; then
    echo "Removing existing .venv directory..."
    rm -rf .venv
    echo "Old .venv removed successfully"
fi

echo ""
echo "Creating fresh virtual environment..."
python3 -m venv .venv

echo "Activating virtual environment..."
source .venv/bin/activate

echo "Python in venv: $(which python3)"
echo "Pip in venv: $(which pip)"
echo ""

echo "Upgrading pip..."
pip install --upgrade pip

echo ""
echo "Installing compatible numpy version..."
pip install --no-cache-dir "numpy>=1.16.5,<1.23.0"

echo ""
echo "Installing other packages..."
pip install --no-cache-dir pandas PyYAML matplotlib seaborn

echo ""
echo "=== Verifying Installation ==="
pip list | grep -E "(numpy|pandas|PyYAML|matplotlib|seaborn)"

echo ""
echo "Testing imports..."
python3 -c "import pandas; import yaml; import matplotlib; import seaborn; import numpy; print('All imports successful!')"

echo ""
echo "=== Python environment setup completed successfully ==="
echo "Virtual environment location: $(pwd)/.venv"
