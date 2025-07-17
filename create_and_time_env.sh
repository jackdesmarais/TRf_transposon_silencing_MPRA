#!/bin/bash

# Script to create conda environment from YAML file, time the process, and delete it
# Usage: ./create_and_time_env.sh [environment_name]

# Set default environment name if not provided

echo ""
echo "=========================================="
echo "=== PART 1: Conda Environment Creation ==="
echo "=========================================="
ENV_NAME=${1:-"conda_test_environment"}
echo "Environment name: $ENV_NAME"
echo "YAML file: graph_env.yml"
echo ""

# Check if YAML file exists
if [ ! -f "graph_env.yml" ]; then
    echo "Error: graph_env.yml file not found!"
    exit 1
fi

# Check if environment already exists
if conda env list | grep -q "^$ENV_NAME "; then
    echo "Warning: Environment '$ENV_NAME' already exists. Removing it first..."
    conda env remove -n "$ENV_NAME" -y
fi

echo "Starting environment creation..."
echo ""

# Record start time
START_TIME=$(date +%s)

# Create environment from YAML file with custom name
conda env create -f graph_env.yml --name "$ENV_NAME"

# Record end time
END_TIME=$(date +%s)

# Calculate duration
DURATION=$((END_TIME - START_TIME))

echo ""
echo "=== Environment Creation Complete ==="
echo "Environment name: $ENV_NAME"
echo "Time taken: ${DURATION} seconds"
echo ""

# List the created environment
echo "Created environment details:"
conda env list | grep "$ENV_NAME"

# Show installed packages in the environment
echo "=== Installed Packages ==="
conda list -n "$ENV_NAME"

echo ""
echo "=== Deleting Environment ==="

# Delete the environment
conda env remove -n "$ENV_NAME" -y

echo "Environment '$ENV_NAME' has been deleted."
echo ""
echo "=== Summary ==="
echo "Total time to install with conda: ${DURATION} seconds"
echo "Script completed successfully!" 

echo ""
echo "=========================================="
echo "=== PART 2: Manual Environment Creation ==="
echo "=========================================="

# Set environment name for part 2
ENV_NAME_2=${2:-"manual_test_env"}

echo "Environment name: $ENV_NAME_2"
echo ""

# Check if environment already exists
if conda env list | grep -q "^$ENV_NAME_2 "; then
    echo "Warning: Environment '$ENV_NAME_2' already exists. Removing it first..."
    conda env remove -n "$ENV_NAME_2" -y
fi

echo "Creating conda environment with Python 3.9 and pip..."
conda create -n "$ENV_NAME_2" python=3.9 pip -y

echo ""
echo "Activating environment..."
source $(conda info --base)/etc/profile.d/conda.sh
conda activate "$ENV_NAME_2"

echo ""
echo "=== Installing library_analyzer ==="
echo "Starting library_analyzer installation..."

# Record start time for library_analyzer installation
START_TIME_2=$(date +%s)

# Install library_analyzer with pip -e
pip install --no-cache-dir -e ".[notebooks]"

# Record end time
END_TIME_2=$(date +%s)

# Calculate duration
DURATION_2=$((END_TIME_2 - START_TIME_2))

echo ""
echo "=== Library Analyzer Installation Complete ==="
echo "Time taken for library_analyzer installation: ${DURATION_2} seconds"
echo ""

# List installed packages
echo "Installed packages:"
pip list

echo ""
echo "=== Deactivating and Deleting Environment ==="

# Deactivate environment
conda deactivate

# Delete the environment
conda env remove -n "$ENV_NAME_2" -y

echo "Environment '$ENV_NAME_2' has been deleted."
echo ""
echo "=== Final Summary ==="
echo "Time for YAML-based environment creation: ${DURATION} seconds"
echo "Time for library_analyzer installation only: ${DURATION_2} seconds"
echo "Script completed successfully!" 


