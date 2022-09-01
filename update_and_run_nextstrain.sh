#!/bin/bash

# Update the conda environment
mamba update -n base conda mamba
source activate nextstrain
mamba update --all
nextstrain update
conda deactivate

# Pull the latest updates from GitHub
cd ncov
git pull

# Uodate the nextstrain cli
python3 -m pip install --upgrade nextstrain-cli

# Start the build
source activate nextstrain
cd ncov
nextstrain build . --configfile my_profiles/omicron/builds.yaml --cores all
conda deactivate
