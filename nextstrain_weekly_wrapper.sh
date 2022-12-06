#!/usr/bin/env bash

# Usage: sh nextstrain_weekly_wrapper.sh <gisaid-metadata.tar> <gisaid-fasta.tar> <today's date "2022.11.23">

# Pull the latest scripts
cd /home/jonr/Prosjekter/Nextstrain_mamba/FHI-scripts/
git pull origin master

# Move the latest build file into ncov. Remember to update the pango list on GitHub first
mv builds.yaml /home/jonr/Prosjekter/Nextstrain_mamba/ncov/my_profiles/omicron/builds.yaml

# Untar Gisaid files
cd /media/jonr/SATA6TB1/Gisaid/
tar -xf $1
rm $1
tar -xf $2
rm $2

# Clean up old files
rm /home/jonr/Prosjekter/Nextstrain_mamba/ncov/data/SC2_weekly/*{tsv,fasta}

# Index the fasta file
Rscript /home/jonr/Prosjekter/Nextstrain_mamba/FHI-scripts/index_fasta.R

# Parse the Gisaid files
Rscript /home/jonr/Prosjekter/Nextstrain_mamba/FHI-scripts/parse_Gisaid_fasta_and_metadata.R

# Get data from BN
Rscript /home/jonr/Prosjekter/Nextstrain_mamba/FHI-scripts/get_data_from_BN.R

# Update the Nextstrain conda environment
mamba update -n base conda mamba

# Update Nextstrain
source activate nextstrain
mamba update --all
nextstrain update
# conda deactivate

# Pull the latest ncov Nextstrain updates from GitHub
cd /home/jonr/Prosjekter/Nextstrain_mamba/ncov
git pull origin master

# Uodate the nextstrain cli
python3 -m pip install --upgrade nextstrain-cli

# Start the build
source activate nextstrain
cd /home/jonr/Prosjekter/Nextstrain_mamba/ncov
nextstrain build . --configfile my_profiles/omicron/builds.yaml --cores all
# conda deactivate

# Rename the final auspice files
mv auspice/ncov_omicron_ba_five_tip-frequencies.json auspice/$3-ncov_omicron_ba_five_tip-frequencies.json
mv auspice/ncov_omicron_ba_five.json auspice/$3-ncov_omicron_ba_five.json
mv auspice/ncov_omicron_ba_five_root-sequence.json auspice/$3-ncov_omicron_ba_five_root-sequence.json
