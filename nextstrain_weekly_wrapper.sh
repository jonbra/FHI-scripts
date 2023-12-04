#!/usr/bin/env bash

# Usage: sh nextstrain_weekly_wrapper.sh <gisaid-metadata.tar> <gisaid-fasta.tar> <today's date "2022.11.23">

# Pull the latest scripts
cd /home/jonr/Prosjekter/Nextstrain_mamba/FHI-scripts/
git pull origin master

# Move the latest build file into ncov. Remember to update the pango list on GitHub first
cp builds.yaml /home/jonr/Prosjekter/Nextstrain_mamba/ncov/my_profiles/omicron/

# Untar Gisaid files
cd /media/jonr/SATA6TB/Gisaid/
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

# Update Nextstrain
nextstrain update

# Pull the latest ncov Nextstrain updates from GitHub
cd /home/jonr/Prosjekter/Nextstrain_mamba/ncov
git pull origin master

# Start the build
# First clean up old files
nextstrain build . clean --configfile my_profiles/omicron/builds.yaml
# And then start the build
nextstrain build . --configfile my_profiles/omicron/builds.yaml --cores all --forceall

# Rename the final auspice files
mv auspice/ncov_omicron_ba_2_86.json auspice/$3-ncov_omicron_ba_2_86.json
mv auspice/ncov_omicron_ba_2_86_tip-frequencies.json auspice/$3-ncov_omicron_ba_2_86_tip-frequencies.json
mv auspice/ncov_omicron_ba_2_86_root-sequence.json auspice/$3-ncov_omicron_ba_2_86_root-sequence.json
#mv auspice/ncov_omicron_ba_five_tip-frequencies.json auspice/$3-ncov_omicron_ba_five_tip-frequencies.json
#mv auspice/ncov_omicron_ba_five.json auspice/$3-ncov_omicron_ba_five.json
#mv auspice/ncov_omicron_ba_five_root-sequence.json auspice/$3-ncov_omicron_ba_five_root-sequence.json
#mv auspice/ncov_omicron_bq_tip-frequencies.json auspice/$3-ncov_omicron_bq_tip-frequencies.json
#mv auspice/ncov_omicron_bq.json auspice/$3-ncov_omicron_bq.json
#mv auspice/ncov_omicron_bq_root-sequence.json auspice/$3-ncov_omicron_bq_root-sequence.json

# Uploading the builds
#nextstrain login
#nextstrain remote upload nextstrain.org/groups/niph auspice/$3-ncov_omicron_xbb*
