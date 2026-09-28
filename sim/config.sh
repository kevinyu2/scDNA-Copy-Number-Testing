#!/bin/bash

# config.sh

PROJECT_DIR="/n/fs/ragr-data/users/ky8418/scDNA_sim/projects/proj_1"
CONDA=/n/fs/ragr-research/users/ky8418/miniconda3/etc/profile.d/conda.sh
ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/scdna_pipeline"

# Simulator
NUM_BATCHES=20
SIM_FA_FOLDER="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/"
TREE_FILE="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/trees/tree1.tree"

# Simulation parameters
ADO_FREQ=0.00001
ADO_MEAN=3
ADO_VAR=1
COVERAGE_MEAN=0.1
COVERAGE_VAR=0.5


# References
REF=" /n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/hg38/Homo_sapiens_assembly38.fasta"

# External tools
DWGSIM="dwgsim"

# Downstream tools — we'll fill these in later
SCAN2_DIR=""
HISCANNER_DIR=""
