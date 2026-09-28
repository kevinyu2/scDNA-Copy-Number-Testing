#!/bin/bash

# config.sh

PROJECT_DIR="/n/fs/ragr-data/users/ky8418/scDNA_sim/projects/proj_COLO_mat"
CONDA=/n/fs/ragr-research/users/ky8418/miniconda3/etc/profile.d/conda.sh
ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/scdna_pipeline"

# Simulator
NUM_BATCHES=20
SIM_MAT_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_mat_chr1.fasta"
SIM_PAT_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_pat_chr1.fasta"
# This one gets created automatically if not present
SIM_FULL_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/test_genomes/HG002_full_chr1.fasta"
TREE_FILE="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/trees/mat_COLO2.tree"

# Simulation parameters
ADO_FREQ=0.000005
ADO_MEAN=3
ADO_VAR=1
COVERAGE_MEAN=0.1
COVERAGE_VAR=0.3


# References
REF="/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/hg38/Homo_sapiens_assembly38.fasta"


# Downstream tools — we'll fill these in later
SCAN2_DIR=""
HISCANNER_DIR=""
