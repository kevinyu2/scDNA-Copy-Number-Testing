#!/bin/bash

# config.sh

PROJECT_DIR="/n/fs/ragr-data/users/ky8418/scDNA_sim/projects/proj_COLO_mat_final"
CONDA=/n/fs/ragr-research/users/ky8418/miniconda3/etc/profile.d/conda.sh
ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/scdna_pipeline"

# Simulator
NUM_BATCHES=20
SIM_MAT_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_mat_chr1.fa"
SIM_PAT_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_pat_chr1.fa"
# This one gets created automatically if not present
SIM_FULL_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_full_chr1.fa"
TREE_FILE="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/trees/mat_COLO2.tree"

# Simulation parameters
ADO_FREQ=0.000005
ADO_MEAN=3
ADO_VAR=1
COVERAGE_MEAN=0.1
COVERAGE_VAR=0.3


# References
REF="/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/hg38/Homo_sapiens_assembly38.fasta"


# Resources for scan2 and hiscanner
RESOURCES_DIR=/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/
WINDOWS_BED="/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/chr1_2000windows.bed"
GENOME="hg38"
CHR="chr1"


# ============================================================
# SLURM 
# ============================================================
SLURM_PARTITION="cs"
SLURM_ACCOUNT="allcs"          # find yours with: sacctmgr show assoc user=$USER format=account


# ============================================================
# HiScanner
# ============================================================
HS_ENV="${ENV}"                                   # change if HiScanner lives in its own env
HS_DIR="${PROJECT_DIR}/hiscanner"                 # project dir: config.yaml, metadata.txt, output/

# Reference files
FASTA_SPLIT_DIR="${RESOURCES_DIR%/}/hg38_split"   # must contain <chrom>.fasta (e.g. chr1.fasta)
# Mappability prefix: HiScanner opens ${MAPPABILITY_STEM}chr1.txt
# Check the real file names after downloading and set the prefix to match.
MAPPABILITY_STEM="${RESOURCES_DIR%/}/hg38_mappability/150_mer."

# Chromosomes to analyze (space-separated, names as in the BAM/VCF)
HS_CHROMS="${CHR}"

# Steps to run, in order (drop ones you've finished; snp/phase/ado are skipped if HS_RDR_ONLY=true)
HS_STEPS="snp phase ado normalize segment cnv"

# Analysis parameters (defaults from HiScanner's default_config.yaml)
HS_RDR_ONLY=false
HS_MULTISAMPLE=false                # use_multisample_segmentation
HS_BINSIZE=500000
HS_LAMBDA=16
HS_LAMBDA_RANGE="[2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]"
HS_MAX_WGD=2
HS_BATCH_SIZE=5
HS_DEPTH_FILTER=0
HS_ADO_THRESHOLD=0.2
HS_ADO_PLOTS=true
HS_AGGREGATE_K=true                 # aggregate_every_k_snp
HS_K=5
HS_THREADS=16

# Run control
HS_USE_CLUSTER=true              
HS_RERUN=false                      # true = redo steps even if outputs exist



# ============================================================
# Universal Genotyping Pipeline
# ============================================================
UGP_ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/genotyping-env"
UGP_REPO="/path/to/Universal-Genotyping-Pipeline"
UGP_CONDA_PREFIX="/n/fs/ragr-data/users/ky8418/scDNA_sim/ugp_conda"   # same as in step 1
UGP_DIR="${PROJECT_DIR}/ugp"
UGP_SAMPLE_ID="HG002sim"

UGP_MODE="bulk"            # bulk | percell
UGP_MAX_CELLS=20           # percell only; 0 = all cells

UGP_CHROMS="1"             # UGP wants bare numbers; it handles the chr prefix itself
UGP_SNP_PANEL="${RESOURCES_DIR%/}/1kGP.chr1.snps.vcf.gz"
UGP_PHASING_PANEL="${RESOURCES_DIR%/}/ugp_phasing_panel"
UGP_GTF="${RESOURCES_DIR%/}/gencode.v38.annotation.gtf.gz"

# Binning (percell). UGP's docs suggest 50–200 for low coverage.
UGP_MIN_SNP_READS="[50, 100, 200]"
UGP_MIN_TOTAL_READS=1000   # read starts per bin per dataset; the default 5000 makes ~7 Mb bins at 0.1x

UGP_CORES=16
UGP_THREADS_STEP=4