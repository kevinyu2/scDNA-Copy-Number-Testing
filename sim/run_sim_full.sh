#!/bin/bash

#SBATCH --job-name=scDNA_pipeline
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --output=logs/scDNA_full_%j.out
#SBATCH --error=logs/scDNA_full_%j.err

# ============================================================
# Usage
#
#   ./run_pipeline.sh START STOP
#
# Examples:
#
#   ./run_pipeline.sh 1 3    # run all stages
#   ./run_pipeline.sh 1 1    # simulation + BAM
#   ./run_pipeline.sh 2 2    # SCAN2
#   ./run_pipeline.sh 3 3    # HiScanner
#   ./run_pipeline.sh 2 3    # SCAN2 + HiScanner
# ============================================================

if [[ $# -ne 2 ]]; then
    echo "Usage: $0 START_STAGE STOP_STAGE"
    echo
    echo "Stages:"
    echo "  1 = simulation + FASTQ → BAM"
    echo "  2 = SCAN2"
    echo "  3 = HiScanner"
    exit 1
fi

START_STAGE="$1"
STOP_STAGE="$2"

if (( START_STAGE < 1 || START_STAGE > 3 )); then
    echo "ERROR: START_STAGE must be 1, 2, or 3"
    exit 1
fi

if (( STOP_STAGE < 1 || STOP_STAGE > 3 )); then
    echo "ERROR: STOP_STAGE must be 1, 2, or 3"
    exit 1
fi

if (( START_STAGE > STOP_STAGE )); then
    echo "ERROR: START_STAGE cannot be greater than STOP_STAGE"
    exit 1
fi

# ============================================================
# Load configuration
# ============================================================

SCRIPT_DIR="./scripts/"
source "./config.sh"

# ============================================================
# Conda
# ============================================================

source ${CONDA}
conda activate ${ENV}

set -euo pipefail

mkdir -p logs

# Dependency flags for the next stage (empty unless an earlier stage ran in this invocation)
DEP_ARGS=()

# ---------------- Stage 1 ----------------
if (( START_STAGE <= 1 && STOP_STAGE >= 1 )); then
    SIM_JOB=$(sbatch --parsable --array=1-"${NUM_BATCHES}" "${SCRIPT_DIR}/run_dwgsim.sh")
    BULK_JOB=$(sbatch --parsable "${SCRIPT_DIR}/make_bulk.sh")
    echo "Submitted simulation array: ${SIM_JOB}, bulk: ${BULK_JOB}"

    # afterok on an array job ID waits for ALL its tasks
    DEP_ARGS=(--dependency=afterok:${SIM_JOB}:${BULK_JOB} --kill-on-invalid-dep=yes)
fi

# ---------------- Stage 2 ----------------
if (( START_STAGE <= 2 && STOP_STAGE >= 2 )); then
    S2_JOB=$(sbatch --parsable ${DEP_ARGS[@]+"${DEP_ARGS[@]}"} \
        --job-name=scan2 --mem=8G --time=72:00:00 \
        --output=logs/scan2_%j.out --error=logs/scan2_%j.err \
        "${SCRIPT_DIR}/scan2_run.sh")
    echo "Submitted SCAN2: ${S2_JOB}"

    DEP_ARGS=(--dependency=afterok:${S2_JOB} --kill-on-invalid-dep=yes)
fi

# ---------------- Stage 3 (same pattern) ----------------
if (( START_STAGE <= 3 && STOP_STAGE >= 3 )); then
    S3_JOB=$(sbatch --parsable ${DEP_ARGS[@]+"${DEP_ARGS[@]}"} "${SCRIPT_DIR}/hiscanner_run.sh")
    echo "Submitted HiScanner: ${S3_JOB}"
fi