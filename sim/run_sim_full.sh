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


# ============================================================
# Stage 1
# ============================================================

if (( START_STAGE <= 1 && STOP_STAGE >= 1 )); then
    echo
    echo "============================================================"
    echo "STAGE 1: Simulation + FASTQ → BAM + Bulk"
    echo "============================================================"

    SIM_JOB=$(sbatch --parsable \
        --array=1-"${NUM_BATCHES}" \
        "${SCRIPT_DIR}/run_dwgsim.sh")

    echo "Submitted simulation array: ${SIM_JOB}"

    BULK_JOB=$(sbatch --parsable \
        "${SCRIPT_DIR}/make_bulk.sh")

    echo "Submitted bulk job: ${BULK_JOB}"
fi


# # ============================================================
# # Stage 2
# # ============================================================

if (( START_STAGE <= 2 && STOP_STAGE >= 2 )); then
    echo
    echo "============================================================"
    echo "STAGE 2: SCAN2"
    echo "============================================================"

    bash "${SCRIPT_DIR}/02_scan2.sh"
fi

# # ============================================================
# # Stage 3
# # ============================================================

# if (( START_STAGE <= 3 && STOP_STAGE >= 3 )); then
#     echo
#     echo "============================================================"
#     echo "STAGE 3: HiScanner"
#     echo "============================================================"

#     bash "${SCRIPT_DIR}/03_hiscanner.sh"
# fi

# echo
# echo "============================================================"
# echo "Pipeline finished"
# echo "Sample: ${SAMPLE_NAME}"
# echo "============================================================"
