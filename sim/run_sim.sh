#!/bin/bash
# run_sim.sh -- simulate the single cells (SLURM array, one task per batch) and
# the 30x matched bulk. Run from sim/ on the login node (it only submits jobs).
#
# Usage: ./run_sim.sh [--cells-only | --bulk-only] [--dry-run] [-- <run_pipeline.sh args>]
#
# Prints the simulation job IDs (colon-separated) on stdout, for --after:
#   SIM=$(./run_sim.sh); ./run_pipeline.sh --phaser ugp --caller chisel --after "$SIM"
#
#   ./run_sim.sh                                      simulate only
#   ./run_sim.sh -- --phaser ugp --caller chisel      simulate, then the pipeline,
#                                                     chained with afterok dependencies
#
# Outputs: ${PROJECT_DIR}/sim/{bams,cn_mat,dwgsim,bulk}

LOG_TAG=run_sim
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
set -euo pipefail

CELLS=true BULK=true DRY=false PIPE_ARGS=()
while (( $# )); do
    case $1 in
        --cells-only) BULK=false; shift ;;
        --bulk-only)  CELLS=false; shift ;;
        --dry-run)    DRY=true; shift ;;
        --)           shift; PIPE_ARGS=("$@"); break ;;
        -h|--help)    sed -n '2,15p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)            die "unknown option '$1' (pipeline options go after --)" ;;
    esac
done

# Validate the pipeline route before spending hours simulating
if (( ${#PIPE_ARGS[@]} )); then
    if ! CHECK=$(./run_pipeline.sh "${PIPE_ARGS[@]}" --after VALIDATE --dry-run 2>&1); then
        echo "${CHECK}" >&2
        die "run_pipeline.sh rejected '${PIPE_ARGS[*]}'; nothing submitted"
    fi
fi

mkdir -p logs

sub() {   # script [extra sbatch args...]
    local script=$1; shift
    if [[ ${DRY} == true ]]; then
        echo "  sbatch --parsable $* ${script}" >&2
        echo "<$(basename "${script}" .sh)>"
    else
        sbatch --parsable "$@" "${script}"
    fi
}

JOBS=()
if [[ ${CELLS} == true ]]; then
    J=$(sub scripts/sim_cells.sh --array=1-"${NUM_BATCHES}")
    note "cells: ${NUM_BATCHES}-task array ${J}"
    JOBS+=("${J}")
fi
if [[ ${BULK} == true ]]; then
    J=$(sub scripts/sim_bulk.sh)
    note "bulk: ${J}"
    JOBS+=("${J}")
fi

# afterok on an array job ID waits for all of its tasks
DEPS=$(IFS=:; echo "${JOBS[*]}")

if (( ${#PIPE_ARGS[@]} )); then
    EXTRA=(); [[ ${DRY} == true ]] && EXTRA=(--dry-run)
    ./run_pipeline.sh "${PIPE_ARGS[@]}" --after "${DEPS}" "${EXTRA[@]+"${EXTRA[@]}"}"
else
    note "Next: ./run_pipeline.sh --phaser <${PHASERS// /|}> --caller <${CALLERS// /|}> --after ${DEPS}"
    echo "${DEPS}"
fi
