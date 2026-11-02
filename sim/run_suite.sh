#!/bin/bash
# run_suite.sh -- simulate once, then every comparison run of the project in config.sh, as chained SLURM jobs:
#   perfect -> CHISEL  (pbins8): plain, --hmm 2-16, --oracle
#   ugp     -> CHISEL  (pbins8): plain, --hmm 2-16, --oracle
#   perfect -> HiScanner
# Run from sim/ on the login node (it only submits jobs). Start from a clean project: delete the old
# simulation and chisel/prep first (or set CHISEL_PREP_FORCE=true), so nothing mixes two simulations.
#
# Usage: ./run_suite.sh [--dry-run] [--no-sim] [--after IDS]
#   --dry-run   print the sbatch commands only
#   --no-sim    skip ./run_sim.sh (the simulation in the project is complete and wanted)
#   --after IDS with --no-sim: wait for these jobs first
#
# Order (A -> B: B starts after A succeeded):
#   simulation -> perfect plain (phase_perfect, chisel_prep, bins_perfect, chisel, eval: lifts the truth)
#   perfect plain -> perfect --hmm 2-16, perfect --oracle (both only re-call on its combo.tsv; side by side)
#   perfect plain -> ugp plain (its phase stage would otherwise build chisel/prep at the same time)
#   ugp plain     -> ugp --hmm 2-16, ugp --oracle
#   perfect plain -> perfect -> HiScanner (from the call stage: reuses phase_perfect's output)
# The --hmm/--oracle runs and HiScanner write their own folders and reuse the lifted truth, so they
# can run side by side. Job IDs of every route are printed at the end and kept in logs/suite_<date>.txt.

set -euo pipefail
cd "$(dirname "$0")"
DRY=() NOSIM=false AFTER=""
while (( $# )); do
    case $1 in
        --dry-run) DRY=(--dry-run); shift ;;
        --no-sim)  NOSIM=true; shift ;;
        --after)   AFTER=${2:?--after needs job IDs}; shift 2 ;;
        -h|--help) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)         echo "ERROR: unknown option '$1'" >&2; exit 1 ;;
    esac
done

rp() { ./run_pipeline.sh "$@" "${DRY[@]+"${DRY[@]}"}"; }      # prints the route's last job ID
CHISEL=(--caller chisel --phaserbins 8)

if [[ ${NOSIM} == true ]]; then
    SIM=${AFTER}
else
    SIM=$(./run_sim.sh "${DRY[@]+"${DRY[@]}"}")
fi
after() { [[ -n $1 ]] && echo "--after $1"; }

# shellcheck disable=SC2046
P=$(rp --phaser perfect "${CHISEL[@]}" $(after "${SIM}"))
PH=$(rp --phaser perfect "${CHISEL[@]}" --hmm 2-16 --from call --after "${P}")
PO=$(rp --phaser perfect "${CHISEL[@]}" --oracle --from call --after "${P}")
U=$(rp --phaser ugp "${CHISEL[@]}" --after "${P}")
UH=$(rp --phaser ugp "${CHISEL[@]}" --hmm 2-16 --from call --after "${U}")
UO=$(rp --phaser ugp "${CHISEL[@]}" --oracle --from call --after "${U}")
H=$(rp --phaser perfect --caller hiscanner --from call --after "${P}")

mkdir -p logs
OUT="logs/suite_$(date +%Y%m%d_%H%M%S).txt"
{
    echo "simulation:                 ${SIM:-none}"
    echo "perfect chisel pbins8:      ${P}"
    echo "perfect chisel hmm 2-16:    ${PH}"
    echo "perfect chisel oracle:      ${PO}"
    echo "ugp chisel pbins8:          ${U}"
    echo "ugp chisel hmm 2-16:        ${UH}"
    echo "ugp chisel oracle:          ${UO}"
    echo "perfect hiscanner:          ${H}"
} | tee "${OUT}" >&2
[[ ${#DRY[@]} -eq 0 ]] && echo "(job IDs in ${OUT}; last jobs above, squeue -u \$USER for the rest)" >&2 || rm -f "${OUT}"
