#!/bin/bash
# One command, overnight: simulate a sensitivity tree as its own project, make UGP bins for each bin size, and for each
# bin size run + eval (detection heatmap only):  plain (no HMM)  |  oracle  |  HMM K=SUGGEST_K  |  HMM find-K  |  HMM find-K + cloner
#
#   bash run_tree_sweep.sh sens_amp_b20_s2          # tree data/trees/sens_amp_b20_s2.tree, project <projects dir>/sens_amp_b20_s2
#   BINS="8 32" SUGGEST_K=40 FIND_RANGE=24-48 bash run_tree_sweep.sh sens_amp_b20_s2
#   bash run_tree_sweep.sh sens_amp_b20_s2 --dry-run   # only writes the config and prints what it would do
# Run it on the LOGIN node (not sbatch). Defaults: BINS="4 8 16 32", SUGGEST_K=40, FIND_RANGE=24-48.
# The config of this tree is written to config_<tree>.sh (a copy of config.sh with PROJECT_DIR and TREE_FILE changed);
# every job gets it through CONFIG_FILE. Needs the patch with EVAL_DETECTION_ONLY (eval_calls.py --detection-only).
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
TREE=${1:?usage: bash run_tree_sweep.sh <tree name, e.g. sens_amp_b20_s2> [--dry-run]}
DRY=false; [[ ${2:-} == --dry-run ]] && DRY=true
BINS=(${BINS:-4 8 16 32})
SUGGEST_K=${SUGGEST_K:-40}
FIND_RANGE=${FIND_RANGE:-24-48}

TREEDIR=$(dirname "$(grep -E '^TREE_FILE=' config.sh | head -1 | cut -d'"' -f2)")
PROJ=$(dirname "$(grep -E '^PROJECT_DIR=' config.sh | head -1 | cut -d'"' -f2)")/${TREE}
[[ -s ${TREEDIR}/${TREE}.tree && -s ${TREEDIR}/${TREE}.design.tsv ]] || { echo "ERROR: no ${TREEDIR}/${TREE}.tree / .design.tsv" >&2; exit 1; }
CFG=$(pwd)/config_${TREE}.sh
sed -E "s|^PROJECT_DIR=.*|PROJECT_DIR=\"${PROJ}\"|; s|^TREE_FILE=.*|TREE_FILE=\"${TREEDIR}/${TREE}.tree\"|" config.sh > "${CFG}"
echo "config ${CFG}: PROJECT_DIR=${PROJ} TREE_FILE=${TREEDIR}/${TREE}.tree  bins: ${BINS[*]}  K=${SUGGEST_K}  find-K ${FIND_RANGE}" >&2
export CONFIG_FILE=${CFG} EVAL_DETECTION_ONLY=${EVAL_DETECTION_ONLY:-true}
[[ ${DRY} == true ]] && { echo "dry run: stopping here" >&2; exit 0; }

SIM=$(./run_sim.sh)
echo "simulation: ${SIM}" >&2
PREV=${SIM}
for n in "${BINS[@]}"; do                       # UGP bins, one at a time (they share one pseudobulk folder)
    PREV=$(./run_pipeline.sh --phaser ugp --to phase --phaserbins "${n}" --after "${PREV}")
    echo "bins ${n}: UGP ${PREV}" >&2
    C=(./run_pipeline.sh --phaser ugp --caller chisel --phaserbins "${n}" --from call --to eval)
    PLAIN=$("${C[@]}" --after "${PREV}")
    echo "bins ${n}: plain+eval ${PLAIN}" >&2
    for variant in "--oracle" "--hmm ${SUGGEST_K}" "--hmm ${FIND_RANGE}" "--hmm ${FIND_RANGE} --cloner"; do
        # shellcheck disable=SC2086
        J=$("${C[@]}" ${variant} --after "${PLAIN}")
        echo "bins ${n}: ${variant} ${J}" >&2
    done
done
echo "submitted; heatmaps: ${PROJ}/eval/*/detection.png (pool trees with scripts/aggregate_detection.py)" >&2
