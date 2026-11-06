#!/bin/bash
# Detection sweep: for every UGP bin size, five CHISEL runs, each followed by its eval (which now also writes
# eval/<run>/detection.png, the subclone size x CNA size heatmap):
#   plain    k-medians, no HMM, no cloner
#   oracle   oracle bin clusters from the truth, no cloner
#   hmmK     HMM with K = SUGGEST_K clusters, no cloner
#   hmmfind  HMM, K chosen by BIC in FIND_RANGE, no cloner
#   cloner   HMM, K chosen by BIC in FIND_RANGE, with the cloner
# Run it on the LOGIN node, not through sbatch:   bash run_detection_sweep.sh   (UGP bins must exist already:
# run_sim_ugp_bins.sh). Override with env vars, e.g.  BINS="16 32" SUGGEST_K=36 FIND_RANGE=24-48 bash run_detection_sweep.sh
# --dry-run is passed on if you give it:  bash run_detection_sweep.sh --dry-run
set -euo pipefail
cd "$(dirname "$(readlink -f "$0")")"
BINS=(${BINS:-4 8 16 32 48 64 80})
SUGGEST_K=${SUGGEST_K:-36}          # 35 event-bearing clades + the normal state
FIND_RANGE=${FIND_RANGE:-24-48}
EXTRA=("$@")
for n in "${BINS[@]}"; do
    C=(./run_pipeline.sh --phaser ugp --caller chisel --phaserbins "${n}" --from call --to eval "${EXTRA[@]+"${EXTRA[@]}"}")
    PLAIN=$("${C[@]}")
    echo "bins ${n}: plain+eval ${PLAIN}" >&2
    for variant in "--oracle" "--hmm ${SUGGEST_K}" "--hmm ${FIND_RANGE}" "--hmm ${FIND_RANGE} --cloner"; do
        # shellcheck disable=SC2086
        J=$("${C[@]}" ${variant} --after "${PLAIN}")
        echo "bins ${n}: ${variant} ${J}" >&2
    done
done
