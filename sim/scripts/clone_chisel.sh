#!/bin/bash
#SBATCH --job-name=clone_chisel
#SBATCH --mem=32G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --output=logs/clone_chisel_%j.out
#SBATCH --error=logs/clone_chisel_%j.err

# CHISEL's cloning step again, on the calls of an existing CHISEL run, with another way of grouping
# cells into clones (run_pipeline.sh --caller chisel --cloning divisive): the chisel_hmm fork's
# Cloner --clustering divisive (top-down splitting on the cells' RDR and BAF, permutation-tested;
# see the fork's doc/chisel-divisive-cloning.md). Works on any CHISEL run, also --hmm and --oracle
# runs: the source is the same run with CHISEL's own cloning.
#
# Input:   $(clone_base_dir PHASER)/calls/calls.tsv   e.g. chisel_pbins8_hmm9/perfect/calls/calls.tsv
# Outputs (in ${PROJECT_DIR}/chisel_<bins>[_hmm<K>|_oracle]_divclone/${PHASER}/):
#   calls/calls.tsv   a copy of the source's calls with CORRECTED_HAP_CN redone by the new clones
#   clones/mapping.tsv  clones/log (the split tree)   plots/   combo -> the source's combo/ (a link)
#   source_run.txt
# The source run is not modified.

LOG_TAG=clone_chisel
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-perfect}

[[ ${CHISEL_CLONING:-hamming} != hamming ]] || die "no cloning method given (run_pipeline.sh --cloning divisive)"
source ${CONDA}

set -euo pipefail
# conda's activate/deactivate scripts read unset variables (e.g. binutils' ADDR2LINE): not under set -u
conda_switch() { set +u; conda "$@"; local rc=$?; set -u; return ${rc}; }
fail() { echo "ERROR: $1" >&2; [[ -s ${2:-} ]] && { echo "--- last lines of $2 ---" >&2; tail -n 20 "$2" >&2; }; exit 1; }

BASE=$(clone_base_dir "${PHASER}")
RUN_DIR=$(call_dir chisel "${PHASER}")
[[ ${BASE} != "${RUN_DIR}" ]] || die "source and output are the same folder: ${RUN_DIR}"
CALLS="${BASE}/calls/calls.tsv"
[[ -s ${CALLS} ]] || die "missing ${CALLS}: run the same command without --cloning first"

conda_switch activate "${HMM_CHISEL_ENV:?HMM_CHISEL_ENV is not set in config.sh}"
command -v chisel_cloning >/dev/null || die "'chisel_cloning' not found in env ${HMM_CHISEL_ENV}"
chisel_cloning -h 2>/dev/null | grep -q -- "--cloning" \
    || die "env ${HMM_CHISEL_ENV}: chisel_cloning has no --cloning (apply the chisel_hmm patch with Cloner --clustering divisive)"

mkdir -p "${RUN_DIR}"
for d in calls clones plots; do rm -rf "${RUN_DIR:?}/${d}"; done
mkdir -p "${RUN_DIR}/calls"
cut -f1-12 "${CALLS}" > "${RUN_DIR}/calls/calls.tsv"      # the source's calls without its corrected states
[[ -e ${BASE}/combo ]] && ln -sfn "$(readlink -f "${BASE}/combo")" "${RUN_DIR}/combo"
WEIGHTED=()
[[ ${CHISEL_VARBINS:-false} == true || ${CHISEL_PHASER_BINS:-false} == true ]] && WEIGHTED=(--weighted)
{ echo "source run: ${BASE}"; echo "calls:      ${CALLS}"; echo "cloning:    ${CHISEL_CLONING} ${WEIGHTED[*]}"; echo "date: $(date)"; } \
    > "${RUN_DIR}/source_run.txt"

note "Source run: ${BASE}"
note "Run:        ${RUN_DIR}  (env ${HMM_CHISEL_ENV})"
note "chisel_cloning calls/calls.tsv -x . --cloning ${CHISEL_CLONING} --seed ${CHISEL_SEED} ${WEIGHTED[*]}"
cd "${RUN_DIR}"
chisel_cloning calls/calls.tsv -x . --cloning "${CHISEL_CLONING}" --seed "${CHISEL_SEED}" \
    "${WEIGHTED[@]+"${WEIGHTED[@]}"}" \
    || true      # CHISEL's drivers do not report failures reliably; decide from the outputs

[[ -s clones/mapping.tsv ]] || fail "no clones/mapping.tsv" clones/log
head -1 calls/calls.tsv | grep -q CORRECTED || fail "calls/calls.tsv has no corrected states" clones/log
note "Split tree (clones/log):"
sed 's/\x1b\[[0-9;]*m//g' clones/log | sed -n '/very noisy\|Divisive splits/,/reassignment/p' | sed 's/^/    /' >&2 || true
note "Clones: $(grep -v '^#' clones/mapping.tsv | cut -f3 | sort | uniq -c | sort -nr | awk '{printf "%s %s; ", $2, $1}')"
note "Done: ${RUN_DIR}"
