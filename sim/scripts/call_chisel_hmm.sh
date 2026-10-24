#!/bin/bash
#SBATCH --job-name=call_chisel_hmm
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --output=logs/call_chisel_hmm_%j.out
#SBATCH --error=logs/call_chisel_hmm_%j.err

# CHISEL's calling step again, with HMM bin clustering instead of k-medians, on the combo.tsv of an
# existing CHISEL run. Submitted by run_pipeline.sh --caller chisel --hmm K (or MIN-MAX).
#
# Input:   $(chisel_base_dir PHASER)/combo/combo.tsv   the same run without --hmm (e.g. chisel_pbins8/perfect/)
# Outputs (in ${PROJECT_DIR}/chisel_<bins>_hmm<K>/${PHASER}/, e.g. chisel_pbins8_hmm7/perfect/):
#   combo -> the source run's combo/ (a link)       source_run.txt   which run it came from
#   calls/calls.tsv  calls/hmm_selection.tsv (per-K table, * = chosen)   clones/mapping.tsv   plots/
# Runs chisel_calling (Caller --clustering hmm, Cloner, Plotter) from HMM_CHISEL_ENV with the run's
# CHISEL settings (CHISEL_MAXPLOIDY, CHISEL_UPPERK, CHISEL_SEED, CHISEL_JOBS); cloning is
# length-weighted for variable-length bins, as chisel does.

LOG_TAG=call_chisel_hmm
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-perfect}

[[ -n ${CHISEL_HMM:-} ]] || die "no HMM clusters given (run_pipeline.sh --hmm K)"
source ${CONDA}
conda activate "${HMM_CHISEL_ENV:?HMM_CHISEL_ENV is not set in config.sh}"

set -euo pipefail

BASE=$(chisel_base_dir "${PHASER}")
RUN_DIR=$(call_dir chisel "${PHASER}")
COMBO="${BASE}/combo/combo.tsv"
[[ -s ${COMBO} ]] || die "missing ${COMBO}: run the same command without --hmm first"
for tool in chisel_calling python2.7; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${HMM_CHISEL_ENV}"
done
chisel_calling -h 2>/dev/null | grep -q -- "--hmmk" || die "env ${HMM_CHISEL_ENV} has no HMM chisel_calling (needs the chisel_hmm hmm branch with chisel_calling --hmm)"

if [[ ${CHISEL_HMM} == *-* ]]; then
    HMM_ARGS=(--hmmmink "${CHISEL_HMM%-*}" --hmmmaxk "${CHISEL_HMM#*-}" --hmmselect "${CHISEL_HMM_SELECT:-bic}")
else
    HMM_ARGS=(--hmmk "${CHISEL_HMM}")
fi
WEIGHTED=()
[[ ${CHISEL_VARBINS:-false} == true || ${CHISEL_PHASER_BINS:-false} == true ]] && WEIGHTED=(--weighted)
# shellcheck disable=SC2206
EXTRA=(${CHISEL_HMM_EXTRA:-})

mkdir -p "${RUN_DIR}"
for d in calls clones plots; do rm -rf "${RUN_DIR:?}/${d}"; done
ln -sfn "${BASE}/combo" "${RUN_DIR}/combo"
{ echo "source run: ${BASE}"; echo "combo.tsv:  ${COMBO}"; echo "HMM: ${CHISEL_HMM} ${HMM_ARGS[*]} ${EXTRA[*]}"; echo "date: $(date)"; } > "${RUN_DIR}/source_run.txt"

note "Source run: ${BASE}"
note "HMM run:    ${RUN_DIR}  (env ${HMM_CHISEL_ENV})"
note "chisel_calling ... --hmm ${HMM_ARGS[*]} --hmmrestarts ${CHISEL_HMM_RESTARTS:-10} ${WEIGHTED[*]} ${EXTRA[*]}"
cd "${RUN_DIR}"
chisel_calling combo/combo.tsv -x . \
    -P "${CHISEL_MAXPLOIDY:-3}" -K "${CHISEL_UPPERK:-100}" --seed "${CHISEL_SEED}" -j "${CHISEL_JOBS}" \
    "${WEIGHTED[@]+"${WEIGHTED[@]}"}" \
    --hmm "${HMM_ARGS[@]}" --hmmrestarts "${CHISEL_HMM_RESTARTS:-10}" \
    "${EXTRA[@]+"${EXTRA[@]}"}" \
    || true      # CHISEL's drivers do not report failures reliably; decide from the outputs

fail() { echo "ERROR: $1" >&2; [[ -s $2 ]] && { echo "--- last lines of $2 ---" >&2; tail -n 20 "$2" >&2; }; exit 1; }
[[ $(grep -vc '^#' calls/calls.tsv 2>/dev/null || true) -gt 0 ]] || fail "no calls in calls/calls.tsv" calls/log
[[ -s calls/hmm_selection.tsv ]] || fail "no calls/hmm_selection.tsv (did the HMM run?)" calls/log
[[ -s clones/mapping.tsv ]] || fail "no clones/mapping.tsv" clones/log
note "Model selection (calls/hmm_selection.tsv):"
sed 's/^/    /' calls/hmm_selection.tsv >&2
note "Bin clusters in calls.tsv: $(grep -v '^#' calls/calls.tsv | cut -f11 | sort -u | wc -l)"
note "Done: ${RUN_DIR}"
