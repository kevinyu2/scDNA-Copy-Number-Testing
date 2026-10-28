#!/bin/bash
#SBATCH --job-name=call_chisel_hmm
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=4:00:00
#SBATCH --output=logs/call_chisel_hmm_%j.out
#SBATCH --error=logs/call_chisel_hmm_%j.err

# CHISEL's calling step again, on the combo.tsv of an existing CHISEL run, with HMM bin clustering
# instead of k-medians (run_pipeline.sh --caller chisel --hmm K or MIN-MAX), or with the oracle bin
# clusters from the simulation truth (--oracle: scripts/oracle_clusters.py writes oracle/clusters.tsv
# and oracle/summary.tsv, then Caller takes them with --clusters).
#
# Input:   $(chisel_base_dir PHASER)/combo/combo.tsv   the same run without --hmm (e.g. chisel_pbins8/perfect/)
# Outputs (in ${PROJECT_DIR}/calls/chisel_<bins>_hmm<K>/${PHASER}/, e.g. chisel_pbins8_hmm7/perfect/):
#   combo -> the source run's combo/ (a link)       source_run.txt   which run it came from
#   calls/calls.tsv  calls/hmm_selection.tsv (per-K table, * = chosen)   clones/mapping.tsv   plots/
# Runs chisel_calling (Caller --clustering hmm, Cloner, Plotter) from HMM_CHISEL_ENV with the run's
# CHISEL settings (CHISEL_MAXPLOIDY, CHISEL_UPPERK, CHISEL_SEED, CHISEL_JOBS); cloning is
# length-weighted for variable-length bins, as chisel does.

LOG_TAG=call_chisel_hmm
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-perfect}

ORACLE=${CHISEL_ORACLE:-false}
SCRIPTS="$(pwd)/scripts"      # jobs run from sim/ (under SLURM $0 is a spool copy)
[[ -n ${CHISEL_HMM:-} || ${ORACLE} == true ]] || die "neither HMM clusters (--hmm K) nor --oracle given"
source ${CONDA}

set -euo pipefail
# conda's activate/deactivate scripts read unset variables (e.g. binutils' ADDR2LINE): not under set -u
conda_switch() { set +u; conda "$@"; local rc=$?; set -u; return ${rc}; }

BASE=$(chisel_base_dir "${PHASER}")
RUN_DIR=$(call_dir chisel "${PHASER}")
COMBO="${BASE}/combo/combo.tsv"
[[ -s ${COMBO} ]] || die "missing ${COMBO}: run the same command without --hmm/--oracle first"
mkdir -p "${RUN_DIR}"
if [[ ${ORACLE} == true ]]; then
    [[ -s ${EVAL_TRUTH} ]] || die "missing ${EVAL_TRUTH}: run any eval of this project first"
    conda_switch activate "${EVAL_ENV}"
    note "Oracle bin clusters from ${EVAL_TRUTH} on the bins of ${BASE}"
    SPLIT_FRAC=$(awk -v p="${CHISEL_ORACLE_SPLIT:-0}" 'BEGIN{print p/100}')
    python3 -B "${SCRIPTS}/oracle_clusters.py" --truth "${EVAL_TRUTH}" --run "${BASE}" --out "${RUN_DIR}/oracle" \
        --split-frac "${SPLIT_FRAC}" \
        || die "oracle_clusters.py failed"
    conda_switch deactivate
fi
conda_switch activate "${HMM_CHISEL_ENV:?HMM_CHISEL_ENV is not set in config.sh}"
for tool in chisel_calling python2.7; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${HMM_CHISEL_ENV}"
done
chisel_calling -h 2>/dev/null | grep -q -- "--hmmk" || die "env ${HMM_CHISEL_ENV} has no HMM chisel_calling (needs the chisel_hmm hmm branch with chisel_calling --hmm)"
MODE_ARGS=($(fork_mode_flags))
if [[ ${CHISEL_MODE:-pseudobulk} != legacy ]]; then
    chisel_calling -h 2>/dev/null | grep -q -- "--nopseudobulk" \
        || die "env ${HMM_CHISEL_ENV}: chisel_calling has no --nopseudobulk (apply the chisel_hmm patch with the pseudobulk defaults, or --mode legacy)"
fi
if [[ ${ORACLE} == true ]]; then
    chisel_calling -h 2>/dev/null | grep -q -- "--clusters" || die "env ${HMM_CHISEL_ENV}: chisel_calling has no --clusters (apply the chisel_hmm patch with Caller --clusters)"
    HMM_ARGS=(--clusters "${RUN_DIR}/oracle/clusters.tsv")
elif [[ ${CHISEL_HMM} == *-* ]]; then
    HMM_ARGS=(--hmmmink "${CHISEL_HMM%-*}" --hmmmaxk "${CHISEL_HMM#*-}" --hmmselect "${CHISEL_HMM_SELECT:-bic}")
else
    HMM_ARGS=(--hmmk "${CHISEL_HMM}")
fi
WEIGHTED=()
[[ ${CHISEL_VARBINS:-false} == true || ${CHISEL_PHASER_BINS:-false} == true ]] && WEIGHTED=(--weighted)
# shellcheck disable=SC2206
EXTRA=(${CHISEL_HMM_EXTRA:-})

for d in calls clones plots; do rm -rf "${RUN_DIR:?}/${d}"; done
rel_link "${BASE}/combo" "${RUN_DIR}/combo"
{ echo "source run: ${BASE}"; echo "combo.tsv:  ${COMBO}"; echo "clusters: ${HMM_ARGS[*]} ${EXTRA[*]}"; echo "mode: ${CHISEL_MODE:-pseudobulk} ${MODE_ARGS[*]}"; echo "date: $(date)"; } > "${RUN_DIR}/source_run.txt"

note "Source run: ${BASE}"
note "Run:        ${RUN_DIR}  (env ${HMM_CHISEL_ENV})"
if [[ ${ORACLE} == true ]]; then
    CLUSTER_ARGS=("${HMM_ARGS[@]}")
    EXTRA=()
else
    CLUSTER_ARGS=(--hmm "${HMM_ARGS[@]}" --hmmrestarts "${CHISEL_HMM_RESTARTS:-10}")
fi
note "chisel_calling ... ${CLUSTER_ARGS[*]} ${MODE_ARGS[*]} ${WEIGHTED[*]} ${EXTRA[*]}"
cd "${RUN_DIR}"
chisel_calling combo/combo.tsv -x . \
    -P "${CHISEL_MAXPLOIDY:-3}" -K "${CHISEL_UPPERK:-100}" --seed "${CHISEL_SEED}" -j "${CHISEL_JOBS}" \
    "${WEIGHTED[@]+"${WEIGHTED[@]}"}" \
    "${CLUSTER_ARGS[@]}" \
    "${MODE_ARGS[@]+"${MODE_ARGS[@]}"}" \
    "${EXTRA[@]+"${EXTRA[@]}"}" \
    || true      # CHISEL's drivers do not report failures reliably; decide from the outputs

fail() { echo "ERROR: $1" >&2; [[ -s $2 ]] && { echo "--- last lines of $2 ---" >&2; tail -n 20 "$2" >&2; }; exit 1; }
[[ $(grep -vc '^#' calls/calls.tsv 2>/dev/null || true) -gt 0 ]] || fail "no calls in calls/calls.tsv" calls/log
[[ -s clones/mapping.tsv ]] || fail "no clones/mapping.tsv" clones/log
if [[ ${ORACLE} == true ]]; then
    note "Oracle clusters (oracle/summary.tsv):"
    cut -f1-7 oracle/summary.tsv | sed 's/^/    /' >&2
else
    [[ -s calls/hmm_selection.tsv ]] || fail "no calls/hmm_selection.tsv (did the HMM run?)" calls/log
    note "Model selection (calls/hmm_selection.tsv):"
    sed 's/^/    /' calls/hmm_selection.tsv >&2
fi
note "Bin clusters in calls.tsv: $(grep -v '^#' calls/calls.tsv | cut -f11 | sort -u | wc -l)"
note "Done: ${RUN_DIR}"
