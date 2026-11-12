#!/bin/bash
#SBATCH --job-name=call_chisel_cellhmm
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=24:00:00
#SBATCH --output=logs/call_chisel_cellhmm_%j.out
#SBATCH --error=logs/call_chisel_cellhmm_%j.err

# The per-cell HMM (run_pipeline.sh --caller chisel --cellhmm K) on the combo.tsv of an existing CHISEL run:
# chisel_hmm's chisel/CellHMM.py, one HMM over every cell and bin emitting (RDR, BAF), K states shared by all
# cells with learned means started from the copy-number grid, fixed transitions (CHISEL_CELLHMM_STAY), uniform
# start. Every cell-bin gets a state, labelled with a copy number a|b. No cloner: the eval scores the per-cell calls.
#
# Input:   $(chisel_base_dir PHASER)/combo/combo.tsv   the same run without --cellhmm (e.g. chisel_pbins8/ugp/)
# Outputs (in ${PROJECT_DIR}/calls/chisel_<bins>_cellhmm<K>/${PHASER}/, e.g. chisel_pbins8_cellhmm12/ugp/):
#   combo -> the source run's combo/ (a link)       source_run.txt   which run and options
#   calls/calls.tsv     CHISEL calls.tsv columns, CLUSTER = state, CN_STATE = the state's a|b
#   calls/cellhmm.tsv   EM history, BIC, per state: label, initial grid point, learned RDR and BAF means, share
#   calls/cells.tsv     per cell: mean copy number, log-likelihood, noise, fraction of bins not 1|1
#   calls/log           the run's messages
#   calls/cellhmm_states.pdf / .png   RDR x BAF of the cell-bins against the states' means, single cells along the
#                       genome (scripts/plot_cellhmm.py, run in EVAL_ENV; a failure there does not fail the job)

LOG_TAG=call_chisel_cellhmm
SCRIPTS="$(pwd)/scripts"      # jobs run from sim/
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-ugp}

[[ -n ${CHISEL_CELLHMM:-} ]] || die "no number of states (--cellhmm K)"
source ${CONDA}

set -euo pipefail
conda_switch() { set +u; conda "$@"; local rc=$?; set -u; return ${rc}; }

BASE=$(chisel_base_dir "${PHASER}")
RUN_DIR=$(call_dir chisel "${PHASER}")
COMBO="${BASE}/combo/combo.tsv"
[[ -s ${COMBO} ]] || die "missing ${COMBO}: run the same command without --cellhmm first"

conda_switch activate "${HMM_CHISEL_ENV:?HMM_CHISEL_ENV is not set in config.sh}"
command -v python2.7 >/dev/null || die "'python2.7' not found in env ${HMM_CHISEL_ENV}"
python2.7 -c "import chisel.CellHMM" 2>/dev/null \
    || die "env ${HMM_CHISEL_ENV}: no chisel.CellHMM (apply the chisel_hmm cellhmm patch to the fork's editable install)"
python2.7 -c "import pandas" 2>/dev/null \
    || note "WARNING: no pandas in ${HMM_CHISEL_ENV}: combo.tsv is read line by line (slow, much more memory)"
# shellcheck disable=SC2206
EXTRA=(${CHISEL_CELLHMM_EXTRA:-})

rm -rf "${RUN_DIR:?}/calls"
mkdir -p "${RUN_DIR}/calls"
rel_link "${BASE}/combo" "${RUN_DIR}/combo"
{ echo "source run: ${BASE}"; echo "combo.tsv:  ${COMBO}"
  echo "cellhmm: -K ${CHISEL_CELLHMM} --stay ${CHISEL_CELLHMM_STAY:-0.99} --seed ${CHISEL_SEED} ${EXTRA[*]}"
  echo "date: $(date)"; } > "${RUN_DIR}/source_run.txt"

note "Source run: ${BASE}"
note "Run:        ${RUN_DIR}  (env ${HMM_CHISEL_ENV})"
note "python2.7 -m chisel.CellHMM combo/combo.tsv -K ${CHISEL_CELLHMM} --stay ${CHISEL_CELLHMM_STAY:-0.99} --seed ${CHISEL_SEED} ${EXTRA[*]}"
cd "${RUN_DIR}"
python2.7 -B -m chisel.CellHMM combo/combo.tsv -K "${CHISEL_CELLHMM}" --stay "${CHISEL_CELLHMM_STAY:-0.99}" \
    --seed "${CHISEL_SEED}" -o calls/calls.tsv --cells calls/cells.tsv --report calls/cellhmm.tsv \
    "${EXTRA[@]+"${EXTRA[@]}"}" 2> >(tee calls/log >&2) \
    || { echo "ERROR: CellHMM failed; last lines of calls/log:" >&2; tail -n 20 calls/log >&2; exit 1; }

[[ $(grep -vc '^#' calls/calls.tsv 2>/dev/null || true) -gt 0 ]] || die "no calls in calls/calls.tsv"
note "States (calls/cellhmm.tsv):"
sed 's/^/    /' calls/cellhmm.tsv >&2
conda_switch activate "${EVAL_ENV}" && python3 -B "${SCRIPTS}/plot_cellhmm.py" --run . >&2 \
    || note "WARNING: plot_cellhmm.py failed (calls are fine); rerun: python3 scripts/plot_cellhmm.py --run ${RUN_DIR}"
note "Done: ${RUN_DIR}"
