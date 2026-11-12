#!/bin/bash
#SBATCH --job-name=call_chisel_bincellhmm
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/call_chisel_bincellhmm_%j.out
#SBATCH --error=logs/call_chisel_bincellhmm_%j.err

# RunpengLuo's bin-cell HMM (run_pipeline.sh --caller chisel --bincellhmm K) on the combo.tsv of an existing CHISEL
# run: bin_cell_hmm/hmm_caller.py from a clone of github.com/RunpengLuo/chisel, branch feat/bin_cell_hmm
# (BINCELLHMM_REPO). One HMM over every cell and bin, K (RDR, BAF) states shared by all cells, plain diagonal
# Gaussians, start from the copy-number grid by constrained k-means, fixed transitions (CHISEL_CELLHMM_STAY), uniform
# start, Viterbi; then CHISEL's own per-cell calling (that clone's src/chisel/Caller.py: base cluster, ploidy by BIC
# up to CHISEL_MAXPLOIDY) gives every cell's states copy numbers. No cloner: the eval scores the per-cell calls.
#
# Input:   $(chisel_base_dir PHASER)/combo/combo.tsv   the same run without --bincellhmm
# Outputs (in ${PROJECT_DIR}/calls/chisel_<bins>_bincellhmm<K>/${PHASER}/):
#   combo -> the source run's combo/ (a link)       source_run.txt   which run, commit and options
#   calls/calls.tsv        CHISEL calls.tsv columns, CLUSTER = HMM state, CN_STATE = the cell's call
#   calls/bincellhmm.npz   the fitted HMM: means, variances, state_cn, init_states, members, loglik, ploidy
#   calls/log              the run's messages
#   calls/cellhmm_states.pdf / .png   RDR x BAF against the states' means, single cells (scripts/plot_cellhmm.py)

LOG_TAG=call_chisel_bincellhmm
SCRIPTS="$(pwd)/scripts"      # jobs run from sim/
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-ugp}

[[ -n ${CHISEL_BINCELLHMM:-} ]] || die "no number of states (--bincellhmm K)"
CALLER_PY="${BINCELLHMM_REPO:-}/bin_cell_hmm/hmm_caller.py"
[[ -s ${CALLER_PY} ]] || die "no ${CALLER_PY}: git clone -b feat/bin_cell_hmm https://github.com/RunpengLuo/chisel.git ${BINCELLHMM_REPO:-<BINCELLHMM_REPO>}"
source ${CONDA}

set -euo pipefail
conda_switch() { set +u; conda "$@"; local rc=$?; set -u; return ${rc}; }

BASE=$(chisel_base_dir "${PHASER}")
RUN_DIR=$(call_dir chisel "${PHASER}")
COMBO="${BASE}/combo/combo.tsv"
[[ -s ${COMBO} ]] || die "missing ${COMBO}: run the same command without --bincellhmm first"

conda_switch activate "${HMM_CHISEL_ENV:?HMM_CHISEL_ENV is not set in config.sh}"
command -v python2.7 >/dev/null || die "'python2.7' not found in env ${HMM_CHISEL_ENV}"
# hmm.py needs numpy >= 1.15 (argsort kind="stable"); CHISEL's Caller needs scipy
python2.7 -c "import numpy as np, scipy.stats; np.argsort([1, 0], kind='stable')" 2>/dev/null \
    || die "env ${HMM_CHISEL_ENV}: needs scipy and numpy >= 1.15 for bin_cell_hmm/hmm.py"
# shellcheck disable=SC2206
EXTRA=(${CHISEL_BINCELLHMM_EXTRA:-})
COMMIT=$(git -C "${BINCELLHMM_REPO}" log -1 --format='%h %ci' 2>/dev/null || echo unknown)

rm -rf "${RUN_DIR:?}/calls"
mkdir -p "${RUN_DIR}/calls"
rel_link "${BASE}/combo" "${RUN_DIR}/combo"
ARGS=(--hmm-states "${CHISEL_BINCELLHMM}" --hmm-stay "${CHISEL_CELLHMM_STAY:-0.99}" --hmm-npz calls/bincellhmm.npz
      -P "${CHISEL_MAXPLOIDY:-3}" -j "${CHISEL_JOBS:-16}" --seed "${CHISEL_SEED}")
{ echo "source run: ${BASE}"; echo "combo.tsv:  ${COMBO}"
  echo "bin_cell_hmm: ${CALLER_PY} (commit ${COMMIT})"
  echo "options: ${ARGS[*]} ${EXTRA[*]}"
  echo "date: $(date)"; } > "${RUN_DIR}/source_run.txt"

note "Source run: ${BASE}"
note "Run:        ${RUN_DIR}  (env ${HMM_CHISEL_ENV}, bin_cell_hmm ${COMMIT})"
note "python2.7 hmm_caller.py combo/combo.tsv ${ARGS[*]} ${EXTRA[*]}"
cd "${RUN_DIR}"
PYTHONDONTWRITEBYTECODE=1 python2.7 "${CALLER_PY}" combo/combo.tsv "${ARGS[@]}" "${EXTRA[@]+"${EXTRA[@]}"}" \
    > calls/calls.tsv 2> >(tee calls/log >&2) \
    || { echo "ERROR: hmm_caller.py failed; last lines of calls/log:" >&2; tail -n 20 calls/log >&2; exit 1; }

[[ $(grep -vc '^#' calls/calls.tsv 2>/dev/null || true) -gt 0 ]] || die "no calls in calls/calls.tsv"
conda_switch activate "${EVAL_ENV}" && python3 -B "${SCRIPTS}/plot_cellhmm.py" --run . >&2 \
    || note "WARNING: plot_cellhmm.py failed (calls are fine); rerun: python3 scripts/plot_cellhmm.py --run ${RUN_DIR}"
note "Done: ${RUN_DIR}"
