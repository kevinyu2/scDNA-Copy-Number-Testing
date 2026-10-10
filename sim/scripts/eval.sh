#!/bin/bash
#SBATCH --job-name=eval
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_%j.out
#SBATCH --error=logs/eval_%j.err

# Evaluate one phaser -> caller route against the simulator's true copy numbers.
# Submitted by run_pipeline.sh (PHASER and CALLER come from the environment).
#
# 1. Align each HG002 haplotype (mat, pat) to hg38 with minimap2 (once; cached
#    in ${EVAL_LIFT_DIR}, reused by every project that uses the same FASTAs).
# 2. lift_truth.py: move cn_mat truth from haplotype coordinates to hg38
#    (${EVAL_DIR}/truth_hg38.tsv, shared by all routes).
# 3. standardize_calls.py: the caller's output -> common calls table.
# 4. eval_calls.py: per-cell track plots, pooled confusion matrices, tables.
#
# Output: ${EVAL_DIR}/<phaser>_<caller>_<binsize bp>/   e.g. eval/ugp_chisel_5000000/
# All settings live in config.sh (EVAL section).

LOG_TAG=eval
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
[[ -n ${PHASER:-} && -n ${CALLER:-} ]] || die "PHASER and CALLER must be set (submit through run_pipeline.sh)"
route_ok "${PHASER}" "${CALLER}" || die "route ${PHASER} -> ${CALLER} is not implemented"

source ${CONDA}
conda activate ${EVAL_ENV}

set -euo pipefail

SCRIPTS="$(pwd)/scripts"
CALL_OUT=$(caller_output "${CALLER}" "${PHASER}")
OUT_DIR=$(eval_out_dir "${PHASER}" "${CALLER}")
TRUTH="${EVAL_TRUTH}"
CALLS_STD="${OUT_DIR}/calls.tsv"
# LEVELS: which haplotype orientations are meaningful for the caller.
#   all              actual / ideal phasing / allele-specific (callers that phase across bins)
#   allele_specific  HiScanner: its CN_A|CN_B is always major|minor (it mirrors BAF per bin),
#                    so only sorted (allele-specific) scoring is meaningful
LEVELS=all
case ${CALLER} in
    hiscanner) LABEL="HiScanner (${PHASER}, $(caller_binsize hiscanner) bp bins)"; LEVELS=allele_specific ;;
    chisel)    CHISEL_EVAL_CN="${CHISEL_EVAL_CN:-corrected}"
               LABEL="CHISEL (${PHASER}, $(caller_binsize chisel) bp bins, ${CHISEL_EVAL_CN})" ;;
    *)         LABEL="${CALLER} (${PHASER})" ;;
esac

for tool in minimap2 samtools python3; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${EVAL_ENV}"
done
python3 -c "import pandas, numpy, matplotlib" 2>/dev/null || die "python needs pandas, numpy, matplotlib in ${EVAL_ENV}"
for f in "${SCRIPTS}/lift_truth.py" "${SCRIPTS}/standardize_calls.py" "${SCRIPTS}/eval_calls.py" "${SIM_MAT_FA}" "${SIM_PAT_FA}" "${REF}"; do
    [[ -e $f ]] || die "missing $f"
done
ls "${CN_MAT_DIR}"/*mat.tsv >/dev/null 2>&1 || die "no truth files in ${CN_MAT_DIR}"
[[ -e ${CALL_OUT} ]] || die "missing ${CALL_OUT} (run the call stage first)"

mkdir -p "${EVAL_LIFT_DIR}" "${OUT_DIR}"

# ============================================================
# 1. Haplotype -> hg38 alignments (cached)
# ============================================================

make_lift_pafs      # scripts/common.sh: reuses the cached PAFs when they are newer than the FASTAs

# Shared files are written to a temp name unique to this job, then renamed into place
# (atomic), so evals running at the same time never write into the same file.
TMP_TAG="tmp.${SLURM_JOB_ID:-nojob}.$(hostname -s).$$"

# ============================================================
# 2. Lift truth to hg38
# ============================================================

# Reused when it is newer than every input and was built with the same settings
# (recorded in truth_hg38.tsv.stamp); otherwise rebuilt. Shared by all routes.
MAT_PAF=$(lift_paf mat)
PAT_PAF=$(lift_paf pat)
STAMP="${TRUTH}.stamp"
STAMP_TEXT="chr=${CHR} merge_gap=${EVAL_MERGE_GAP:-1000} cn_mat=${CN_MAT_DIR} mat_paf=${MAT_PAF} pat_paf=${PAT_PAF}"

truth_is_current() {
    [[ -s ${TRUTH} && -s ${STAMP} ]] || return 1
    [[ $(cat "${STAMP}") == "${STAMP_TEXT}" ]] || return 1
    # any input newer than the truth file -> stale
    [[ -z $(find "${CN_MAT_DIR}" "${MAT_PAF}" "${PAT_PAF}" "${SCRIPTS}/lift_truth.py" \
               -newer "${TRUTH}" -print -quit) ]]
}

if truth_is_current; then
    note "Reusing ${TRUTH} (inputs unchanged)"
else
    note "Lifting truth to hg38"
    python3 -B "${SCRIPTS}/lift_truth.py" \
        --cn-mat-dir "${CN_MAT_DIR}" \
        --mat-paf "${MAT_PAF}" \
        --pat-paf "${PAT_PAF}" \
        --query-chrom "${CHR}" \
        --target-chrom "${CHR}" \
        --merge-gap "${EVAL_MERGE_GAP:-1000}" \
        --out "${TRUTH}.${TMP_TAG}"
    echo "${STAMP_TEXT}" > "${STAMP}.${TMP_TAG}"
    # atomic renames: evals of other routes may read these concurrently
    mv "${TRUTH}.${TMP_TAG}" "${TRUTH}"
    mv "${STAMP}.${TMP_TAG}" "${STAMP}"
fi

# ============================================================
# 3. Caller output -> common calls table
# ============================================================

case ${CALLER} in
    hiscanner)
        python3 -B "${SCRIPTS}/standardize_calls.py" hiscanner \
            --calls-dir "${CALL_OUT}" --out "${CALLS_STD}" ;;
    chisel)
        python3 -B "${SCRIPTS}/standardize_calls.py" chisel \
            --calls "${CALL_OUT}" --barcodes "${CHISEL_BARCODES}" \
            --clones "$(call_dir chisel "${PHASER}")/clones/mapping.tsv" \
            --cn "${CHISEL_EVAL_CN}" --out "${CALLS_STD}" ;;
esac

# ============================================================
# 4. Evaluate
# ============================================================

python3 -B "${SCRIPTS}/eval_calls.py" \
    --truth "${TRUTH}" \
    --calls "${CALLS_STD}" \
    --caller-name "${LABEL}" \
    --out-dir "${OUT_DIR}" \
    --purity "${EVAL_PURITY}" \
    --min-covered "${EVAL_MIN_COVERED:-0.5}" \
    --y-linear-max "${EVAL_Y_LINEAR_MAX:-8}" \
    --merge-gap "${EVAL_MERGE_GAP:-1000}" \
    --max-mean-dev "${EVAL_MAX_MEAN_DEV:-0.5}" \
    --orientation "${EVAL_ORIENTATION}" \
    --cells-per-page "${EVAL_CELLS_PER_PAGE}" \
    --max-cells-plot "${EVAL_MAX_CELLS_PLOT}" \
    --conf-unit "${EVAL_CONF_UNIT:-bp}" \
    --levels "${LEVELS}"

note "Done: ${OUT_DIR}"