#!/bin/bash
#SBATCH --job-name=scDNA_sim_eval
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=12:00:00
#SBATCH --output=logs/eval_scDNA_sim_%j.out
#SBATCH --error=logs/eval_scDNA_sim_%j.err

# Evaluate HiScanner against the simulator's true copy numbers.
#
# 1. Align each HG002 haplotype (mat, pat) to hg38 with minimap2 (once; cached
#    in ${EVAL_LIFT_DIR}, reused by every project that uses the same FASTAs).
# 2. lift_truth.py: move cn_mat truth from haplotype coordinates to hg38.
# 3. eval_hiscanner.py: per-cell track plots, pooled confusion matrices, tables.
#
# All settings live in config.sh (EVAL section).

source "${CONFIG_FILE:-./config.sh}"

source ${CONDA}
conda activate ${EVAL_ENV}

set -euo pipefail

SCRIPTS="$(pwd)/scripts"
CN_MAT_DIR="${PROJECT_DIR}/sim/cn_mat"
CALLS_DIR="${HS_DIR}/output/final_calls"
OUT_DIR="${EVAL_DIR}/hiscanner"
TRUTH="${EVAL_DIR}/truth_hg38.tsv"

die()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[eval_run] $*" >&2; }

for tool in minimap2 samtools python3; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${EVAL_ENV}"
done
python3 -c "import pandas, numpy, matplotlib" 2>/dev/null || die "python needs pandas, numpy, matplotlib in ${EVAL_ENV}"
for f in "${SCRIPTS}/lift_truth.py" "${SCRIPTS}/eval_hiscanner.py" "${SIM_MAT_FA}" "${SIM_PAT_FA}" "${REF}"; do
    [[ -e $f ]] || die "missing $f"
done
ls "${CN_MAT_DIR}"/*mat.tsv >/dev/null 2>&1 || die "no truth files in ${CN_MAT_DIR}"
[[ -d ${CALLS_DIR} ]] || die "missing ${CALLS_DIR} (run HiScanner first)"

mkdir -p "${EVAL_LIFT_DIR}" "${OUT_DIR}"

# ============================================================
# 1. Haplotype -> hg38 alignments (cached)
# ============================================================

HG38_CHR="${EVAL_LIFT_DIR}/hg38.${CHR}.fa"
if [[ ! -s ${HG38_CHR} ]]; then
    note "Extracting ${CHR} from ${REF}"
    samtools faidx "${REF}" "${CHR}" > "${HG38_CHR}.tmp" && mv "${HG38_CHR}.tmp" "${HG38_CHR}"
fi

for hap in mat pat; do
    if [[ ${hap} == mat ]]; then FA="${SIM_MAT_FA}"; else FA="${SIM_PAT_FA}"; fi
    PAF="${EVAL_LIFT_DIR}/${hap}_to_hg38.${CHR}.paf"
    if [[ -s ${PAF} && ${PAF} -nt ${FA} ]]; then
        note "Reusing ${PAF}"
        continue
    fi
    note "Aligning ${hap} haplotype to hg38 ${CHR} (minimap2 asm5; tens of minutes)"
    minimap2 -c -x asm5 --cs --secondary=no -t "${SLURM_CPUS_PER_TASK:-8}" \
        "${HG38_CHR}" "${FA}" > "${PAF}.tmp"
    mv "${PAF}.tmp" "${PAF}"
done

# ============================================================
# 2. Lift truth to hg38
# ============================================================

python3 -B "${SCRIPTS}/lift_truth.py" \
    --cn-mat-dir "${CN_MAT_DIR}" \
    --mat-paf "${EVAL_LIFT_DIR}/mat_to_hg38.${CHR}.paf" \
    --pat-paf "${EVAL_LIFT_DIR}/pat_to_hg38.${CHR}.paf" \
    --query-chrom "${CHR}" \
    --target-chrom "${CHR}" \
    --merge-gap "${EVAL_MERGE_GAP:-1000}" \
    --out "${TRUTH}"

# ============================================================
# 3. Evaluate HiScanner
# ============================================================

python3 -B "${SCRIPTS}/eval_hiscanner.py" \
    --truth "${TRUTH}" \
    --calls-dir "${CALLS_DIR}" \
    --out-dir "${OUT_DIR}" \
    --purity "${EVAL_PURITY}" \
    --min-covered "${EVAL_MIN_COVERED:-0.5}" \
    --y-linear-max "${EVAL_Y_LINEAR_MAX:-8}" \
    --merge-gap "${EVAL_MERGE_GAP:-1000}" \
    --orientation "${EVAL_ORIENTATION}" \
    --cells-per-page "${EVAL_CELLS_PER_PAGE}" \
    --max-cells-plot "${EVAL_MAX_CELLS_PLOT}" \
    --max-cn "${EVAL_MAX_CN}"

note "Done: ${OUT_DIR}"