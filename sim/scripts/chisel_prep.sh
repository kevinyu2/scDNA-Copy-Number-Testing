#!/bin/bash
#SBATCH --job-name=chisel_prep
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/chisel_prep_%j.out
#SBATCH --error=logs/chisel_prep_%j.err

# Submitted by run_pipeline.sh --caller chisel (alongside the phaser job; it
# only needs the cell BAMs). Skips itself if the barcoded BAM already exists.
#
# Convert the per-cell BAMs in ${PROJECT_DIR}/sim/bams into the single
# barcoded BAM that CHISEL expects, using CHISEL's own `chisel_prep`.
#
# Outputs (in ${PROJECT_DIR}/chisel/prep, shared by every phaser):
#   barcodedcells.bam(.bai)   all cells, each read tagged with its cell barcode
#   barcodedcells.info.tsv    CELL -> BARCODE map (use this to match CHISEL
#                             results back to your simulated cells)
#   inputs.tsv                the BAM -> cell-name table given to chisel_prep
#
# All settings live in config.sh (CHISEL section).

LOG_TAG=chisel_prep
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh

source ${CONDA}
conda activate ${CHISEL_ENV}

set -euo pipefail

BAM_DIR="${CELL_BAM_DIR}"
PREP_DIR="${CHISEL_PREP_DIR}"
OUT_BAM="barcodedcells.bam"


# ============================================================
# Checks
# ============================================================

for tool in chisel_prep samtools bcftools; do
    # chisel_prep refuses to start without bcftools on PATH, even for BAM input
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${CHISEL_ENV}"
done

mapfile -t BAMS < <(ls "${BAM_DIR}"/*.bam | LC_ALL=C sort)
N=${#BAMS[@]}
(( N > 0 )) || die "no BAMs in ${BAM_DIR}"
for b in "${BAMS[@]}"; do
    [[ -e ${b}.bai || -e ${b%.bam}.bai ]] || die "missing index for $b"
done
note "${N} cell BAMs found"

mkdir -p "${PREP_DIR}"
cd "${PREP_DIR}"   # chisel_prep writes -o relative to the working directory

# ============================================================
# Lock: several CHISEL routes (e.g. ugp and perfect) submit their own chisel_prep
# at the same time; only one may build the shared BAM, the others wait and reuse it.
# mkdir is atomic on NFS. A lock left by a job that no longer runs is taken over.
# ============================================================

LOCK="${PREP_DIR}/.chisel_prep.lock"
ME="${SLURM_JOB_ID:-pid$$}"
lock_owner_gone() {   # true when the job named in the lock is no longer queued or running
    local owner out
    owner=$(cat "${LOCK}/owner" 2>/dev/null || true)
    [[ -n ${owner} ]] || return 1                          # just created; owner not written yet
    [[ ${owner} == pid* ]] && { kill -0 "${owner#pid}" 2>/dev/null && return 1 || return 0; }
    out=$(squeue -h -j "${owner}" -o %T 2>&1) || { [[ ${out} == *"Invalid job id"* ]]; return; }
    [[ ! ${out} =~ (PENDING|RUNNING|CONFIGURING|COMPLETING|SUSPENDED|REQUEUE) ]]
}
WAITED=0
until mkdir "${LOCK}" 2>/dev/null; do
    if lock_owner_gone; then
        note "Removing stale lock from job $(cat "${LOCK}/owner" 2>/dev/null)"
        rm -rf "${LOCK}"
        continue
    fi
    (( WAITED % 600 == 0 )) && note "Waiting for job $(cat "${LOCK}/owner" 2>/dev/null || echo '?') to finish building the barcoded BAM"
    sleep 30; WAITED=$(( WAITED + 30 ))
done
echo "${ME}" > "${LOCK}/owner"
trap 'rm -rf "${LOCK}"' EXIT

if [[ -s ${OUT_BAM} && -s ${OUT_BAM}.bai && ${CHISEL_PREP_FORCE:-false} != true ]]; then
    note "${PREP_DIR}/${OUT_BAM} already exists; set CHISEL_PREP_FORCE=true to rebuild"
    exit 0
fi

# chisel_prep refuses to run if a previous (failed) run left these behind
for d in _TMP_CHISEL_PREP _ERR_CHISEL_PREP; do
    if [[ -e $d ]]; then
        note "Removing leftover $d from a previous run"
        rm -rf "$d"
    fi
done
rm -f "${OUT_BAM}" "${OUT_BAM}.bai"

# samtools merge opens every cell BAM at once; ~1000 files can exceed the
# default open-file limit (often 1024)
NEED_FD=$(( N + 256 ))
SOFT_FD=$(ulimit -Sn)
HARD_FD=$(ulimit -Hn)
if [[ ${SOFT_FD} != unlimited ]] && (( SOFT_FD < NEED_FD )); then
    if [[ ${HARD_FD} == unlimited ]] || (( HARD_FD >= NEED_FD )); then
        ulimit -n "${NEED_FD}"
        note "Raised open-file limit to ${NEED_FD}"
    else
        die "open-file limit (hard ${HARD_FD}) is below ${NEED_FD} needed to merge ${N} BAMs"
    fi
fi

# ============================================================
# Input table: BAM <tab> cell name  (cell name = BAM basename, same as SCAN2/HiScanner)
# ============================================================

{
    printf '#FILE\tCELL\n'
    for b in "${BAMS[@]}"; do
        printf '%s\t%s\n' "$b" "$(basename "$b" .bam)"
    done
} > inputs.tsv

# ============================================================
# Run chisel_prep
#   --noduplicates: dwgsim reads have no PCR duplicates, so skip markdup
# ============================================================

chisel_prep \
    -x "${PREP_DIR}" \
    -o "${OUT_BAM}" \
    --noduplicates \
    --barcodelength "${CHISEL_BARCODE_LENGTH}" \
    --seed "${CHISEL_SEED}" \
    -j "${CHISEL_JOBS}" \
    inputs.tsv

# ============================================================
# Verify
# ============================================================

[[ -s ${OUT_BAM} && -s ${OUT_BAM}.bai ]] || die "chisel_prep finished but ${OUT_BAM} or its index is missing"

N_INFO=$(grep -vc '^#' barcodedcells.info.tsv || true)
(( N_INFO == N )) || die "barcodedcells.info.tsv lists ${N_INFO} cells, expected ${N}"

# CHISEL finds barcodes with the regex CB:Z:[ACGT]+ anywhere in the read line
FIRST=$(samtools view "${OUT_BAM}" | head -n 1000 | grep -m1 -oE 'CB:Z:[ACGT]+' || true)
[[ -n ${FIRST} ]] || die "no CB:Z:<barcode> found in the first reads of ${OUT_BAM}"

N_RG=$(samtools view -H "${OUT_BAM}" | grep -c '^@RG.*CB:Z:' || true)
note "Done: ${PREP_DIR}/${OUT_BAM}"
note "  cells in info.tsv: ${N_INFO}; barcoded read groups in header: ${N_RG}; example tag: ${FIRST}"
note "  cell -> barcode map: ${PREP_DIR}/barcodedcells.info.tsv"