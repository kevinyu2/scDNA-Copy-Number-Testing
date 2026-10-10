#!/bin/bash
#SBATCH --job-name=call_chisel
#SBATCH --mem=64G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/call_chisel_%j.out
#SBATCH --error=logs/call_chisel_%j.err

# CHISEL on the barcoded BAM from chisel_prep.sh, using the phased het SNPs
# from ${PHASER}. Submitted by run_pipeline.sh --caller chisel.
#
# Outputs (in ${PROJECT_DIR}/chisel_<binsize bp>/${PHASER}/, e.g. chisel_5000000/ugp/):
#   phased_snps.tsv   CHISEL's -l input ("#CHR POS PHASE", PHASE = 0|1 or 1|0),
#                     made from the phaser's VCF
#   rdr/ baf/ combo/ calls/ clones/ plots/   CHISEL's own folders
#   calls/calls.tsv   one row per bin x cell; CELL is the barcode
#                     (barcodedcells.info.tsv maps it back), CN_STATE = "A|B"
#   clones/mapping.tsv  cell -> cluster -> clone
#
# CHISEL's driver ignores the exit status of its sub-steps, and its error()
# exits 0, so success is decided here from the output files.
#
# All settings live in config.sh (CHISEL section).

LOG_TAG=call_chisel
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-ugp}

source ${CONDA}
conda activate ${CHISEL_ENV}

set -euo pipefail

# Defaults (override in config.sh)
CHISEL_CHROMS="${CHISEL_CHROMS:-${CHR}}"
CHISEL_BINSIZE="${CHISEL_BINSIZE:-5Mb}"
BIN_BP=$(to_bp "${CHISEL_BINSIZE}")          # also names the run folder
CHISEL_BLOCKSIZE="${CHISEL_BLOCKSIZE:-50kb}"
CHISEL_MINREADS="${CHISEL_MINREADS:-10000}"
CHISEL_MAXPLOIDY="${CHISEL_MAXPLOIDY:-3}"
CHISEL_UPPERK="${CHISEL_UPPERK:-100}"

# ============================================================
# Phaser -> phased VCF (and which sample column holds the phased genotype)
# ============================================================

case ${PHASER} in
    ugp)     VCF="${UGP_PHASED_VCF}"; SAMPLE="" ;;              # one sample (the bulk)
    perfect) VCF="${PERFECT_PHASED_VCF}"; SAMPLE="phasedgt" ;;  # true phase (phase_perfect.sh)
    *)   die "CHISEL with phaser '${PHASER}' is not implemented (${ROUTE_TODO[${PHASER}:chisel]:-})" ;;
esac

RUN_DIR=$(call_dir chisel "${PHASER}")
SNPS="${RUN_DIR}/phased_snps.tsv"
REF_DICT="${REF%.*}.dict"     # CHISEL reads chromosome lengths from <ref minus extension>.dict

# ============================================================
# Checks
# ============================================================

for tool in chisel bcftools samtools python2.7; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${CHISEL_ENV}"
done

for f in "${VCF}" "${CHISEL_BARCODED_BAM}" "${CHISEL_BARCODED_BAM}.bai" "${CHISEL_BARCODES}" \
         "${BULK_BAM}" "${REF}" "${REF}.fai"; do
    [[ -e $f ]] || die "missing $f"
done
[[ -e ${REF_DICT} ]] || die "missing ${REF_DICT} (CHISEL needs it; make it with: samtools dict ${REF} -o ${REF_DICT})"
[[ -e ${BULK_BAM}.bai || -e ${BULK_BAM%.bam}.bai ]] || die "missing index for ${BULK_BAM}"

N_PREP=$(grep -vc '^#' "${CHISEL_BARCODES}" || true)
note "Phaser: ${PHASER} (${VCF})"
note "Cells in barcoded BAM: ${N_PREP}"

mkdir -p "${RUN_DIR}"

# ============================================================
# Phased SNP list for CHISEL
#   biallelic SNPs, phased het genotype, chromosome names as in the BAM,
#   only CHISEL_CHROMS, one record per position
# ============================================================

if [[ -z ${SAMPLE} ]]; then
    SAMPLE=$(bcftools query -l "${VCF}" | awk 'NR == 1')
fi
[[ -n ${SAMPLE} ]] || die "no sample column in ${VCF}"

{
    printf '#CHR\tPOS\tPHASE\n'
    bcftools view -s "${SAMPLE}" -v snps -m2 -M2 -Ou "${VCF}" \
        | bcftools query -i 'GT="0|1" || GT="1|0"' -f '%CHROM\t%POS\t[%GT]\n' \
        | awk -F'\t' -v OFS='\t' -v want="${CHISEL_CHROMS}" '
            BEGIN { n = split(want, w, " ")
                    for (i = 1; i <= n; i++) { keep[w[i]] = 1; if (w[i] ~ /^chr/) usechr = 1 } }
            {
                c = $1
                if (usechr && c !~ /^chr/) c = "chr" c
                else if (!usechr) sub(/^chr/, "", c)
                if (!(c in keep)) next
                k = c ":" $2
                if (k in seen) next
                seen[k] = 1
                print c, $2, $3
            }'
} > "${SNPS}.tmp"
mv "${SNPS}.tmp" "${SNPS}"

N_SNPS=$(grep -vc '^#' "${SNPS}" || true)
(( N_SNPS > 0 )) || die "no phased het SNPs on '${CHISEL_CHROMS}' in ${VCF} (sample ${SAMPLE})"
note "Phased het SNPs for CHISEL: ${N_SNPS} -> ${SNPS}"

# ============================================================
# Run CHISEL
# ============================================================

# Clear earlier results so a failed step cannot leave an old calls.tsv behind
for d in rdr baf combo calls clones plots; do
    rm -rf "${RUN_DIR:?}/${d}"
done

note "Running CHISEL: bins ${BIN_BP} bp, blocks ${CHISEL_BLOCKSIZE}, minreads ${CHISEL_MINREADS}, maxploidy ${CHISEL_MAXPLOIDY}"
RC=0
chisel \
    -x "${RUN_DIR}" \
    -t "${CHISEL_BARCODED_BAM}" \
    -n "${BULK_BAM}" \
    -r "${REF}" \
    -l "${SNPS}" \
    -b "${BIN_BP}" \
    -k "${CHISEL_BLOCKSIZE}" \
    -c "${CHISEL_CHROMS}" \
    -m "${CHISEL_MINREADS}" \
    -p "${CHISEL_MAXPLOIDY}" \
    -K "${CHISEL_UPPERK}" \
    --seed "${CHISEL_SEED}" \
    -j "${CHISEL_JOBS}" \
    || RC=$?

# ============================================================
# Verify (step by step, so the error points at the step that failed)
# ============================================================

step_failed() {   # step-folder message
    local log="${RUN_DIR}/$1/log"
    echo "ERROR: CHISEL $1 step: $2" >&2
    if [[ -s ${log} ]]; then
        echo "--- last lines of ${log} ---" >&2
        tail -n 20 "${log}" >&2
    fi
    exit 1
}

# rdr/total.tsv: "normal" line + one line per selected cell (last line has no newline)
TOTAL="${RUN_DIR}/rdr/total.tsv"
[[ -s ${TOTAL} ]] || step_failed rdr "no total.tsv"
N_SEL=$(awk '$1 != "normal"' "${TOTAL}" | awk 'END { print NR }')
(( N_SEL > 0 )) || step_failed rdr "no cell has >= ${CHISEL_MINREADS} reads (CHISEL_MINREADS); lower it"
note "Cells passing minreads: ${N_SEL}/${N_PREP}"

[[ -s ${RUN_DIR}/baf/baf.tsv ]]     || step_failed baf "empty baf.tsv"
[[ -s ${RUN_DIR}/combo/combo.tsv ]] || step_failed combo "empty combo.tsv"

CALLS="${RUN_DIR}/calls/calls.tsv"
N_ROWS=$(awk 'END { print NR - 1 }' "${CALLS}" 2>/dev/null || echo 0)
(( N_ROWS > 0 )) || step_failed calls "no calls in ${CALLS}"
N_CALLED=$(awk -F'\t' 'NR > 1 { c[$4] = 1 } END { print length(c) }' "${CALLS}")
note "Calls: ${N_ROWS} bin x cell rows, ${N_CALLED} cells -> ${CALLS}"

# Cloning and plotting are useful but not needed for evaluation
[[ -s ${RUN_DIR}/clones/mapping.tsv ]] || note "WARNING: clones/mapping.tsv missing or empty (see ${RUN_DIR}/clones/log)"
ls "${RUN_DIR}"/plots/* >/dev/null 2>&1 || note "WARNING: no plots in ${RUN_DIR}/plots"
(( RC == 0 )) || note "WARNING: chisel exited ${RC}, but the calls are complete"

note "Done: ${RUN_DIR}"
