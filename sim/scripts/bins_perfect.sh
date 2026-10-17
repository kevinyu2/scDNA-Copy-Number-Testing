#!/bin/bash
#SBATCH --job-name=bins_perfect
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --time=8:00:00
#SBATCH --output=logs/bins_perfect_%j.out
#SBATCH --error=logs/bins_perfect_%j.err
# Variable-length bins from the TRUE het SNPs, made the way UGP's pseudobulk mode makes its bins,
# so CHISEL can be run on them like on UGP's (run_pipeline.sh --phaser perfect --phaserbins N).
# Submitted by run_pipeline.sh after phase_perfect.sh and chisel_prep.sh.
#
# 1. Pooled SNP-covering reads per true het SNP: REF+ALT depths of every cell (not the bulk),
#    summed from phase_perfect's ${PERFECT_AD_VCF}       -> ${PERFECT_DIR}/bb/pooled_snp_depth.tsv.gz
# 2. For every N in PERFECT_BINS_MSR: the CHISEL fork's chisel_bins in pooled mode, i.e. UGP's
#    window merge (1 kb windows within chromosome arms minus PERFECT_BINS_BLACKLIST) until the
#    pooled SNP reads >= N, het SNPs >= UGP_MIN_SNP_PER_BIN and read starts >= UGP_MIN_TOTAL_READS
#    in both chisel_prep's pooled BAM and the bulk. Window counts are cached, so only the first
#    threshold reads the BAMs.                           -> ${PERFECT_DIR}/bb/MSR<N>/bins.tsv
#    (+ bins.summary.tsv: #bins, mean and median size; bins.pdf, windows.tsv.gz, log_binning).
#    Existing outputs are kept.
# Differences to UGP's bins: the true het SNPs instead of UGP's, no gene-aware boundaries.
LOG_TAG=bins_perfect
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
source ${CONDA}
conda activate "${NEW_CHISEL_ENV:?NEW_CHISEL_ENV is not set in config.sh}"   # chisel_bins comes with the fork
set -euo pipefail
THREADS=${SLURM_CPUS_PER_TASK:-8}
CHISEL_CHROMS="${CHISEL_CHROMS:-${CHR}}"

for f in "${PERFECT_AD_VCF}" "${BULK_BAM}" "${CHISEL_BARCODED_BAM}" "${REF}" "${REF%.fasta}.dict"; do
    [[ -s $f ]] || die "missing input: $f"
done
for tool in bcftools samtools chisel_bins chisel_binsplot; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${NEW_CHISEL_ENV} (chisel_bins comes with the chisel_hmm fork: pip install -e .)"
done
chisel_bins -h | grep -q -- "--snps" || die "chisel_bins has no pooled mode (--snps); update the chisel_hmm fork"

BB_DIR="${PERFECT_DIR}/bb"
mkdir -p "${BB_DIR}"
SNPS="${BB_DIR}/pooled_snp_depth.tsv.gz"
NCELLS=$(bcftools query -l "${PERFECT_AD_VCF}" | grep -vcx bulk)

if [[ ! -s ${SNPS} ]]; then
    note "Summing allele depths of ${NCELLS} cells at the true het SNPs"
    bcftools query -s ^bulk -f '%CHROM\t%POS[\t%AD]\n' "${PERFECT_AD_VCF}" \
        | awk -F'\t' -v OFS='\t' '{ d = 0; for (i = 3; i <= NF; i++) { n = split($i, a, ","); for (j = 1; j <= n; j++) if (a[j] != ".") d += a[j] } print $1, $2, d }' \
        | gzip > "${SNPS}.tmp"
    mv "${SNPS}.tmp" "${SNPS}"
fi
note "Pooled SNP depths: $(zcat "${SNPS}" | wc -l) SNPs, $(zcat "${SNPS}" | awk '{ s += $3 } END { print s }') reads (${SNPS})"

BL_ARGS=()
if [[ -n ${PERFECT_BINS_BLACKLIST:-} ]]; then
    [[ -s ${PERFECT_BINS_BLACKLIST} ]] || die "missing PERFECT_BINS_BLACKLIST ${PERFECT_BINS_BLACKLIST}"
    BL_ARGS=(--blacklist "${PERFECT_BINS_BLACKLIST}")
fi

for MSR in $(echo "${PERFECT_BINS_MSR}" | tr -d '[],'); do
    OUT="${BB_DIR}/MSR${MSR}"
    if [[ -s ${OUT}/bins.tsv ]]; then
        note "MSR${MSR}: ${OUT}/bins.tsv exists, kept"
        continue
    fi
    mkdir -p "${OUT}"
    note "MSR${MSR}: binning (pooled SNP reads >= ${MSR}, i.e. $(awk -v m="${MSR}" -v n="${NCELLS}" 'BEGIN { printf "%.1f", m / n }') per cell)"
    chisel_bins -t "${CHISEL_BARCODED_BAM}" -n "${BULK_BAM}" -r "${REF}" -c "${CHISEL_CHROMS}" \
        --snps "${SNPS}" --ncells "${NCELLS}" \
        --minsnpreads "${MSR}" --minsnpreadspercell 0 \
        --mintotalreads "${UGP_MIN_TOTAL_READS}" --mintotalreadspercell 0 --minsnps "${UGP_MIN_SNP_PER_BIN}" \
        "${BL_ARGS[@]+"${BL_ARGS[@]}"}" \
        --counts "${BB_DIR}/window_counts.npz" --windows "${OUT}/windows.tsv.gz" --summary "${OUT}/bins.summary.tsv" \
        -j "${THREADS}" --outdir "${OUT}" > "${OUT}/bins.tsv.tmp" 2> "${OUT}/log_binning" \
        || { tail -n 20 "${OUT}/log_binning" >&2; die "chisel_bins failed for MSR${MSR}"; }
    mv "${OUT}/bins.tsv.tmp" "${OUT}/bins.tsv"
    note "MSR${MSR}: $(grep -vc '^#' "${OUT}/bins.tsv") bins -> ${OUT}/bins.tsv"
    chisel_binsplot "${OUT}/bins.tsv" "${OUT}/windows.tsv.gz" -o "${OUT}/bins.pdf" 2>> "${OUT}/log_binning" \
        || note "WARNING: plotting MSR${MSR} failed (see ${OUT}/log_binning)"
done
note "Done: ${BB_DIR}/MSR*/bins.tsv"
