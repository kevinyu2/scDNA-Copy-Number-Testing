#!/bin/bash
#SBATCH --job-name=call_hiscanner
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=48:00:00
#SBATCH --output=logs/call_hiscanner_%j.out
#SBATCH --error=logs/call_hiscanner_%j.err

# Submitted by run_pipeline.sh --caller hiscanner (only from --phaser scan2).
# Output: ${PROJECT_DIR}/hiscanner_<binsize bp>/${PHASER}/output/final_calls
#
# All settings live in config.sh (HiScanner section).
# This script writes HiScanner's metadata.txt, config.yaml and cluster.yaml
# from those settings, checks inputs, then runs the requested steps.
#
# Normalize-step recovery: if BICseq normalization fails for some cells,
# those cells are retried (HS_NORM_RETRIES times), then dropped from
# metadata.txt and logged to ${HS_DIR}/excluded_cells.tsv, and the run
# continues. Exclusions persist across reruns until HS_RESET_EXCLUSIONS=true.

LOG_TAG=call_hiscanner
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
PHASER=${PHASER:-scan2}
[[ ${PHASER} == scan2 ]] || die "HiScanner from phaser '${PHASER}' is not implemented (${ROUTE_TODO[${PHASER}:hiscanner]:-})"
HS_DIR=$(call_dir hiscanner "${PHASER}")

source ${CONDA}
conda activate ${HS_ENV}

set -euo pipefail

SCAN2_OUT="${SCAN2_DIR}"
BAM_DIR="${CELL_BAM_DIR}"
BULK_SAMPLE="bulk"
HS_OUT="${HS_DIR}/output"
META="${HS_DIR}/metadata.txt"
EXCLUDED_TSV="${HS_DIR}/excluded_cells.tsv"

# Defaults for the recovery settings (override in config.sh)
HS_EXCLUDE_CELLS="${HS_EXCLUDE_CELLS:-}"
HS_MIN_CHR1_READS="${HS_MIN_CHR1_READS:-0}"
HS_AUTO_EXCLUDE="${HS_AUTO_EXCLUDE:-true}"
HS_NORM_RETRIES="${HS_NORM_RETRIES:-1}"
HS_MAX_EXCLUDE_PCT="${HS_MAX_EXCLUDE_PCT:-5}"
HS_RESET_EXCLUSIONS="${HS_RESET_EXCLUSIONS:-false}"

mkdir -p "${HS_DIR}"


# ============================================================
# Checks
# ============================================================

for tool in hiscanner bcftools samtools snakemake Rscript; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${HS_ENV}"
done

for f in \
    "${SCAN2_OUT}/gatk/hc_raw.mmq60.vcf.gz" \
    "${SCAN2_OUT}/gatk/hc_raw.mmq60.vcf.gz.tbi" \
    "${SCAN2_OUT}/shapeit/phased_hets.vcf.gz" \
    "${SCAN2_OUT}/shapeit/phased_hets.vcf.gz.tbi" \
    "${BULK_BAM}"; do
    [[ -f $f ]] || die "missing $f"
done

for c in ${HS_CHROMS}; do
    # HiScanner's Snakefile hard-codes the .fasta extension
    [[ -f "${FASTA_SPLIT_DIR}/${c}.fasta" ]] || die "missing ${FASTA_SPLIT_DIR}/${c}.fasta"
    [[ -f "${MAPPABILITY_STEM}${c}.txt" ]]   || die "missing mappability file ${MAPPABILITY_STEM}${c}.txt (check MAPPABILITY_STEM)"
done

# bcftools reads only the header; avoids the zcat|grep -m1 SIGPIPE problem under pipefail
LAST_SAMPLE=$(bcftools query -l "${SCAN2_OUT}/shapeit/phased_hets.vcf.gz" | tail -n 1)
[[ ${LAST_SAMPLE} == phasedgt ]] \
    || die "phased_hets.vcf.gz sample column is '${LAST_SAMPLE}', not phasedgt"

if [[ ${HS_USE_CLUSTER} == true && -z ${SLURM_ACCOUNT} ]]; then
    die "HS_USE_CLUSTER=true but SLURM_ACCOUNT is empty (HiScanner always passes --account to sbatch)"
fi

# ============================================================
# Exclusions carried over from earlier runs
# ============================================================

if [[ ${HS_RESET_EXCLUSIONS} == true && -f ${EXCLUDED_TSV} ]]; then
    mv "${EXCLUDED_TSV}" "${EXCLUDED_TSV}.$(date +%Y%m%d_%H%M%S).bak"
    note "HS_RESET_EXCLUSIONS=true: previous exclusions archived"
fi
if [[ ! -f ${EXCLUDED_TSV} ]]; then
    printf 'cell\tstep\treads_used\treason\ttime\n' > "${EXCLUDED_TSV}"
fi

declare -A AUTO_EXCLUDED=()
while IFS=$'\t' read -r cell _rest; do
    [[ -n ${cell} && ${cell} != cell ]] && AUTO_EXCLUDED["${cell}"]=1
done < "${EXCLUDED_TSV}"

# ============================================================
# metadata.txt  (bamID must equal the sample name in the SCAN2 VCF)
# ============================================================

FIRST_CHROM="${HS_CHROMS%% *}"
{
    printf 'bamID\tbam\tsinglecell\n'
    printf '%s\t%s\tN\n' "${BULK_SAMPLE}" "${BULK_BAM}"
    for b in "${BAM_DIR}"/*.bam; do
        cell=$(basename "$b" .bam)
        if [[ " ${HS_EXCLUDE_CELLS} " == *" ${cell} "* ]]; then
            note "Excluding ${cell} (HS_EXCLUDE_CELLS)"; continue
        fi
        if [[ -n ${AUTO_EXCLUDED[${cell}]+x} ]]; then
            note "Excluding ${cell} (listed in excluded_cells.tsv)"; continue
        fi
        if (( HS_MIN_CHR1_READS > 0 )); then
            n=$(samtools idxstats "$b" | awk -v c="${FIRST_CHROM}" '$1==c{print $3}')
            if (( ${n:-0} < HS_MIN_CHR1_READS )); then
                note "Excluding ${cell} (${n:-0} reads < HS_MIN_CHR1_READS=${HS_MIN_CHR1_READS})"; continue
            fi
        fi
        printf '%s\t%s\tY\n' "${cell}" "$b"
    done
} > "${META}"

N_CELLS=$(( $(wc -l < "${META}") - 2 ))
(( N_CELLS > 0 )) || die "no cells left in ${META}"
note "Wrote ${META}: ${N_CELLS} cells + bulk"

# Every bamID must be a sample in hc_raw.mmq60.vcf.gz, or HiScanner's bcftools -s fails
MISSING=$(comm -23 \
    <(tail -n +2 "${META}" | cut -f1 | LC_ALL=C sort) \
    <(bcftools query -l "${SCAN2_OUT}/gatk/hc_raw.mmq60.vcf.gz" | LC_ALL=C sort))
[[ -z ${MISSING} ]] || die "bamIDs not found in SCAN2 VCF (read-group SM mismatch?):
${MISSING}"

# ============================================================
# config.yaml
# ============================================================

CHROM_YAML=$(printf '"%s", ' ${HS_CHROMS})
CHROM_YAML="[${CHROM_YAML%, }]"

cat > "${HS_DIR}/config.yaml" <<EOF
# Generated by scripts/call_hiscanner.sh from config.sh -- edit config.sh, not this file
scan2_output: ${SCAN2_OUT}
metadata_path: ${META}
outdir: ${HS_OUT}
use_multisample_segmentation: ${HS_MULTISAMPLE}

fasta_folder: ${FASTA_SPLIT_DIR}
mappability_folder_stem: ${MAPPABILITY_STEM}

rdr_only: ${HS_RDR_ONLY}
binsize: $(to_bp "${HS_BINSIZE}")
chrom_list: ${CHROM_YAML}
lambda_range: ${HS_LAMBDA_RANGE}
lambda_value: ${HS_LAMBDA}
max_wgd: ${HS_MAX_WGD}
batch_size: ${HS_BATCH_SIZE}
depth_filter: ${HS_DEPTH_FILTER}
ado_threshold: ${HS_ADO_THRESHOLD}
threads: ${HS_THREADS}
ado_plot_baf_distribution: ${HS_ADO_PLOTS}

aggregate_every_k_snp: ${HS_AGGREGATE_K}
k: ${HS_K}

keep_raw_files: false
rerun: ${HS_RERUN}
keep_temp_files: false

use_cluster: ${HS_USE_CLUSTER}
cluster_queue: "${SLURM_PARTITION}"
max_cluster_jobs: 100
EOF

# ============================================================
# cluster.yaml  (only used by: hiscanner run --step normalize --use-cluster)
# ============================================================

cat > "${HS_DIR}/cluster.yaml" <<EOF
__default__:
  account: "${SLURM_ACCOUNT}"
  partition: "${SLURM_PARTITION}"
  time: "12:00:00"
  nodes: 1
  ntasks: 1
  mem: "4G"
  job-name: "hiscanner_{rule}"
  output: "cluster_logs/{rule}-%j.out"
  error: "cluster_logs/{rule}-%j.err"

get_read_pos:
  time: "4:00:00"
  mem: "4G"
  ntasks: 1
  partition: "${SLURM_PARTITION}"

run_bicseq_norm:
  time: "12:00:00"
  mem: "16G"
  ntasks: 4
  partition: "${SLURM_PARTITION}"
EOF

# ============================================================
# Helpers for normalize recovery
# ============================================================

# Cells in metadata.txt that are missing a non-empty bin file for any chromosome.
# Snakemake deletes the outputs of failed jobs, so these are the cells that failed.
cells_missing_bins() {
    tail -n +2 "${META}" | awk -F'\t' '$3=="Y"{print $1}' | while read -r cell; do
        for c in ${HS_CHROMS}; do
            if [[ ! -s "${HS_OUT}/bins/${cell}/${c}.bin" ]]; then
                echo "${cell}"; break
            fi
        done
    done
}

# Last error line from the cell's BICseq log, tabs/newlines flattened.
failure_reason() {
    local log="${HS_OUT}/logs/bicseq_norm_${1}.log" r=""
    if [[ -f ${log} ]]; then
        r=$(grep -E 'Error|error|Execution halted|Killed|CANCELLED|TIME LIMIT' "${log}" | grep -v '^Warning' | head -n 1 || true)
        [[ -z ${r} ]] && r=$(tail -n 1 "${log}" || true)
    else
        r="no BICseq log (job may not have started; check ${HS_OUT}/.workflow/cluster_logs)"
    fi
    echo "${r}" | tr '\t\n' '  '
}

reads_used() {
    local total=0 f
    for c in ${HS_CHROMS}; do
        f="${HS_OUT}/readpos/${1}/${c}.readpos.seq"
        [[ -f $f ]] && total=$(( total + $(wc -l < "$f") ))
    done
    echo "${total}"
}

exclude_cell() {   # cell step
    local cell=$1 step=$2
    printf '%s\t%s\t%s\t%s\t%s\n' "${cell}" "${step}" "$(reads_used "${cell}")" \
        "$(failure_reason "${cell}")" "$(date '+%F %T')" >> "${EXCLUDED_TSV}"
    awk -F'\t' -v c="${cell}" '$1!=c' "${META}" > "${META}.tmp" && mv "${META}.tmp" "${META}"
    note "Excluded ${cell} -> ${EXCLUDED_TSV}"
}

run_hiscanner() {  # step [extra flags...]
    local step=$1; shift
    note "=== hiscanner: ${step} $* ==="
    hiscanner --config config.yaml run --step "${step}" "$@"
}

run_normalize_with_recovery() {
    local flags=("$@") retries=0 rounds=0 n_failed n_meta
    local -a failed

    while :; do
        rounds=$(( rounds + 1 ))
        (( rounds <= HS_NORM_RETRIES + 3 )) || die "normalize still failing after ${rounds} rounds; stopping"

        if run_hiscanner normalize "${flags[@]+"${flags[@]}"}"; then
            return 0
        fi

        mapfile -t failed < <(cells_missing_bins)
        n_failed=${#failed[@]}
        (( n_failed > 0 )) || die "normalize failed but every cell has bin files; see ${HS_OUT}/.workflow logs"

        note "normalize failed for ${n_failed} cell(s): ${failed[*]}"

        # Never pass --rerun on recovery rounds: only the failed cells should be redone
        flags=()
        [[ ${HS_USE_CLUSTER} == true ]] && flags+=(--use-cluster)

        if (( retries < HS_NORM_RETRIES )); then
            retries=$(( retries + 1 ))
            note "Retry ${retries}/${HS_NORM_RETRIES} for the failed cells (BICseq sampling is random, a retry can pass)"
            continue
        fi

        [[ ${HS_AUTO_EXCLUDE} == true ]] \
            || die "normalize failed for: ${failed[*]} (HS_AUTO_EXCLUDE=false; add them to HS_EXCLUDE_CELLS to skip)"

        n_meta=$(tail -n +2 "${META}" | awk -F'\t' '$3=="Y"' | wc -l)
        if (( n_failed * 100 > HS_MAX_EXCLUDE_PCT * n_meta )); then
            die "${n_failed}/${n_meta} cells failed (> HS_MAX_EXCLUDE_PCT=${HS_MAX_EXCLUDE_PCT}%). Looks systemic, not low-coverage cells; not excluding. Example log: ${HS_OUT}/logs/bicseq_norm_${failed[0]}.log"
        fi

        for cell in "${failed[@]}"; do
            exclude_cell "${cell}" normalize
        done
        note "Re-running normalize without the excluded cells"
        retries=0
    done
}

# ============================================================
# Run
# ============================================================

# HiScanner looks for config.yaml / cluster.yaml in the working directory
cd "${HS_DIR}"

hiscanner validate config.yaml

FLAGS=()
[[ ${HS_RERUN} == true ]] && FLAGS+=(--rerun)

for step in ${HS_STEPS}; do
    if [[ ${HS_RDR_ONLY} == true && ${step} =~ ^(snp|phase|ado)$ ]]; then
        note "Skipping ${step} (RDR-only mode)"
        continue
    fi

    STEP_FLAGS=("${FLAGS[@]+"${FLAGS[@]}"}")

    if [[ ${step} == normalize ]]; then
        [[ ${HS_USE_CLUSTER} == true ]] && STEP_FLAGS+=(--use-cluster)
        run_normalize_with_recovery "${STEP_FLAGS[@]+"${STEP_FLAGS[@]}"}"
    else
        run_hiscanner "${step}" "${STEP_FLAGS[@]+"${STEP_FLAGS[@]}"}"
    fi
done

N_EXCL=$(( $(wc -l < "${EXCLUDED_TSV}") - 1 ))
note "HiScanner finished. Results: ${HS_OUT}/final_calls"
if (( N_EXCL > 0 )); then
    note "${N_EXCL} cell(s) excluded, see ${EXCLUDED_TSV}:"
    if command -v column >/dev/null; then
        column -t -s $'\t' "${EXCLUDED_TSV}" >&2
    else
        cat "${EXCLUDED_TSV}" >&2
    fi
fi