#!/bin/bash
#SBATCH --job-name=scDNA_sim_bulk
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/bulk_scDNA_sim_%j.out
#SBATCH --error=logs/bulk_scDNA_sim_%j.err

# Initialize conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scdna_pipeline

set -euo pipefail

source ../config.sh



FULL_REF="${SIM_FA_FOLDER}/full.fa"

BULK_DIR="${PROJECT_DIR}/sim/bulk"
BULK_FASTQ_DIR="${BULK_DIR}/fastq"

mkdir -p "${BULK_FASTQ_DIR}"

# ============================================================
# Make combined maternal + paternal reference if needed
# ============================================================

if [[ ! -f "${FULL_REF}" ]]; then

    echo "Creating combined maternal/paternal reference:"
    echo "${FULL_REF}"

    cat \
        "${SIM_FA_FOLDER}/mat.fa" \
        "${SIM_FA_FOLDER}/pat.fa" \
        > "${FULL_REF}"

else

    echo "Combined reference already exists:"
    echo "${FULL_REF}"

fi

# ============================================================
# Run bulk DWGSIM
# ============================================================

echo "Running bulk DWGSIM..."

dwgsim \
    -H \
    -C 30 \
    -1 150 \
    -2 150 \
    -e 0 \
    -E 0 \
    -r 0 \
    "${FULL_REF}" \
    "${BULK_FASTQ_DIR}/bulk"

echo "Bulk DWGSIM finished."




# ============================================================
# Bulk FASTQ → BAM
# ============================================================

BULK_BAM_DIR="${BULK_DIR}/bams"

mkdir -p "${BULK_BAM_DIR}"

BULK_READ1="${BULK_FASTQ_DIR}/bulk.bwa.read1.fastq.gz"
BULK_READ2="${BULK_FASTQ_DIR}/bulk.bwa.read2.fastq.gz"
BULK_BAM="${BULK_BAM_DIR}/bulk.bam"

echo "Aligning bulk reads..."

bwa mem \
    "${FULL_REF}" \
    "${BULK_READ1}" \
    "${BULK_READ2}" \
    | samtools sort \
        -o "${BULK_BAM}" \
        -

echo "Indexing bulk BAM..."

samtools index "${BULK_BAM}"

echo
echo "Bulk BAM complete"
