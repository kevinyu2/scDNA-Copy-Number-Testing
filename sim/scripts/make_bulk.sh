#!/bin/bash
#SBATCH --job-name=scDNA_sim_bulk
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/bulk_scDNA_sim_%j.out
#SBATCH --error=logs/bulk_scDNA_sim_%j.err

# Initialize conda
source "./config.sh"


source ${CONDA}
conda activate ${ENV}

set -euo pipefail



BULK_DIR="${PROJECT_DIR}/sim/bulk"
BULK_FASTQ_DIR="${BULK_DIR}/fastq"

mkdir -p "${BULK_FASTQ_DIR}"

# ============================================================
# Make combined maternal + paternal reference if needed
# ============================================================

if [[ ! -f "${SIM_FULL_FA}" ]]; then

    echo "Creating combined maternal/paternal reference:"
    echo "${SIM_FULL_FA}"

    cat \
        "${SIM_MAT_FA}" \
        "${SIM_PAT_FA}" \
        > "${SIM_FULL_FA}"

else

    echo "Combined reference already exists:"
    echo "${SIM_FULL_FA}"

fi

# ============================================================
# Run bulk DWGSIM
# ============================================================

echo "Running bulk DWGSIM..."

dwgsim -H -C 15 -1 150 -2 150 -e 0 -E 0 -r 0 "${SIM_FULL_FA}" "${BULK_FASTQ_DIR}/bulk"   # 15 per haplotype ≈ 30x total


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
    "${REF}" \
    "${BULK_READ1}" \
    "${BULK_READ2}" \
    | samtools sort \
        -o "${BULK_BAM}" \
        -

bwa mem -t 16 -R "@RG\tID:bulk\tSM:bulk\tLB:bulk\tPL:ILLUMINA" "${REF}" "${BULK_READ1}" "${BULK_READ2}" \
  | samtools sort -@ 4 -m 2G -o "${BULK_BAM}" -


echo "Indexing bulk BAM..."

samtools index "${BULK_BAM}"

echo
echo "Bulk BAM complete"
