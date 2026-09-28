#!/bin/bash
#SBATCH --job-name=sc_bam
#SBATCH --output=logs/sc_bam_%A_%a.out
#SBATCH --error=logs/sc_bam_%A_%a.err
#SBATCH --array=1-20
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00

# Initialize conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scDNA_sim


set -euo pipefail


FASTQ_DIR="/orcd/data/mgjones/001/kyu06/scDNA_test/COLO_mat_test/dwgsim"
BAM_DIR="/orcd/data/mgjones/001/kyu06/scDNA_test/COLO_mat_test/bams"
REF="/orcd/data/mgjones/001/kyu06/scDNA_test/hg38/Homo_sapiens_assembly38.fasta"

mkdir -p "$BAM_DIR" logs

# One batch per array task
BATCH="batch${SLURM_ARRAY_TASK_ID}"


echo "Processing ${BATCH}"

for R1 in "$FASTQ_DIR"/"${BATCH}"_cell*_sim.bwa.read1.fastq.gz; do

    BASE=$(basename "$R1" .bwa.read1.fastq.gz)
    R2="$FASTQ_DIR/${BASE}.bwa.read2.fastq.gz"

    echo "Processing ${BASE}"

    if [[ ! -f "$R2" ]]; then
        echo "ERROR: missing R2: $R2"
        exit 1
    fi

    bwa mem -t "$SLURM_CPUS_PER_TASK" \
        -R "@RG\tID:${BASE}\tSM:${BASE}\tPL:ILLUMINA" \
        "$REF" \
        "$R1" \
        "$R2" |
        samtools sort \
            -@ "$SLURM_CPUS_PER_TASK" \
            -o "$BAM_DIR/${BASE}.bam" -

    samtools index "$BAM_DIR/${BASE}.bam"

done

echo "Finished ${BATCH}"