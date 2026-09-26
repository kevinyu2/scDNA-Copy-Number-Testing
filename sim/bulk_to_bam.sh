#!/bin/bash
#SBATCH --job-name=bulk_bam
#SBATCH --output=logs/bulk_bam_%j.out
#SBATCH --error=logs/bulk_bam_%j.err
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=02:00:00



eval "$(conda shell.bash hook)"
conda activate scDNA_sim

set -euo pipefail

FASTQ_DIR="/orcd/data/mgjones/001/kyu06/scDNA_test/bulk/fastq/"
BAM_DIR="/orcd/data/mgjones/001/kyu06/scDNA_test/bulk/bams/"
REF="/orcd/data/mgjones/001/kyu06/scDNA_test/hg38/Homo_sapiens_assembly38.fasta"





R1="$FASTQ_DIR/bulk.bwa.read1.fastq.gz"
R2="$FASTQ_DIR/bulk.bwa.read2.fastq.gz"

mkdir -p "$BAM_DIR"

echo "Aligning bulk..."

bwa mem -t "$SLURM_CPUS_PER_TASK" \
    -R "@RG\tID:bulk\tSM:bulk\tPL:ILLUMINA" \
    "$REF" \
    "$R1" \
    "$R2" |
    samtools sort \
        -@ "$SLURM_CPUS_PER_TASK" \
        -o "$BAM_DIR/bulk.bam" -

samtools index "$BAM_DIR/bulk.bam"

echo "Finished bulk."
