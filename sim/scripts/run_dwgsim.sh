#!/bin/bash
#SBATCH --job-name=scDNA_sim
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/scDNA_sim_%A_%a.out
#SBATCH --error=logs/scDNA_sim_%A_%a.err

# Initialize conda

source "./config.sh"


source ${CONDA}
conda activate ${ENV}

set -euo pipefail


# Sample name comes from the SLURM array
SAMPLE_NAME="batch${SLURM_ARRAY_TASK_ID}"

python "./scDNA_sim.py" \
    --mat "${SIM_MAT_FA}" \
    --pat "${SIM_PAT_FA}" \
    --cnv "${TREE_FILE}" \
    --out "${PROJECT_DIR}/sim" \
    --sample-name "${SAMPLE_NAME}" \
    --ado-freq "${ADO_FREQ}" \
    --ado-mean "${ADO_MEAN}" \
    --ado-var "${ADO_VAR}" \
    --coverage-mean "${COVERAGE_MEAN}" \
    --coverage-var "${COVERAGE_VAR}"

mkdir -p "${PROJECT_DIR}/sim/bams"

# Now find this batch's FASTQs and make BAMs
for READ1 in "${PROJECT_DIR}/sim/dwgsim/${SAMPLE_NAME}"*_sim.bwa.read1.fastq.gz; do

    READ2="${READ1/.bwa.read1.fastq.gz/.bwa.read2.fastq.gz}"

    CELL_NAME="$(basename "${READ1}" .bwa.read1.fastq.gz)"
    BAM="${PROJECT_DIR}/sim/bams/${CELL_NAME}.bam"

    
    bwa mem \
        "${REF}" \
        "${READ1}" \
        "${READ2}" \
        | samtools sort \
            -o "${BAM}" \
            -

    samtools index "${BAM}"

done