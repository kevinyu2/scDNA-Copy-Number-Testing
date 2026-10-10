#!/bin/bash
#SBATCH --job-name=sim_cells
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=08:00:00
#SBATCH --output=logs/sim_cells_%A_%a.out
#SBATCH --error=logs/sim_cells_%A_%a.err

# One SLURM array task per batch: simulate the batch's cells (scDNA_sim.py ->
# dwgsim), then align each cell with its own read group. Submitted by run_sim.sh.

LOG_TAG=sim_cells
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh


source ${CONDA}
conda activate ${ENV}

set -euo pipefail


# Sample name comes from the SLURM array
SAMPLE_NAME="batch${SLURM_ARRAY_TASK_ID}"

python "./scripts/scDNA_sim.py" \
    --mat "${SIM_MAT_FA}" \
    --pat "${SIM_PAT_FA}" \
    --cnv "${TREE_FILE}" \
    --out "${SIM_DIR}" \
    --sample-name "${SAMPLE_NAME}" \
    --ado-freq "${ADO_FREQ}" \
    --ado-mean "${ADO_MEAN}" \
    --ado-var "${ADO_VAR}" \
    --coverage-mean "${COVERAGE_MEAN}" \
    --coverage-var "${COVERAGE_VAR}" \
    --read-length "${SIM_READ_LEN}"

mkdir -p "${CELL_BAM_DIR}"

# Now find this batch's FASTQs and make BAMs
for READ1 in "${SIM_DIR}/dwgsim/${SAMPLE_NAME}_"*_sim.bwa.read1.fastq.gz; do
    READ2="${READ1/.bwa.read1.fastq.gz/.bwa.read2.fastq.gz}"

    CELL_NAME="$(basename "${READ1}" .bwa.read1.fastq.gz)"
    BAM="${CELL_BAM_DIR}/${CELL_NAME}.bam"

    
    bwa mem \
        -t 16 \
        -R "@RG\tID:${CELL_NAME}\tSM:${CELL_NAME}\tLB:${CELL_NAME}\tPL:ILLUMINA" \
        "${REF}" \
        "${READ1}" \
        "${READ2}" \
        | samtools sort \
            -o "${BAM}" \
            -

    samtools index "${BAM}"

done