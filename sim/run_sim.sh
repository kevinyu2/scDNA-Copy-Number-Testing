#!/bin/bash
#SBATCH --job-name=scDNA_sim
#SBATCH --array=1-20
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=logs/scDNA_%A_%a.out
#SBATCH --error=logs/scDNA_%A_%a.err

# Initialize conda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate scDNA_sim

# Run simulation
python ./scDNA_sim.py \
    --mat /orcd/data/mgjones/001/kyu06/scDNA_test/hg002/HG002_mat_chr1.fa \
    --pat /orcd/data/mgjones/001/kyu06/scDNA_test/hg002/HG002_pat_chr1.fa \
    --cnv ./data/trees/mat_COLO.tree \
    --out /orcd/data/mgjones/001/kyu06/scDNA_test/COLO_mat_test \
    --sample-name "batch${SLURM_ARRAY_TASK_ID}"
