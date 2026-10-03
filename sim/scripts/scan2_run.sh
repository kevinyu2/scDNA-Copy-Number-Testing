#!/bin/bash
#SBATCH --job-name=scDNA_sim_scan2
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=logs/scan2_scDNA_sim_%j.out
#SBATCH --error=logs/scan2_scDNA_sim_%j.err

# Initialize conda
source "./config.sh"


source ${CONDA}
conda activate ${ENV}

set -euo pipefail

################################################
# Fix RG mistake

THREADS=${THREADS:-4}

add_rg() {
  bam=$1; thr=$2
  s=$(basename "$bam" .bam)                       # sample name = filename stem
  [[ $s == bulk ]] && s=bulk
  if [[ $(samtools view -H "$bam" | grep -c '^@RG') -gt 0 ]]; then
    echo "skip (has RG): $bam"; return
  fi
  tmp=$bam.rg.tmp.bam
  samtools addreplacerg -@ "$thr" \
    -r "@RG\tID:$s\tSM:$s\tLB:$s\tPL:ILLUMINA" \
    -o "$tmp" "$bam"
  samtools quickcheck "$tmp"
  mv "$tmp" "$bam"
  samtools index -@ "$thr" "$bam"
  echo "done: $bam"
}
export -f add_rg

{ find "$PROJECT_DIR/sim/bams" -maxdepth 1 -name '*.bam'
  echo "$PROJECT_DIR/sim/bulk/bams/bulk.bam"; } \
  | xargs -P 8 -I{} bash -c 'add_rg {} 2'

###################################################


scan2 -d "${PROJECT_DIR}/scan2_out" init



OUT_DIR=$PROJECT_DIR/scan2_out
scan2 -d "$OUT_DIR" init
cd "$OUT_DIR"

SC_ARGS=()
for b in "$BAM_DIR"/*.bam; do SC_ARGS+=(--sc-bam "$b"); done

scan2 config \
  --verbose \
  --sex male \
  --analysis call_mutations \
  --gatk gatk3_joint \
  --ref "$REF" \
  --genome hg38 \
  --dbsnp "$DBSNP" \
  --phaser eagle \
  --eagle-refpanel "$EAGLE_PANEL_DIR" \
  --eagle-genmap "$EAGLE_GENMAP" \
  --regions-file "$WINDOWS_BED" \
  --bulk-bam "$BULK_BAM" \
  "${SC_ARGS[@]}"

scan2 validate

mkdir -p logs
scan2 run \
  --joblimit 5000 \
  --cluster "sbatch --partition=mit_normal --cpus-per-task={threads} --mem={resources.mem_mb}M --time=72:00:00 --output=logs/%j.out" \
  --snakemake-args ' --until eagle_scatter --latency-wait 120'

