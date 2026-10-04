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




OUT_DIR=$PROJECT_DIR/scan2_out
scan2 -d "$OUT_DIR" init
cd "$OUT_DIR"

RES=${RESOURCES_DIR%/}
BAM_DIR=$PROJECT_DIR/sim/bams
BULK_BAM=$PROJECT_DIR/sim/bulk/bams/bulk.bam
DBSNP=$RES/common_all_20180418.chrprefix.vcf
EAGLE_GENMAP=$RES/genetic_map_hg38_withX.txt.gz
EAGLE_PANEL_DIR=$RES/eagle_1000g_panel

SC_ARGS=()
for b in "$BAM_DIR"/*.bam; do SC_ARGS+=(--sc-bam "$b"); done

scan2 config \
  --verbose \
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

# This install's Snakefiles read keys its `scan2 config` doesn't write
add_key() {
  grep -q "^$1:" scan.yaml || { echo "$1: $2" >> scan.yaml; echo "added to scan.yaml -> $1: $2"; }
}
add_key phased_hsnps "''"          # needed: read when the workflow loads
add_key phased_hsnps_n_cores 20    # needed: read when the workflow loads
add_key fdr 0.01                   # old name for target_fdr
add_key min_bulk_dp 11             # old name for snv_min_bulk_dp
add_key min_sc_alt 2               # old name for snv_min_sc_alt
add_key min_sc_dp 6                # old name for snv_min_sc_dp

scan2 validate

mkdir -p logs
scan2 run \
  --joblimit 5000 \
  --cluster "sbatch --partition=mit_normal --cpus-per-task={threads} --mem={resources.mem_mb}M --time=72:00:00 --output=logs/%j.out" \
  --snakemake-args ' --until phasing_gather --latency-wait 120'



# ---- Prepare SCAN2 output for HiScanner ----
cd "$OUT_DIR"

# Compressed copy of the raw calls; the original .vcf stays for Snakemake
if [[ ! -f gatk/hc_raw.mmq60.vcf.gz ]]; then
  bgzip -c gatk/hc_raw.mmq60.vcf > gatk/hc_raw.mmq60.vcf.gz
fi
tabix -f -p vcf gatk/hc_raw.mmq60.vcf.gz

# HiScanner expects shapeit/; point it at the Eagle results
if [[ -d shapeit && ! -L shapeit ]]; then
  # a real shapeit/ dir exists: link the files individually
  ln -sf ../eagle/phased_hets.vcf.gz     shapeit/phased_hets.vcf.gz
  ln -sf ../eagle/phased_hets.vcf.gz.tbi shapeit/phased_hets.vcf.gz.tbi
else
  ln -sfn eagle shapeit
fi

# silences HiScanner's missing-md5 warnings
md5sum gatk/hc_raw.mmq60.vcf.gz    > gatk/hc_raw.mmq60.vcf.gz.md5
md5sum eagle/phased_hets.vcf.gz    > eagle/phased_hets.vcf.gz.md5


# ---- Remove GATK per-chunk files once the merged VCF is in place ----
if [[ ${SCAN2_CLEAN_CHUNKS:-true} == true ]]; then
  if [[ -s gatk/hc_raw.mmq60.vcf.gz && -s gatk/hc_raw.mmq60.vcf.gz.tbi ]]; then
    find gatk -maxdepth 1 \( -name 'hc_raw.mmq60_chunk*.vcf' \
                          -o -name 'hc_raw.mmq60_chunk*.vcf.idx' \
                          -o -name 'scatter_benchmark.mmq60_chunk*.tsv' \) -delete
    echo "Removed GATK chunk files"
  else
    echo "Merged VCF missing; keeping chunk files" >&2
  fi
fi
