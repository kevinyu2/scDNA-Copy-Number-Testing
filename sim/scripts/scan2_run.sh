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



scan2 -d "${PROJECT_DIR}/scan2_out" init



# ---- paths ----
RES=${RESOURCES_DIR%/}
BAM_DIR=$PROJECT_DIR/sim/bams
BULK_BAM=$PROJECT_DIR/sim/bulk/bams/bulk.bam
BULK_SAMPLE=bulk
DBSNP=$RES/common_all_20180418.chrprefix.vcf
EAGLE_GENMAP=$RES/genetic_map_hg38_withX.txt.gz
EAGLE_PANEL_DIR=$RES/eagle_1000g_panel
SCAN2_LIB=$ENV/lib/scan2

OUT_DIR=$PROJECT_DIR/scan2
OUT=$OUT_DIR/scan.yaml

# ---- checks ----
for f in "$BULK_BAM" "$WINDOWS_BED" "$DBSNP" "$EAGLE_GENMAP" "$REF"; do
  [[ -f $f ]] || { echo "Missing file: $f" >&2; exit 1; }
done
[[ -d $BAM_DIR ]] || { echo "Missing dir: $BAM_DIR" >&2; exit 1; }

  [[ -f $DBSNP.idx ]] || { echo "Missing Tribble index: $DBSNP.idx (make with igvtools)" >&2; exit 1; }

# ---- metadata ----
UUID=$(uuidgen 2>/dev/null || python3 -c 'import uuid; print(uuid.uuid4())')
CREATE_DATE=$(date '+%Y-%m-%d %H:%M %Z')
CREATE_TIME=$(python3 -c 'import time; print(time.time())')
CREATOR=${USER:-$(whoami)}

# ---- single-cell BAMs: every *.bam in the folder, sample = basename minus .bam ----
mapfile -t SAMPLES < <(find "$BAM_DIR" -maxdepth 1 -name '*.bam' -printf '%f\n' \
                       | sed 's/\.bam$//' | LC_ALL=C sort)
[[ ${#SAMPLES[@]} -gt 0 ]] || { echo "No BAMs in $BAM_DIR" >&2; exit 1; }

sc_block() {
  for s in "${SAMPLES[@]}"; do printf '  %s: %s/%s.bam\n' "$s" "$BAM_DIR" "$s"; done
}
# bam_map = single cells + bulk (LC_ALL=C sorted, so "bulk" lands after batch9_*)
bam_map_block() {
  { sc_block; printf '  %s: %s\n' "$BULK_SAMPLE" "$BULK_BAM"; } | LC_ALL=C sort -t: -k1,1
}
eagle_block() {
  for c in $(printf '%s\n' {1..22} X | LC_ALL=C sort); do
    printf '  chr%s: %s/chr%s_GRCh38.genotypes.bcf\n' "$c" "$EAGLE_PANEL_DIR" "$c"
  done
}
# BED 0-based half-open -> 1-based inclusive
regions_block() {
  awk 'NF>=3 {printf "- %s:%d-%d\n", $1, $2+1, $3}' "$WINDOWS_BED"
}

# ---- write (temp file, then atomic replace) ----
mkdir -p "$OUT_DIR"
TMP=$(mktemp "$OUT_DIR/.scan2.yaml.XXXXXX")
trap 'rm -f "$TMP"' EXIT

cat > "$TMP" <<EOF
abmodel_chunk_strategy: false
abmodel_chunks: 10
abmodel_hsnp_n_tiles: 250
abmodel_hsnp_tile_size: 100
abmodel_n_cores: 20
abmodel_refine_steps: 4
abmodel_samples_per_step: 20000
abmodel_use_fit: {}
add_muts: null
analysis: call_mutations
analysis_uuid: $UUID
analyze_indels: true
analyze_mosaic_snvs: false
analyze_snvs: true
bam_map:
$(bam_map_block)
bams: {}
bulk_bam: $BULK_BAM
bulk_sample: $BULK_SAMPLE
callable_regions: true
chr_prefix: chr
chrs:
- $CHR
compute_sensitivity: false
create_date: $CREATE_DATE
create_time: $CREATE_TIME
creator: $CREATOR
cross_sample_panel: null
dbsnp: $DBSNP
digest_depth_n_cores: 20
eagle_genmap: $EAGLE_GENMAP
eagle_refpanel:
$(eagle_block)
gatk: gatk3_joint
gatk_chunks: 2000
gatk_regions:
$(regions_block)
gatk_vcf: null
genome: $GENOME
genotype_n_cores: 20
indel_max_bulk_af: 0
indel_max_bulk_alt: 0
indel_min_bulk_dp: 11
indel_min_sc_alt: 2
indel_min_sc_dp: 10
integrate_table_n_cores: 20
is_started: false
makepanel_metadata: null
makepanel_n_cores: 20
mimic_legacy: false
min_base_quality_score: 20
parsimony_phasing: false
permtool_bedtools_genome_file: null
permtool_callable_bed_n_cores: 20
permtool_combine_permutations_n_cores: 20
permtool_config_map: {}
permtool_indel_generation_param: 0.02
permtool_make_permutations_n_cores: 20
permtool_matrix_map: {}
permtool_muts: null
permtool_n_permutations: 10000
permtool_snv_generation_param: 100000
phased_hsnps: ''
phased_hsnps_n_cores: 20
phaser: eagle
ref: $REF
resample_M: 20
rescue_n_cores: 12
rescue_target_fdr: 0.01
sc_bams:
$(sc_block)
scan2_objects: {}
scripts: $SCAN2_LIB
sensitivity_n_cores: 20
shapeit_refpanel: null
snakefile: $SCAN2_LIB/Snakefile
snv_max_bulk_af: 0
snv_max_bulk_alt: 0
snv_min_bulk_dp: 11
snv_min_sc_alt: 2
snv_min_sc_dp: 6
target_fdr: 0.01
EOF

mv "$TMP" "$OUT"
trap - EXIT
echo "Wrote $OUT: ${#SAMPLES[@]} cells, $(grep -c '^- chr[0-9X]*:' "$OUT") regions"


scan2 run \
  --joblimit 5000 \
  --cluster "sbatch --partition=mit_normal --cpus-per-task={threads} --mem={resources.mem_mb}M --time=72:00:00 --output=%logdir/%j.out" \
  --snakemake-args ' --until eagle_scatter --latency-wait 120'