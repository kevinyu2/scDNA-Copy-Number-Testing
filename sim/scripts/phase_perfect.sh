#!/bin/bash
#SBATCH --job-name=phase_perfect
#SBATCH --mem=32G
#SBATCH --cpus-per-task=16
#SBATCH --time=24:00:00
#SBATCH --output=logs/phase_perfect_%j.out
#SBATCH --error=logs/phase_perfect_%j.err

# "Perfect phasing": the simulated genome's TRUE phased het SNPs, so the callers
# can be tested without phasing error. Submitted by run_pipeline.sh --phaser perfect.
# One run serves both callers (HiScanner and CHISEL).
#
# 1. Align the HG002 mat/pat haplotypes to hg38 (minimap2; cached in EVAL_LIFT_DIR and
#    shared with eval.sh, so usually just reused).
# 2. perfect_phase.py: het SNPs = positions where exactly one haplotype differs from
#    hg38, both haplotypes align there exactly once, and no indel within
#    PERFECT_INDEL_PAD bp. GT = maternal|paternal.
#      -> ${PERFECT_DIR}/phased_hets.vcf.gz (sample "phasedgt")    CHISEL + HiScanner
#      -> ${PERFECT_DIR}/het_sites.tsv.gz
# 3. HiScanner also needs per-cell allele depths at those sites (it normally takes them
#    from SCAN2's GATK file). bcftools mpileup counts them in every cell BAM and the
#    bulk (chunks of PERFECT_CHUNK_CELLS BAMs in parallel), perfect_ad.py keeps the
#    true REF/ALT depths, and the chunks are merged:
#      -> ${PERFECT_DIR}/hiscanner_input/gatk/hc_raw.mmq60.vcf.gz   (written last = done)
#         ${PERFECT_DIR}/hiscanner_input/shapeit/phased_hets.vcf.gz -> ../../phased_hets.vcf.gz
#    i.e. the file names HiScanner expects from a SCAN2 folder.
#
# All settings live in config.sh (Perfect phasing section).

LOG_TAG=phase_perfect
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh

source ${CONDA}
conda activate ${PERFECT_ENV:-${ENV}}

set -euo pipefail

SCRIPTS="$(pwd)/scripts"
THREADS=${SLURM_CPUS_PER_TASK:-8}
PERFECT_MIN_MAPQ_ASM=${PERFECT_MIN_MAPQ_ASM:-5}
PERFECT_INDEL_PAD=${PERFECT_INDEL_PAD:-10}
PERFECT_MIN_MAPQ=${PERFECT_MIN_MAPQ:-60}
PERFECT_MIN_BASEQ=${PERFECT_MIN_BASEQ:-13}
PERFECT_CHUNK_CELLS=${PERFECT_CHUNK_CELLS:-50}

# ============================================================
# Checks
# ============================================================

for tool in minimap2 samtools bcftools bgzip tabix python3 md5sum; do
    command -v "$tool" >/dev/null || die "'$tool' not found in env ${PERFECT_ENV:-${ENV}}"
done
python3 -c "import numpy" 2>/dev/null || die "python needs numpy"
for f in "${SIM_MAT_FA}" "${SIM_PAT_FA}" "${REF}" "${REF}.fai" "${BULK_BAM}" \
         "${SCRIPTS}/perfect_phase.py" "${SCRIPTS}/perfect_ad.py"; do
    [[ -e $f ]] || die "missing $f"
done
mapfile -t BAMS < <(ls "${CELL_BAM_DIR}"/*.bam | LC_ALL=C sort)
(( ${#BAMS[@]} > 0 )) || die "no cell BAMs in ${CELL_BAM_DIR}"
note "${#BAMS[@]} cell BAMs + bulk"

CONTIG_LEN=$(awk -v c="${CHR}" '$1 == c { print $2 }' "${REF}.fai")
[[ -n ${CONTIG_LEN} ]] || die "${CHR} not in ${REF}.fai"

# Start clean: the completion marker goes first, so a failed run never looks finished
rm -f "${PERFECT_AD_VCF}" "${PERFECT_AD_VCF}.tbi" "${PERFECT_AD_VCF}.md5"
mkdir -p "${PERFECT_DIR}" "${PERFECT_HS_INPUT}/gatk" "${PERFECT_HS_INPUT}/shapeit"
WORK="${PERFECT_DIR}/work"
rm -rf "${WORK}"
mkdir -p "${WORK}"

# ============================================================
# 1. Haplotype -> hg38 alignments (cached, shared with eval.sh)
# ============================================================

make_lift_pafs

# ============================================================
# 2. True phased het SNPs
# ============================================================

note "Extracting true het SNPs (min alignment MAPQ ${PERFECT_MIN_MAPQ_ASM}, indel pad ${PERFECT_INDEL_PAD} bp)"
python3 -B "${SCRIPTS}/perfect_phase.py" \
    --mat-paf "$(lift_paf mat)" \
    --pat-paf "$(lift_paf pat)" \
    --target-chrom "${CHR}" \
    --query-chrom "${CHR}" \
    --min-mapq "${PERFECT_MIN_MAPQ_ASM}" \
    --indel-pad "${PERFECT_INDEL_PAD}" \
    --sample phasedgt \
    --out-vcf "${WORK}/phased_hets.vcf" \
    --out-sites "${WORK}/het_sites.tsv"

bgzip -f "${WORK}/phased_hets.vcf"
tabix -f -p vcf "${WORK}/phased_hets.vcf.gz"
bgzip -f "${WORK}/het_sites.tsv"
tabix -f -s1 -b2 -e2 "${WORK}/het_sites.tsv.gz"
mv "${WORK}/phased_hets.vcf.gz.tbi" "${PERFECT_PHASED_VCF}.tbi"
mv "${WORK}/phased_hets.vcf.gz"     "${PERFECT_PHASED_VCF}"
mv "${WORK}/het_sites.tsv.gz.tbi"   "${PERFECT_SITES}.tbi"
mv "${WORK}/het_sites.tsv.gz"       "${PERFECT_SITES}"
N_HETS=$(bcftools index -n "${PERFECT_PHASED_VCF}" 2>/dev/null || bcftools view -H "${PERFECT_PHASED_VCF}" | wc -l)
note "True het SNPs: ${N_HETS} -> ${PERFECT_PHASED_VCF}"

# HiScanner's file names (like phase_scan2.sh links eagle/ -> shapeit/)
ln -sfn ../../phased_hets.vcf.gz     "${PERFECT_HS_INPUT}/shapeit/phased_hets.vcf.gz"
ln -sfn ../../phased_hets.vcf.gz.tbi "${PERFECT_HS_INPUT}/shapeit/phased_hets.vcf.gz.tbi"
md5sum "${PERFECT_PHASED_VCF}" | awk '{ print $1 }' > "${PERFECT_HS_INPUT}/shapeit/phased_hets.vcf.gz.md5"

# ============================================================
# 3. Allele depths per cell at the true het sites (HiScanner's hc_raw.mmq60.vcf.gz)
# ============================================================

# Chunks of BAM paths; the bulk is its own chunk (sample names come from the @RG SM tags)
split -l "${PERFECT_CHUNK_CELLS}" -d -a 4 <(printf '%s\n' "${BAMS[@]}") "${WORK}/chunk_"
echo "${BULK_BAM}" > "${WORK}/chunk_bulk"
CHUNKS=("${WORK}"/chunk_[0-9]* "${WORK}/chunk_bulk")
note "Counting alleles: ${#CHUNKS[@]} chunks (${PERFECT_CHUNK_CELLS} BAMs each), ${THREADS} at a time;" \
     "MAPQ >= ${PERFECT_MIN_MAPQ}, base quality >= ${PERFECT_MIN_BASEQ}"

count_chunk() {   # chunk-list-file
    local list=$1 bcf="$1.bcf"
    bcftools mpileup -f "${REF}" -r "${CHR}" -T "${PERFECT_SITES}" -b "${list}" \
        -a FORMAT/AD -q "${PERFECT_MIN_MAPQ}" -Q "${PERFECT_MIN_BASEQ}" -B -I -d 100000 \
        -Ob -o "${bcf}" 2> "${list}.mpileup.log"
    bcftools query -l "${bcf}" > "${list}.samples"
    # every BAM must give its own sample column (duplicate SM tags would be merged)
    [[ $(wc -l < "${list}.samples") -eq $(wc -l < "${list}") ]] \
        || { echo "ERROR: ${list}: $(wc -l < "${list}.samples") samples from $(wc -l < "${list}") BAMs (missing or duplicate @RG SM?)" >&2; return 1; }
    bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%AD]\n' "${bcf}" \
        | python3 -B "${SCRIPTS}/perfect_ad.py" --sites "${PERFECT_SITES}" --samples "${list}.samples" \
                  --contig-length "${CONTIG_LEN}" 2> "${list}.ad.log" \
        | bgzip -c > "${list}.vcf.gz"
    tabix -f -p vcf "${list}.vcf.gz"
    rm -f "${bcf}"
}
export -f count_chunk
export REF CHR PERFECT_SITES PERFECT_MIN_MAPQ PERFECT_MIN_BASEQ SCRIPTS CONTIG_LEN

printf '%s\n' "${CHUNKS[@]}" | xargs -P "${THREADS}" -I{} bash -c 'set -euo pipefail; count_chunk "$1"' _ {} \
    || die "allele counting failed; see ${WORK}/chunk_*.log"
cat "${WORK}"/chunk_bulk.ad.log >&2

note "Merging ${#CHUNKS[@]} chunks"
printf '%s.vcf.gz\n' "${CHUNKS[@]}" > "${WORK}/merge.list"
bcftools merge -m none -l "${WORK}/merge.list" --threads "${THREADS}" -Oz -o "${WORK}/hc_raw.vcf.gz"
tabix -f -p vcf "${WORK}/hc_raw.vcf.gz"

N_SAMPLES=$(bcftools query -l "${WORK}/hc_raw.vcf.gz" | wc -l)
(( N_SAMPLES == ${#BAMS[@]} + 1 )) || die "merged VCF has ${N_SAMPLES} samples, expected $(( ${#BAMS[@]} + 1 ))"
bcftools query -l "${WORK}/hc_raw.vcf.gz" | grep -qx bulk || die "no 'bulk' sample in the merged VCF (bulk BAM @RG SM?)"

# Sanity check: at true het sites the 30x bulk should be ~50% ALT
bcftools query -s bulk -f '[%AD]\n' "${WORK}/hc_raw.vcf.gz" | awk -F, '
    $1 != "." { r = $1; a = $2; if (r + a >= 10) { n++; f += a / (r + a); if (a == 0 || r == 0) mono++ } }
    END { if (n) printf "[phase_perfect] bulk at true het sites (depth >= 10): %d sites, mean ALT fraction %.3f, one allele only %.2f%%\n", n, f / n, 100 * mono / n
          else print "[phase_perfect] WARNING: no bulk sites with depth >= 10" }' >&2

md5sum "${WORK}/hc_raw.vcf.gz" | awk '{ print $1 }' > "${PERFECT_AD_VCF}.md5"
mv "${WORK}/hc_raw.vcf.gz.tbi" "${PERFECT_AD_VCF}.tbi"
mv "${WORK}/hc_raw.vcf.gz"     "${PERFECT_AD_VCF}"      # completion marker (phaser_output perfect)

rm -rf "${WORK}"
note "Done: ${PERFECT_PHASED_VCF} (CHISEL, HiScanner) and ${PERFECT_AD_VCF} (HiScanner)"
