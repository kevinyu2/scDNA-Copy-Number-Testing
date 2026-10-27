#!/bin/bash

# config.sh -- every setting lives here. Scripts source it when a job STARTS,
# so edits affect queued and dependent jobs too.
# Derived paths (where each stage reads/writes) are in scripts/common.sh:
#   ${PROJECT_DIR}/sim/                      simulation
#   ${PROJECT_DIR}/scan2_out/, ${UGP_DIR}/   phasers
#   ${PROJECT_DIR}/<caller>_<binsize bp>/<phaser>/   callers, e.g. chisel_5000000/ugp/
#                                            (chisel/prep/ is shared by all)
#   ${EVAL_DIR}/<phaser>_<caller>_<binsize bp>/      evaluation, e.g. eval/ugp_chisel_5000000/

PROJECT_DIR="/n/fs/ragr-data/users/ky8418/scDNA_sim/projects/no_ecDNA_small"
CONDA=/n/fs/ragr-research/users/ky8418/miniconda3/etc/profile.d/conda.sh
ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/scdna_pipeline"

# Simulator
NUM_BATCHES=20
SIM_MAT_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_mat_chr1.fa"
SIM_PAT_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_pat_chr1.fa"
# This one gets created automatically if not present
SIM_FULL_FA="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/hg002/HG002_full_chr1.fa"
TREE_FILE="/n/fs/ragr-research/users/ky8418/scDNA_cn/scDNA-Copy-Number-Testing/sim/data/trees/small_no_ecDNA.tree"

# Simulation parameters
ADO_FREQ=0.000005
ADO_MEAN=3
ADO_VAR=1
COVERAGE_MEAN=0.1
COVERAGE_VAR=0.3
SIM_READ_LEN=150          # bp per read (paired), cells and bulk. Must be >= the HiScanner
                          # mappability k-mer (HS_MAPPABILITY_K; the hg38 track is 150mer)


# References
REF="/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/hg38/Homo_sapiens_assembly38.fasta"


# Resources for scan2 and hiscanner
RESOURCES_DIR=/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/
WINDOWS_BED="/n/fs/ragr-data/users/ky8418/scDNA_sim/sim_resources/chr1_2000windows.bed"
GENOME="hg38"
CHR="chr1"

# Delete GATK per-chunk files once the merged VCF exists
SCAN2_CLEAN_CHUNKS=true

# ============================================================
# SLURM 
# ============================================================
SLURM_PARTITION="cs"
SLURM_ACCOUNT="allcs"


# ============================================================
# HiScanner
# ============================================================
HS_ENV="${ENV}"                                   # change if HiScanner lives in its own env
# Run dir is ${PROJECT_DIR}/hiscanner_<HS_BINSIZE>/<phaser>/ (config.yaml, metadata.txt, output/)

# Reference files
FASTA_SPLIT_DIR="${RESOURCES_DIR%/}/hg38_split"   # must contain <chrom>.fasta (e.g. chr1.fasta)
MAPPABILITY_STEM="${RESOURCES_DIR%/}/hg38_mappability/150mer."
HS_MAPPABILITY_K=150                 # k-mer length of that track; call_hiscanner.sh refuses shorter reads

# Chromosomes to analyze (space-separated, names as in the BAM/VCF)
HS_CHROMS="${CHR}"

# Steps to run, in order (drop ones you've finished; snp/phase/ado are skipped if HS_RDR_ONLY=true)
HS_STEPS="snp phase ado normalize segment cnv"

# Analysis parameters (defaults from HiScanner's default_config.yaml)
HS_RDR_ONLY=false
HS_MULTISAMPLE=false                # use_multisample_segmentation
HS_BINSIZE=500000
HS_LAMBDA=4                         # BIC-seq penalty used for the calls (must be in HS_LAMBDA_RANGE).
                                    # HiScanner's 16 is tuned for deeper data: at 0.1x with 500 kb bins it
                                    # misses 5 Mb single-copy gains even in clean simulated cells; 4 finds
                                    # 2-5 Mb gains/losses. All lambdas in the range are segmented anyway, so
                                    # trying another value only needs HS_STEPS="cnv" + eval
HS_SEG_NOSCALE=false                # false (HiScanner default): BIC-seq scales lambda per cell by a noise factor;
                                    # high-copy ecDNA inflates it ~6x+, which can hide smaller events elsewhere
                                    # in those cells. true: --noscale. Changing this needs HS_STEPS="segment cnv"
HS_LAMBDA_RANGE="[2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]"
HS_MAX_WGD=1 # Maximum whole genome duplication it tests. Set at 1 for no WGD
HS_BATCH_SIZE=5
HS_DEPTH_FILTER=0
HS_ADO_THRESHOLD=0.2
HS_ADO_PLOTS=true
HS_AGGREGATE_K=true                 # aggregate_every_k_snp
HS_K=5
HS_THREADS=16

# BIC-seq normalization (patched into a copy of HiScanner's Snakefile by call_hiscanner.sh)
HS_NORM_P=0.01                       # fraction of positions sampled to fit the GC/mappability model.
                                     # HiScanner's 0.0001 leaves ~15 reads in the fit at 0.1x on chr1
                                     # (unstable "expected", random normalize failures). 0.01 -> ~1.5k reads,
                                     # ~2M rows in the R fit (fits in run_bicseq_norm's 16G)
HS_NORM_READLEN=150                  # BIC-seq -l (its default is 50); = SIM_READ_LEN
HS_NORM_FRAGSIZE=500                 # BIC-seq -s (default 300); dwgsim's default outer distance is 500

# Run control
HS_USE_CLUSTER=true              
HS_RERUN=true                      # true = redo steps even if outputs exist



# ============================================================
# Universal Genotyping Pipeline
# ============================================================

UGP_ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/genotyping-env"
UGP_REPO="/n/fs/ragr-research/users/ky8418/Universal-Genotyping-Pipeline"
UGP_CONDA_PREFIX="/n/fs/ragr-data/users/ky8418/scDNA_sim/ugp_conda"   # same as in step 1
UGP_DIR="${PROJECT_DIR}/ugp_pseudobulk"
UGP_SAMPLE_ID="HG002sim"

UGP_MODE="pseudobulk"      # bulk | pseudobulk | percell
#   bulk:       bulk normal only; stops after phasing (all CHISEL with fixed/--varbins bins needs)
#   pseudobulk: bulk normal + ONE scDNA dataset, chisel_prep's barcodedcells.bam (all cells pooled;
#               UGP piles it up with --ignore-RG); makes bins. chisel_prep runs before UGP.
#   percell:    bulk normal + every cell as its own dataset. UGP keeps only SNPs covered in EVERY
#               dataset (phase_and_concat_bulk, min_depth floored at 1), so at ~0.1-0.3x per cell
#               all SNPs are dropped and binning fails. Only usable with deeply sequenced cells.
UGP_MAX_CELLS=20           # percell only; 0 = all cells

UGP_CHROMS="1"             # UGP wants bare numbers; it handles the chr prefix itself
UGP_SNP_PANEL="${RESOURCES_DIR%/}/1kGP.chr1.snps.vcf.gz"
UGP_PHASING_PANEL="${RESOURCES_DIR%/}/ugp_phasing_panel"
UGP_GTF="${RESOURCES_DIR%/}/gencode.v38.annotation.gtf.gz"
EAGLE_GENMAP="${RESOURCES_DIR%/}/genetic_map_hg38_withX.txt.gz"

# Adaptive binning (pseudobulk/percell; bulk mode stops before binning)
# Thresholds are set as MEAN SNP-covering reads per cell per bin. In pseudobulk mode UGP sees one
# tumor dataset, the pool of all cells, so it (and bins_perfect.sh) gets the POOLED count
# value x PHASER_BINS_NCELLS, and its outputs keep pooled names: 40/cell x 1000 cells -> MSR40000/.
# 40/cell is roughly 500 kb mean bins on chr1 here. run_pipeline.sh --phaserbins N sets N per cell, and
# the phase stage makes exactly that threshold (UGP min_snp_reads: [N x PHASER_BINS_NCELLS]); without
# --phaserbins it makes CHISEL_PHASER_PER_CELL (below). Other thresholds already on disk are kept.
# Cells in the pool, from the simulation: NUM_BATCHES x the cells of one tree (the sum of its leaves,
# first line of TREE_FILE, e.g. ((20[..],10[..])[..],20[..]) = 50). 0 if TREE_FILE cannot be read.
tree_cells() { [[ -r $1 ]] || { echo 0; return; }; awk 'NF && $0 !~ /^#/ { print; exit }' "$1" | sed 's/\[[^]]*\]//g' | grep -o '[0-9]\+' | awk '{ s += $1 } END { print s + 0 }'; }
PHASER_BINS_NCELLS=$(( NUM_BATCHES * $(tree_cells "${TREE_FILE}") ))
UGP_MIN_SNP_PER_BIN=1                # het SNPs per bin (UGP default: 1)
UGP_MIN_TOTAL_READS=1000             # read starts per bin in every dataset, incl. the bulk; 0 disables (UGP default: 5000)
UGP_GENE_AWARE_BINNING=true          # never cut inside a gene (UGP default: true)


UGP_CORES=16
UGP_THREADS_GENOTYPE=1     # bcftools mpileup is single-threaded; more doesn't help
UGP_THREADS_PHASE=16       # Eagle uses all of these
UGP_THREADS_PILEUP=1       # per-cell jobs (percell mode): 1 each lets 16 cells run at once
UGP_THREADS_MOSDEPTH=1
UGP_EXECUTOR=local         # local = run inside this job (best for bulk mode); slurm = one cluster job per step

# ============================================================
# Perfect phasing (--phaser perfect): true phased hets from the HG002 assemblies
# ============================================================
# Output: ${PROJECT_DIR}/perfect/ (phased_hets.vcf.gz for both callers, hiscanner_input/ for HiScanner)
PERFECT_ENV="${ENV}"         # needs minimap2, samtools, bcftools, bgzip/tabix, python3 + numpy
PERFECT_MIN_MAPQ_ASM=5       # haplotype -> hg38 alignments (eval's cached PAFs) used for het calling
PERFECT_INDEL_PAD=10         # skip SNPs within this many bp of an indel in either haplotype
PERFECT_MIN_MAPQ=60          # read MAPQ for HiScanner's allele depths (SCAN2's file is "mmq60")
PERFECT_MIN_BASEQ=13         # base quality for allele depths (bcftools default)
PERFECT_CHUNK_CELLS=50       # BAMs per bcftools mpileup job (they run in parallel)

# ============================================================
# CHISEL
# ============================================================
# Barcoded BAM: ${PROJECT_DIR}/chisel/prep/ (shared). Runs: ${PROJECT_DIR}/chisel_<binsize bp>/<phaser>/
CHISEL_ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/chisel"           # original CHISEL
NEW_CHISEL_ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/new_chisel"   # chisel_hmm fork
# The fork is used automatically for --varbins, --phaserbins and bins_perfect.sh; everything else
# (chisel_prep, fixed bins) uses CHISEL_ENV. CHISEL_USE_NEW=true also runs fixed bins with the fork,
# into chisel_new<bp>/ folders so the original CHISEL's chisel_<bp>/ results are kept.
CHISEL_USE_NEW=false
# HMM bin clustering (the chisel_hmm fork's hmm branch, HATCHet2-style; run_pipeline.sh --hmm K or
# --hmm MIN-MAX). An HMM run only re-does the calling (Caller with the HMM, Cloner, Plotter) on the
# combo.tsv of the same run without --hmm, which must exist; folders chisel_<bins>_hmm<K>/<phaser>/.
HMM_CHISEL_ENV="/n/fs/ragr-research/users/ky8418/miniconda3/envs/chisel_hmmenv"
CHISEL_HMM=""                 # set per run with --hmm (K, or MIN-MAX to let the criterion choose)
CHISEL_HMM_SELECT="bic"       # bic | silhouette (only matters with a MIN-MAX range)
CHISEL_HMM_RESTARTS=10        # restarts per K
CHISEL_HMM_EXTRA=""           # any other chisel_calling --hmm* options, e.g. "--hmmtransmat fixed --hmmdecode viterbi"
CHISEL_ORACLE_SPLIT=0         # set per run with --oracle-split P: boundary bins (< P% in one oracle cluster) get a
                              # cluster each instead of joining their majority; folders chisel_<bins>_oracle_split<P>/
CHISEL_CLONING=hamming        # set per run with --cloning divisive: re-do only the cloning of the run with the fork's
                              # Cloner --clustering divisive (HMM_CHISEL_ENV); folders chisel_<bins>[_hmm<K>|_oracle]_divclone/
CHISEL_ORACLE=false           # set per run with --oracle: the calling with the oracle bin clusters (scripts/oracle_clusters.py,
                              # from eval/truth_hg38.tsv) instead of computed ones; folders chisel_<bins>_oracle/<phaser>/
CHISEL_SEED=12
CHISEL_JOBS=16
CHISEL_BARCODE_LENGTH=12
CHISEL_PREP_FORCE=false     # true = rebuild barcodedcells.bam even if it exists

# chisel run (needs ${REF%.fasta}.dict next to the reference)
CHISEL_CHROMS="${CHR}"      # space-separated, names as in the BAM
CHISEL_BINSIZE="5Mb"        # whole number + optional kb/Mb; names the run folder in bp (chisel_5000000/)
                            # run_pipeline.sh --binsize overrides it per run. ~3k reads per cell per 5 Mb bin at 0.1x
CHISEL_BLOCKSIZE="50kb"     # haplotype block size for BAF (CHISEL default; 0 disables)
CHISEL_MINREADS=10000       # cells with fewer reads (MAPQ>=13, on CHISEL_CHROMS) are dropped.
                            # CHISEL's default (300000) would drop every cell here: 0.1x on chr1
                            # is ~165k reads per cell
CHISEL_MAXPLOIDY=3          # base ploidies tried are 2, 4, ... up to this: 3 = diploid only
                            # (no WGD, like HS_MAX_WGD=1); 4 = also test a WGD
CHISEL_UPPERK=100           # max bin clusters
CHISEL_EVAL_CN="corrected"  # corrected = clone consensus (CORRECTED_HAP_CN, CHISEL's final output); raw = each cell's own call (HAP_CN)

# Variable-length bins (chisel --varbins, needs the chisel_hmm fork). CHISEL_BINSIZE is then only
# used for the first, fixed-bin RDR pass that selects cells. run_pipeline.sh --varbins N overrides
# CHISEL_VARBINS/CHISEL_VB_MINSNPREADSPERCELL per run; the run folder is chisel_var<N>/<phaser>/
CHISEL_VARBINS=false              # true = variable-length bins
CHISEL_VB_MINSNPREADSPERCELL=25   # mean SNP-covering reads per cell per bin (x #cells = pseudobulk threshold)
CHISEL_VB_MINTOTALREADS=5000      # read starts per bin in the cell pseudobulk and in the bulk normal
CHISEL_VB_MINSNPS=1               # het SNPs per bin
CHISEL_VB_WINDOW="1kb"            # windows that are merged into bins
CHISEL_VB_REFVERSION="auto"       # chromosome arms bound the bins: auto (from chr1 length) | hg19 | hg38 | none

# Bins made by the phaser (chisel --bins); run_pipeline.sh --phaserbins N (alias --ugpbins N) sets
# these per run; the run folder is chisel_pbins<N>/<phaser>/. N = SNP reads PER CELL per bin; the
# bins files use the pooled count P = N x PHASER_BINS_NCELLS. --phaserbins N makes them if missing:
#   ugp:     UGP's bins, ${UGP_DIR}/out/bb/MSR<P>/bulk/bb.tsv.gz. Needs UGP_MODE=pseudobulk.
#   perfect: the same window merge on the TRUE het SNPs (scripts/bins_perfect.sh, CHISEL fork's
#            chisel_bins in pooled mode): pooled SNP reads of all cells at the true sites (summed
#            from PERFECT_AD_VCF), read starts of chisel_prep's pooled BAM and of the bulk, UGP's
#            UGP_MIN_SNP_PER_BIN, UGP_MIN_TOTAL_READS and blacklist, chromosome arms. Not UGP's
#            gene-aware cuts. -> ${PERFECT_DIR}/bb/MSR<P>/bins.tsv
CHISEL_PHASER_BINS=false
CHISEL_PHASER_PER_CELL=40        # N when --phaserbins is not given (e.g. --phaser ugp --to phase)
PERFECT_BINS_BLACKLIST="${UGP_REPO}/resources/data/hg38-blacklist.v2.bed.gz"   # "" = none


# ============================================================
# Evaluation
# ============================================================
EVAL_ENV="${ENV}"
EVAL_DIR="${PROJECT_DIR}/eval"
EVAL_LIFT_DIR="${RESOURCES_DIR%/}/liftover"   # haplotype->hg38 alignments, shared across projects
EVAL_PURITY=0.9           # min fraction of a bin in one true segment to score it
EVAL_ORIENTATION=global   # global | cell
EVAL_CELLS_PER_PAGE=4
EVAL_MAX_CELLS_PLOT=0     # 0 = plot every cell
EVAL_SEG_MIN_FRAC=0.75    # threshold for confusion_segments_*_{bp,cells}_min75.pdf: a bin counts toward a segment
                          # only if >= this fraction of it (its aligned truth) lies inside. The unthresholded
                          # _bp / _cells versions are always made too.

EVAL_MIN_COVERED=0.5      # min fraction of a bin with aligned truth to score it
EVAL_Y_LINEAR_MAX=8       # track plots: CN above this goes above an axis break
EVAL_MERGE_GAP=1000       # bridge indels up to this size (bp) when lifting truth
EVAL_MAX_MEAN_DEV=0.5     # skip bins where small high-CN truth pieces (ecDNA) dominate read depth (when majority is this much different from the mean)
                          # NOTE: if using small sections, may want to set this to inf to turn it off