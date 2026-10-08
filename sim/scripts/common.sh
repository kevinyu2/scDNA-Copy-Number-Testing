#!/bin/bash
# scripts/common.sh -- sourced by every script after config.sh.
#
# Holds the things that are derived, not chosen: where each stage reads and
# writes, which phaser -> caller routes exist, and two log helpers.
# Settings belong in config.sh, not here.
#
# Scripts set LOG_TAG before sourcing this (under SLURM, $0 is "slurm_script").

die()  { echo "ERROR: $*" >&2; exit 1; }
note() { echo "[${LOG_TAG:-pipeline}] $*" >&2; }

# Bin sizes pinned by run_pipeline.sh at submit time win over config.sh, so editing
# config.sh while jobs are queued cannot change the bin size of a run in flight.
[[ -n ${PIN_CHISEL_BINSIZE:-} ]] && CHISEL_BINSIZE=${PIN_CHISEL_BINSIZE}
[[ -n ${PIN_HS_BINSIZE:-} ]] && HS_BINSIZE=${PIN_HS_BINSIZE}

# "5Mb" / "500kb" / "5000000" -> bp. Same rules as CHISEL's -b: whole numbers,
# optional suffix kb or Mb (case-sensitive).
to_bp() {
    local v=$1
    case $v in
        *Mb) v=${v%Mb}; [[ $v =~ ^[0-9]+$ ]] && { echo $(( v * 1000000 )); return; } ;;
        *kb) v=${v%kb}; [[ $v =~ ^[0-9]+$ ]] && { echo $(( v * 1000 )); return; } ;;
        *)   [[ $v =~ ^[0-9]+$ ]] && { echo "$v"; return; } ;;
    esac
    die "bad bin size '$1' (use a whole number, optionally ending in kb or Mb, e.g. 5Mb, 500kb, 2000000)"
}

caller_binsize() {                                   # caller -> bin size in bp
    case $1 in
        hiscanner) to_bp "${HS_BINSIZE}" ;;
        chisel)    to_bp "${CHISEL_BINSIZE}" ;;
        *)         die "unknown caller '$1'" ;;
    esac
}

# ============================================================
# Routes
# ============================================================

PHASERS="scan2 ugp"
CALLERS="hiscanner chisel"
ROUTES="scan2:hiscanner ugp:chisel"          # phaser:caller pairs that are wired up

# Why a route is missing (shown to the user when they ask for it)
declare -A ROUTE_TODO=(
    [ugp:hiscanner]="HiScanner also reads SCAN2's raw GATK VCF (gatk/hc_raw.mmq60.vcf.gz), which UGP does not produce"
    [scan2:chisel]="call_chisel.sh would need to read SCAN2's eagle/phased_hets.vcf.gz (sample 'phasedgt'); not wired up yet"
)

route_ok() { [[ " ${ROUTES} " == *" $1:$2 "* ]]; }

# ============================================================
# Simulation outputs
# ============================================================

SIM_DIR="${PROJECT_DIR}/sim"
CELL_BAM_DIR="${SIM_DIR}/bams"
BULK_BAM="${SIM_DIR}/bulk/bams/bulk.bam"
CN_MAT_DIR="${SIM_DIR}/cn_mat"

# ============================================================
# Phaser outputs
# ============================================================

SCAN2_DIR="${PROJECT_DIR}/scan2_out"
# UGP_DIR is set in config.sh
UGP_PHASED_VCF="${UGP_DIR}/out/phase/phased_het_snps.vcf.gz"

# The file whose presence means "this phaser finished"
phaser_output() {
    case $1 in
        scan2) echo "${SCAN2_DIR}/shapeit/phased_hets.vcf.gz" ;;   # written by phase_scan2.sh's last steps
        ugp)   echo "${UGP_PHASED_VCF}" ;;
        *)     die "unknown phaser '$1'" ;;
    esac
}

# ============================================================
# Caller outputs:  ${PROJECT_DIR}/<caller>_<binsize bp>/<phaser>/
#   e.g. chisel_5000000/ugp/, hiscanner_500000/scan2/
# ============================================================

call_dir() { echo "${PROJECT_DIR}/$1_$(caller_binsize "$1")/$2"; }     # caller phaser

# CHISEL's barcoded BAM depends on neither the phaser nor the bin size, so it is shared
CHISEL_PREP_DIR="${PROJECT_DIR}/chisel/prep"
CHISEL_BARCODED_BAM="${CHISEL_PREP_DIR}/barcodedcells.bam"
CHISEL_BARCODES="${CHISEL_PREP_DIR}/barcodedcells.info.tsv"

caller_output() {                                    # caller phaser
    case $1 in
        hiscanner) echo "$(call_dir hiscanner "$2")/output/final_calls" ;;
        chisel)    echo "$(call_dir chisel "$2")/calls/calls.tsv" ;;
        *)         die "unknown caller '$1'" ;;
    esac
}

# ============================================================
# Evaluation:  ${EVAL_DIR}/<phaser>_<caller>_<binsize bp>/   (lifted truth shared)
#   e.g. eval/ugp_chisel_5000000/
# ============================================================

EVAL_TRUTH="${EVAL_DIR}/truth_hg38.tsv"
eval_out_dir() { echo "${EVAL_DIR}/$1_$2_$(caller_binsize "$2")"; }      # phaser caller
