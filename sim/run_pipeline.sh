#!/bin/bash
# run_pipeline.sh -- phasing -> CN calling -> evaluation, as chained SLURM jobs.
# Run from sim/ on the login node (it only submits jobs). See usage() below.

LOG_TAG=run_pipeline
source "${CONFIG_FILE:-./config.sh}"
source ./scripts/common.sh
set -euo pipefail

STAGES=(phase call eval)

usage() {
    cat <<EOF
Usage: ./run_pipeline.sh --phaser P [--caller C] [options]

  --phaser P     ${PHASERS// / | }
  --caller C     ${CALLERS// / | }    (not needed with --to phase)
  --from S       first stage: phase | call | eval     (default: phase)
  --to S         last stage:  phase | call | eval     (default: eval)
  --only S       same as --from S --to S
  --binsize B    caller bin size for this run (e.g. 1Mb, 500kb, 2000000); overrides
                 CHISEL_BINSIZE / HS_BINSIZE from config.sh. Output folders are
                 named by it: <caller>_<bp>/ and eval/<phaser>_<caller>_<bp>/
  --varbins N    CHISEL only: variable-length bins with N mean SNP-covering reads per
                 cell per bin (chisel --varbins); output folders are chisel_var<N>/
  --phaserbins N CHISEL only: variable-length bins made by the phaser, N mean SNP reads per
                 cell per bin: --phaser ugp has UGP make them (UGP_MODE=pseudobulk; the phase
                 stage builds threshold N x cells), --phaser perfect makes them from the true
                 het SNPs (scripts/bins_perfect.sh). Folders chisel_pbins<N>/. With --to phase
                 (no caller) it only sets which UGP bins to build
  --ugpbins N    same as --phaserbins N
  --hmm K        CHISEL only: re-do the calling with HMM bin clustering (chisel_hmm's hmm branch,
                 HMM_CHISEL_ENV) on the combo.tsv of the same run without --hmm, with exactly K
                 clusters, or MIN-MAX to let CHISEL_HMM_SELECT choose; needs --from call.
                 Folders chisel_<bins>_hmm<K>/ (e.g. --phaserbins 8 --hmm 7: chisel_pbins8_hmm7/)
  --oracle       CHISEL only, for diagnosis: re-do the calling of the run without --oracle with the
                 oracle bin clusters (same copy numbers in the same cells = one cluster, from
                 eval/truth_hg38.tsv; scripts/oracle_clusters.py) instead of computed ones; needs
                 --from call. Folders chisel_<bins>_oracle/
  --oracle-split P  as --oracle, but bins with less than P% of their bp in one oracle cluster get a
                 cluster each (boundary bins kept out of the clusters they dilute). Folders
                 chisel_<bins>_oracle_split<P>/
  --mode M       CHISEL --hmm / --oracle runs and --reclone: the fork's mode (default CHISEL_MODE,
                 config.sh): pseudobulk = divisive cell groups before calling, HMM on their pooled
                 reads, clone calls from pooled reads (the fork's default; folders ..._pb/);
                 divisive = divisive cloning without pseudobulk (..._div/); legacy = CHISEL's own
                 cloning (no suffix: the folders made before these modes)
  --reclone      CHISEL only: re-do only the cloning of an existing legacy run (also --hmm or
                 --oracle) in the --mode, on a copy of its calls; needs --from call (or --from eval
                 to re-score). Folders ..._pbclone/ (pseudobulk) or ..._divclone/ (divisive)
                 (--hmm and --oracle re-do only the calling on the plain CHISEL run's combo.tsv; from the
                 phase stage, or when that run does not exist yet, it is queued first, with --oracle also
                 its eval for the lifted truth)
  --hmmnoise N   with --hmm: the HMM's noise model (default CHISEL_HMM_NOISE, config.sh): counting =
                 known read-count variances of every bin + overdispersion and extra variance shared
                 by all states (the fork's default; folders ..._hmm<K>_cn...); legacy = every state its
                 own variances (no suffix: the HMM folders made before)
  --cloning divisive  the same as --reclone --mode divisive (the older option)
  --after IDS    wait for these SLURM job IDs (colon- or comma-separated) before
                 starting; run_sim.sh passes its job IDs here
  --dry-run      print the sbatch commands without submitting
  --list         show implemented routes and exit

Prints the last submitted job ID on stdout (everything else goes to stderr), so
  P=$(./run_pipeline.sh --phaser perfect --to phase)
  ./run_pipeline.sh --phaser perfect --caller chisel --from call --after "$P"

Stages
  phase   scan2: SCAN2 (GATK + Eagle)       ugp: Universal Genotyping Pipeline (bulk mode)
          perfect: true phase from the HG002 assemblies (+ per-cell allele depths for
                   HiScanner); one run serves both callers
  call    hiscanner: HiScanner              chisel: chisel_prep (barcoded BAM) + CHISEL
  eval    lift truth to hg38, convert calls to a common table, score them

Examples
  ./run_pipeline.sh --phaser ugp --caller chisel                # all three stages
  ./run_pipeline.sh --phaser ugp --caller chisel --from call    # UGP already done
  ./run_pipeline.sh --phaser scan2 --caller hiscanner --only eval
  ./run_pipeline.sh --phaser ugp --to phase                     # just phasing
  ./run_pipeline.sh --phaser ugp --caller chisel --from call --binsize 1Mb   # a second bin size
  ./run_sim.sh -- --phaser ugp --caller chisel                  # simulate, then all of this

Outputs: ${PROJECT_DIR}/calls/<caller>_<binsize bp>/<phaser>/ and ${EVAL_DIR}/<phaser>_<caller>_<binsize bp>/
EOF
}

list_routes() {
    echo "Implemented phaser -> caller routes:"
    for r in ${ROUTES}; do echo "  ${r%%:*} -> ${r##*:}"; done
    echo "Not implemented yet:"
    for p in ${PHASERS}; do for c in ${CALLERS}; do
        route_ok "$p" "$c" || echo "  $p -> $c: ${ROUTE_TODO[$p:$c]:-no reason recorded}"
    done; done
}

stage_index() {
    local i
    for i in "${!STAGES[@]}"; do [[ ${STAGES[$i]} == "$1" ]] && { echo "$i"; return; }; done
    die "unknown stage '$1' (use: ${STAGES[*]})"
}

# ============================================================
# Arguments
# ============================================================

PHASER="" CALLER="" FROM=phase TO=eval AFTER="" DRY=false BINSIZE="" VARBINS="" UGPBINS=""
while (( $# )); do
    case $1 in
        --phaser)  PHASER=${2:?--phaser needs a value}; shift 2 ;;
        --caller)  CALLER=${2:?--caller needs a value}; shift 2 ;;
        --from)    FROM=${2:?--from needs a value};     shift 2 ;;
        --to)      TO=${2:?--to needs a value};         shift 2 ;;
        --only)    FROM=${2:?--only needs a value}; TO=$2; shift 2 ;;
        --after)   AFTER=${2:?--after needs a value};   shift 2 ;;
        --binsize) BINSIZE=${2:?--binsize needs a value}; shift 2 ;;
        --varbins) VARBINS=${2:?--varbins needs a value}; shift 2 ;;
        --ugpbins|--phaserbins) UGPBINS=${2:?$1 needs a value}; shift 2 ;;
        --hmm)     HMMK=${2:?--hmm needs a value}; shift 2 ;;
        --oracle)  ORACLE=true; shift ;;
        --oracle-split) ORACLE=true; ORACLE_SPLIT=${2:?--oracle-split needs a percentage}; shift 2 ;;
        --cloning) CLONING=${2:?--cloning needs a method}; shift 2 ;;
        --mode)    MODE=${2:?--mode needs pseudobulk, divisive or legacy}; shift 2 ;;
        --hmmnoise) HMMNOISE=${2:?--hmmnoise needs counting or legacy}; shift 2 ;;
        --reclone) RECLONE=true; shift ;;
        --dry-run) DRY=true; shift ;;
        --list)    list_routes; exit 0 ;;
        -h|--help) usage; exit 0 ;;
        *)         usage >&2; die "unknown option '$1'" ;;
    esac
done
AFTER=${AFTER//,/:}

FROM_I=$(stage_index "${FROM}")
TO_I=$(stage_index "${TO}")
(( FROM_I <= TO_I )) || die "--from ${FROM} comes after --to ${TO}"
runs() { local i; i=$(stage_index "$1"); (( FROM_I <= i && i <= TO_I )); }

[[ -n ${PHASER} ]] || { usage >&2; die "--phaser is required"; }
[[ " ${PHASERS} " == *" ${PHASER} "* ]] || die "unknown phaser '${PHASER}' (use: ${PHASERS})"

if (( TO_I >= 1 )); then
    [[ -n ${CALLER} ]] || die "--caller is required when running past the phase stage"
    [[ " ${CALLERS} " == *" ${CALLER} "* ]] || die "unknown caller '${CALLER}' (use: ${CALLERS})"
    if ! route_ok "${PHASER}" "${CALLER}"; then
        echo >&2
        echo "The ${PHASER} -> ${CALLER} route is not implemented yet:" >&2
        echo "  ${ROUTE_TODO[${PHASER}:${CALLER}]:-no reason recorded}" >&2
        echo >&2
        list_routes >&2
        exit 1
    fi
elif [[ -n ${CALLER} ]]; then
    note "--to phase: ignoring --caller ${CALLER}"
    CALLER=""
fi

# Bin size: --binsize overrides config.sh; either way it is pinned into the jobs
if [[ -n ${BINSIZE} ]]; then
    [[ -n ${CALLER} ]] || die "--binsize needs --caller"
    case ${CALLER} in
        chisel)    CHISEL_BINSIZE=${BINSIZE} ;;
        hiscanner) HS_BINSIZE=${BINSIZE} ;;
    esac
fi
# (plain assignments, so a bad size stops the launcher before anything is submitted)
if [[ -n ${VARBINS} ]]; then
    [[ ${CALLER} == chisel ]] || die "--varbins needs --caller chisel"
    [[ ${VARBINS} =~ ^[0-9]+([.][0-9]+)?$ ]] || die "--varbins needs a number of SNP reads per cell, e.g. 25"
    CHISEL_VARBINS=true
    CHISEL_VB_MINSNPREADSPERCELL=${VARBINS}
fi
if [[ -n ${UGPBINS} ]]; then
    [[ ${CALLER} == chisel || -z ${CALLER} ]] || die "--phaserbins needs --caller chisel (or --to phase)"
    [[ -z ${VARBINS} ]] || die "use either --varbins or --phaserbins"
    [[ ${UGPBINS} =~ ^[0-9]+$ ]] && (( UGPBINS > 0 )) || die "--phaserbins needs a whole number of SNP reads per cell, e.g. 40"
    CHISEL_PHASER_BINS=true
    CHISEL_PHASER_PER_CELL=${UGPBINS}
fi
if [[ ${CHISEL_PHASER_BINS:-false} == true && ( ${CALLER} == chisel || -z ${CALLER} ) ]]; then
    case ${PHASER} in
        ugp)
            [[ ${UGP_MODE} == pseudobulk || ${UGP_MODE} == percell ]] \
                || die "UGP writes bins only with UGP_MODE=pseudobulk (config.sh has '${UGP_MODE}')" ;;
        perfect) ;;
        *)  die "phaser '${PHASER}' makes no bins; --phaserbins works with --phaser ugp or perfect" ;;
    esac
    (( PHASER_BINS_NCELLS > 0 )) || die "cannot count cells: TREE_FILE=${TREE_FILE} is not readable"
    NBAMS=$(compgen -G "${CELL_BAM_DIR}/*.bam" | wc -l || true)   # compgen fails when nothing matches (pipefail)
    (( NBAMS == 0 || NBAMS == PHASER_BINS_NCELLS )) \
        || note "WARNING: ${NBAMS} cell BAMs in ${CELL_BAM_DIR}, but NUM_BATCHES x tree cells = ${PHASER_BINS_NCELLS}"
    [[ ${CHISEL_VARBINS:-false} != true ]] || die "CHISEL_VARBINS and CHISEL_PHASER_BINS are both true"
fi
if runs phase && [[ ${PHASER} == ugp && ( ${UGP_MODE} == pseudobulk || ${UGP_MODE} == percell ) ]]; then
    (( PHASER_BINS_NCELLS > 0 )) || die "cannot count cells: TREE_FILE=${TREE_FILE} is not readable"
    note "UGP bins: ${CHISEL_PHASER_PER_CELL} SNP reads per cell -> MSR$(( CHISEL_PHASER_PER_CELL * PHASER_BINS_NCELLS ))$([[ -n ${UGPBINS} ]] || echo " (CHISEL_PHASER_PER_CELL; pass --phaserbins N for another)")"
fi
if [[ -n ${HMMK:-} ]]; then
    [[ ${CALLER} == chisel ]] || die "--hmm needs --caller chisel"
    [[ ${HMMK} =~ ^[0-9]+(-[0-9]+)?$ ]] || die "--hmm needs a number of clusters K or a range MIN-MAX, e.g. 7 or 5-10"
    if [[ ${HMMK} == *-* ]]; then
        (( ${HMMK%-*} >= 1 && ${HMMK%-*} <= ${HMMK#*-} )) || die "--hmm ${HMMK}: need 1 <= MIN <= MAX"
    else
        (( HMMK >= 1 )) || die "--hmm needs at least 1 cluster"
    fi
    CHISEL_HMM=${HMMK}
    note "HMM: ${HMMK} clusters, on the combo.tsv of $(chisel_base_dir "${PHASER}")"
fi
if [[ ${ORACLE:-false} == true ]]; then
    [[ ${CALLER} == chisel ]] || die "--oracle needs --caller chisel"
    [[ -z ${HMMK:-} ]] || die "use either --hmm or --oracle"
    CHISEL_ORACLE=true
    if [[ -n ${ORACLE_SPLIT:-} ]]; then
        [[ ${ORACLE_SPLIT} =~ ^[0-9]+$ ]] && (( ORACLE_SPLIT >= 1 && ORACLE_SPLIT <= 100 )) \
            || die "--oracle-split needs a whole percentage 1-100, e.g. 75"
        CHISEL_ORACLE_SPLIT=${ORACLE_SPLIT}
        note "Boundary bins (< ${ORACLE_SPLIT}% in one oracle cluster) get a cluster each"
    fi
    note "Oracle bin clusters (from ${EVAL_TRUTH}), on the combo.tsv of $(chisel_base_dir "${PHASER}")"
fi
if [[ -n ${CLONING:-} ]]; then                     # the older option: --cloning divisive = --reclone --mode divisive
    case ${CLONING} in
        divisive) RECLONE=true; MODE=divisive ;;
        hamming)  ;;
        *)        die "--cloning: unknown method '${CLONING}' (divisive; or use --reclone --mode M)" ;;
    esac
fi
if [[ -n ${MODE:-} ]]; then
    [[ ${MODE} =~ ^(pseudobulk|divisive|legacy)$ ]] || die "--mode needs pseudobulk, divisive or legacy (not '${MODE}')"
    CHISEL_MODE=${MODE}
fi
if [[ -n ${HMMNOISE:-} ]]; then
    [[ ${HMMNOISE} =~ ^(counting|legacy)$ ]] || die "--hmmnoise needs counting or legacy (not '${HMMNOISE}')"
    [[ -n ${HMMK:-} ]] || die "--hmmnoise goes with --hmm (it is the HMM's noise model)"
    CHISEL_HMM_NOISE=${HMMNOISE}
fi
if [[ ${RECLONE:-false} == true ]]; then
    [[ ${CALLER} == chisel ]] || die "--reclone needs --caller chisel"
    [[ ${CHISEL_MODE:-pseudobulk} != legacy ]] || die "--reclone needs --mode pseudobulk or divisive (the run itself is the legacy one)"
    (( FROM_I >= 1 )) || die "--reclone re-does only the cloning: start at the call stage (--from call)"
    CHISEL_RECLONE=true
    note "Re-cloning (${CHISEL_MODE:-pseudobulk}) the calls of $(clone_base_dir "${PHASER}")"
elif [[ ${CALLER} == chisel ]] && { [[ -n ${CHISEL_HMM:-} ]] || [[ ${CHISEL_ORACLE:-false} == true ]]; }; then
    note "Fork mode: ${CHISEL_MODE:-pseudobulk}$([[ -n ${CHISEL_HMM:-} ]] && echo ", HMM noise: ${CHISEL_HMM_NOISE:-counting}")"
fi
PIN_C=$(to_bp "${CHISEL_BINSIZE}")
PIN_H=$(to_bp "${HS_BINSIZE}")
PINS="PIN_CHISEL_BINSIZE=${PIN_C},PIN_HS_BINSIZE=${PIN_H},PIN_CHISEL_VARBINS=${CHISEL_VARBINS:-false},PIN_CHISEL_VB_MSR=${CHISEL_VB_MINSNPREADSPERCELL:-25},PIN_CHISEL_PHASER_BINS=${CHISEL_PHASER_BINS:-false},PIN_CHISEL_PHASER_PER_CELL=${CHISEL_PHASER_PER_CELL:-100},PIN_CHISEL_HMM=${CHISEL_HMM:-none},PIN_CHISEL_ORACLE=${CHISEL_ORACLE:-false},PIN_CHISEL_ORACLE_SPLIT=${CHISEL_ORACLE_SPLIT:-0},PIN_CHISEL_MODE=${CHISEL_MODE:-pseudobulk},PIN_CHISEL_HMM_NOISE=${CHISEL_HMM_NOISE:-counting},PIN_CHISEL_RECLONE=${CHISEL_RECLONE:-false}"
# An --hmm or --oracle run re-does only the calling on the combo.tsv of the plain CHISEL run. When the route
# starts before the call stage, or that combo.tsv does not exist yet, the plain run is queued first (pinned
# without --hmm/--oracle, so it writes its own folder), and for --oracle also its eval (the lifted truth).
BASE_FIRST=false
if [[ ${CALLER} == chisel ]] && chisel_recall && runs call; then
    if runs phase; then
        BASE_FIRST=true
    elif [[ -z ${AFTER} && ! -s $(chisel_base_dir "${PHASER}")/combo/combo.tsv ]]; then
        BASE_FIRST=true
    fi
    if [[ ${BASE_FIRST} == true ]]; then
        BASE_PINS=$(echo "${PINS}" | sed -E 's/PIN_CHISEL_HMM=[^,]*/PIN_CHISEL_HMM=none/; s/PIN_CHISEL_ORACLE=[^,]*/PIN_CHISEL_ORACLE=false/; s/PIN_CHISEL_ORACLE_SPLIT=[^,]*/PIN_CHISEL_ORACLE_SPLIT=0/; s/PIN_CHISEL_RECLONE=[^,]*/PIN_CHISEL_RECLONE=false/')
        note "The plain CHISEL run $(chisel_base_dir "${PHASER}") is queued first$([[ ${CHISEL_ORACLE:-false} == true ]] && echo ", with its eval (the lifted truth for the oracle)")"
    fi
fi
if [[ -n ${CALLER} ]]; then
    CALLER_BP=$(caller_binsize "${CALLER}")
    CALLER_BP=${CALLER_BP%%_divclone}
    CALLER_BP=${CALLER_BP%%_pbclone}
    CALLER_BP=${CALLER_BP%%_hmm*}            # the bins, without an HMM run's suffix
    CALLER_BP=${CALLER_BP%%_oracle*}
    case ${CALLER_BP} in
        var*)    note "${CALLER} bins: variable (${CALLER_BP#var} SNP reads per cell)" ;;
        pbins*)  note "${CALLER} bins: ${PHASER}'s, ${CALLER_BP#pbins} SNP reads per cell ($(phaser_bins_file "${PHASER}" "${CALLER_BP#pbins}"))" ;;
        *)       note "${CALLER} bin size: ${CALLER_BP} bp" ;;
    esac
fi

# ============================================================
# Inputs for the first stage must exist unless we are waiting on --after jobs
# ============================================================

need() {   # path description fix
    [[ -e $1 ]] && return
    die "$2 not found: $1
  $3 (or pass --after JOBID if the job that makes it is still running)"
}

if [[ -z ${AFTER} ]]; then
    if chisel_reclone && runs call; then
        need "$(clone_base_dir "${PHASER}")/calls/calls.tsv" "calls of the legacy run" \
            "run it first: the same command with --mode legacy and without --reclone"
    elif chisel_recall && runs call && [[ ${BASE_FIRST} != true ]]; then
        need "$(chisel_base_dir "${PHASER}")/combo/combo.tsv" "combo.tsv of the run without --hmm/--oracle" \
            "run it first: the same command without --hmm/--oracle (--from call)"
        [[ ${CHISEL_ORACLE:-false} != true ]] || need "${EVAL_TRUTH}" "lifted truth" \
            "run any eval of this project first (it writes ${EVAL_TRUTH})"
    elif runs phase || { runs call && [[ ${CALLER} == chisel ]]; }; then
        ls "${CELL_BAM_DIR}"/*.bam >/dev/null 2>&1 \
            || die "no cell BAMs in ${CELL_BAM_DIR}; run ./run_sim.sh first"
        need "${BULK_BAM}" "bulk BAM" "run ./run_sim.sh first"
    fi
    if (( FROM_I >= 1 )) && ! chisel_reclone && { ! chisel_recall || [[ ${BASE_FIRST} == true ]]; }; then
        need "$(phaser_output "${PHASER}")" "${PHASER} output" "start from the phase stage (--from phase)"
        if [[ ${CALLER} == chisel && ${CHISEL_PHASER_BINS:-false} == true && ${PHASER} == ugp ]]; then
            need "$(phaser_bins_file ugp "${CHISEL_PHASER_PER_CELL}")" "UGP bins" "make them: --phaser ugp --to phase --phaserbins ${CHISEL_PHASER_PER_CELL} (UGP_MODE=pseudobulk), or start this run --from phase"
        fi
    fi
    if (( FROM_I >= 2 )); then
        need "$(caller_output "${CALLER}" "${PHASER}")" "${CALLER} output for ${PHASER}" "start from the call stage (--from call)"
    fi
fi

# ============================================================
# Submit
# ============================================================

mkdir -p logs   # the #SBATCH --output paths are relative to sim/

join_deps() { local out="" x; for x in "$@"; do [[ -n $x ]] && out+="${out:+:}$x"; done; echo "${out}"; }

submit() {   # label deps script [extra sbatch args...]; the jobs' settings are SUBMIT_PINS if set, else PINS
    local label=$1 deps=$2 script=$3; shift 3
    local args=(--parsable --export="ALL,PHASER=${PHASER},CALLER=${CALLER},${SUBMIT_PINS:-${PINS}}")
    [[ -n ${deps} ]] && args+=(--dependency="afterok:${deps}" --kill-on-invalid-dep=yes)
    if [[ ${DRY} == true ]]; then
        echo "  sbatch ${args[*]} $* ${script}" >&2
        echo "<${label}>"
    else
        sbatch "${args[@]}" "$@" "${script}"
    fi
}

ROUTE_DESC="${PHASER}${CALLER:+ -> ${CALLER}}, stages ${FROM}..${TO}"
[[ ${DRY} == true ]] && note "Dry run (${ROUTE_DESC}); would submit:" || note "Submitting ${ROUTE_DESC}"

PREV="${AFTER}"
SUMMARY=()

PREP=""
if runs phase; then
    if [[ ${PHASER} == ugp && ${UGP_MODE} == pseudobulk ]]; then
        # UGP pseudobulk mode reads chisel_prep's pooled BAM
        PREP=$(submit chisel_prep "${PREV}" scripts/chisel_prep.sh)
        SUMMARY+=("phase (chisel_prep for UGP pseudobulk): ${PREP}")
        PREV=$(join_deps "${PREV}" "${PREP}")
    fi
    J=$(submit "phase_${PHASER}" "${PREV}" "scripts/phase_${PHASER}.sh")
    SUMMARY+=("phase (${PHASER}): ${J}")
    PREV=${J}
fi

if runs call; then
    case ${CALLER} in
        hiscanner)
            J=$(submit call_hiscanner "${PREV}" scripts/call_hiscanner.sh)
            SUMMARY+=("call (hiscanner): ${J}")
            ;;
        chisel)
          if chisel_reclone; then
            # re-cloning: only the cloning, on a copy of the calls of the legacy run
            J=$(submit clone_chisel "${PREV}" scripts/clone_chisel.sh)
            SUMMARY+=("call (chisel re-cloning, ${CHISEL_MODE:-pseudobulk}): ${J}")
          elif chisel_recall; then
            # HMM or oracle clusters: only the calling, on the combo.tsv of the run without them
            if [[ ${BASE_FIRST} == true ]]; then
                SUBMIT_PINS=${BASE_PINS}
                P=${PREP:-$(submit chisel_prep "${AFTER}" scripts/chisel_prep.sh)}
                SUMMARY+=("call (chisel_prep): ${P}")
                DEPS=$(join_deps "${PREV}" "${P}")
                if [[ ${CHISEL_PHASER_BINS:-false} == true && ${PHASER} == perfect ]]; then
                    B=$(submit bins_perfect "${DEPS}" scripts/bins_perfect.sh)
                    SUMMARY+=("call (bins_perfect): ${B}")
                    DEPS=$(join_deps "${DEPS}" "${B}")
                fi
                BJ=$(submit call_chisel "${DEPS}" scripts/call_chisel.sh)
                SUMMARY+=("call (plain chisel, makes combo.tsv): ${BJ}")
                PREV=${BJ}
                if [[ ${CHISEL_ORACLE:-false} == true ]]; then
                    BE=$(submit eval "${BJ}" scripts/eval.sh)
                    SUMMARY+=("eval (plain chisel, writes the lifted truth): ${BE}")
                    PREV=$(join_deps "${BJ}" "${BE}")
                fi
                unset SUBMIT_PINS
            fi
            J=$(submit call_chisel_hmm "${PREV}" scripts/call_chisel_hmm.sh)
            SUMMARY+=("call (chisel $([[ ${CHISEL_ORACLE:-false} == true ]] && echo "oracle clusters" || echo "HMM ${CHISEL_HMM}")): ${J}")
          else
            # The barcoded BAM needs only the cell BAMs, so it runs alongside phasing
            P=${PREP:-$(submit chisel_prep "${AFTER}" scripts/chisel_prep.sh)}
            SUMMARY+=("call (chisel_prep): ${P}")
            DEPS=$(join_deps "${PREV}" "${P}")
            if [[ ${CHISEL_PHASER_BINS:-false} == true && ${PHASER} == perfect ]]; then
                # bins from the true het SNPs: needs phase_perfect's allele depths and the pooled BAM
                B=$(submit bins_perfect "${DEPS}" scripts/bins_perfect.sh)
                SUMMARY+=("call (bins_perfect): ${B}")
                DEPS=$(join_deps "${DEPS}" "${B}")
            fi
            J=$(submit call_chisel "${DEPS}" scripts/call_chisel.sh)
            SUMMARY+=("call (chisel): ${J}")
          fi
            ;;
    esac
    PREV=${J}
fi

if runs eval; then
    J=$(submit eval "${PREV}" scripts/eval.sh)
    SUMMARY+=("eval: ${J}")
    PREV=${J}
fi

for s in "${SUMMARY[@]}"; do note "  ${s}"; done
if [[ ${DRY} != true ]]; then
    runs phase && note "Phasing: $(dirname "$(phaser_output "${PHASER}")")"
    runs call && note "Calls:   $(call_dir "${CALLER}" "${PHASER}")"
    runs eval && note "Results: $(eval_out_dir "${PHASER}" "${CALLER}")"
fi
echo "${PREV}"      # last job ID, for --after in a later command
