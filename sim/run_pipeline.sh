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
  --after IDS    wait for these SLURM job IDs (colon- or comma-separated) before
                 starting; run_sim.sh passes its job IDs here
  --dry-run      print the sbatch commands without submitting
  --list         show implemented routes and exit

Stages
  phase   scan2: SCAN2 (GATK + Eagle)       ugp: Universal Genotyping Pipeline (bulk mode)
  call    hiscanner: HiScanner              chisel: chisel_prep (barcoded BAM) + CHISEL
  eval    lift truth to hg38, convert calls to a common table, score them

Examples
  ./run_pipeline.sh --phaser ugp --caller chisel                # all three stages
  ./run_pipeline.sh --phaser ugp --caller chisel --from call    # UGP already done
  ./run_pipeline.sh --phaser scan2 --caller hiscanner --only eval
  ./run_pipeline.sh --phaser ugp --to phase                     # just phasing
  ./run_sim.sh -- --phaser ugp --caller chisel                  # simulate, then all of this

Outputs: ${PROJECT_DIR}/<caller>/<phaser>/ and ${EVAL_DIR}/<phaser>_<caller>/
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

PHASER="" CALLER="" FROM=phase TO=eval AFTER="" DRY=false
while (( $# )); do
    case $1 in
        --phaser)  PHASER=${2:?--phaser needs a value}; shift 2 ;;
        --caller)  CALLER=${2:?--caller needs a value}; shift 2 ;;
        --from)    FROM=${2:?--from needs a value};     shift 2 ;;
        --to)      TO=${2:?--to needs a value};         shift 2 ;;
        --only)    FROM=${2:?--only needs a value}; TO=$2; shift 2 ;;
        --after)   AFTER=${2:?--after needs a value};   shift 2 ;;
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

# ============================================================
# Inputs for the first stage must exist unless we are waiting on --after jobs
# ============================================================

need() {   # path description fix
    [[ -e $1 ]] && return
    die "$2 not found: $1
  $3 (or pass --after JOBID if the job that makes it is still running)"
}

if [[ -z ${AFTER} ]]; then
    if runs phase || { runs call && [[ ${CALLER} == chisel ]]; }; then
        ls "${CELL_BAM_DIR}"/*.bam >/dev/null 2>&1 \
            || die "no cell BAMs in ${CELL_BAM_DIR}; run ./run_sim.sh first"
        need "${BULK_BAM}" "bulk BAM" "run ./run_sim.sh first"
    fi
    if (( FROM_I >= 1 )); then
        need "$(phaser_output "${PHASER}")" "${PHASER} output" "start from the phase stage (--from phase)"
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

submit() {   # label deps script [extra sbatch args...]
    local label=$1 deps=$2 script=$3; shift 3
    local args=(--parsable --export="ALL,PHASER=${PHASER},CALLER=${CALLER}")
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

if runs phase; then
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
            # The barcoded BAM needs only the cell BAMs, so it runs alongside phasing
            P=$(submit chisel_prep "${AFTER}" scripts/chisel_prep.sh)
            J=$(submit call_chisel "$(join_deps "${PREV}" "${P}")" scripts/call_chisel.sh)
            SUMMARY+=("call (chisel_prep): ${P}" "call (chisel): ${J}")
            ;;
    esac
    PREV=${J}
fi

if runs eval; then
    J=$(submit eval "${PREV}" scripts/eval.sh)
    SUMMARY+=("eval: ${J}")
fi

for s in "${SUMMARY[@]}"; do note "  ${s}"; done
if [[ ${DRY} != true ]]; then
    runs call && note "Calls:   $(call_dir "${CALLER}" "${PHASER}")"
    runs eval && note "Results: $(eval_out_dir "${PHASER}" "${CALLER}")"
fi
