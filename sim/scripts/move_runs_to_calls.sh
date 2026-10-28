#!/bin/bash
# move_runs_to_calls.sh -- move a project's run folders (chisel_*, hiscanner_*) into ${PROJECT_DIR}/calls/.
#
#   scripts/move_runs_to_calls.sh            show what would move (dry run)
#   scripts/move_runs_to_calls.sh --apply    move them
#
# Run folders used to sit at the top of the project (chisel_pbins8_hmm9/, hiscanner_500000/, ...); the
# pipeline now writes them to ${CALLS_DIR} (default ${PROJECT_DIR}/calls). This moves existing ones there.
# Symlinks inside the moved folders that point into another moved folder (the combo/ link of --hmm,
# --oracle and re-cloned runs) are rewritten as relative links to the new place. eval/ and chisel/prep/
# stay where they are. Refuses while you have SLURM jobs queued or running (they would look for the
# folders mid-move); --force skips that check. Uses CONFIG_FILE like the pipeline (default ./config.sh).
set -euo pipefail
cd "$(dirname "$0")/.."
APPLY=false FORCE=false
for a in "$@"; do
    case $a in
        --apply) APPLY=true ;;
        --force) FORCE=true ;;
        -h|--help) sed -n 2,13p "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) echo "unknown option $a (--apply, --force)" >&2; exit 1 ;;
    esac
done
source "${CONFIG_FILE:-config.sh}"
source scripts/common.sh
[[ -d ${PROJECT_DIR} ]] || die "no project folder ${PROJECT_DIR}"
mkdir -p "${CALLS_DIR}"

runs=()
for d in "${PROJECT_DIR}"/chisel_* "${PROJECT_DIR}"/hiscanner_*; do
    [[ -d $d && ! -L $d ]] || continue
    runs+=("$(basename "$d")")
done
if (( ${#runs[@]} == 0 )); then
    note "Nothing to move: no chisel_* or hiscanner_* folders at the top of ${PROJECT_DIR}"
    exit 0
fi
clash=()
for r in "${runs[@]}"; do [[ -e ${CALLS_DIR}/${r} ]] && clash+=("${r}"); done
(( ${#clash[@]} == 0 )) || die "already in ${CALLS_DIR}: ${clash[*]} (move or remove those by hand first)"

note "${#runs[@]} run folders -> ${CALLS_DIR}/:"
printf '    %s\n' "${runs[@]}" >&2
if [[ ${APPLY} != true ]]; then
    note "Dry run: nothing moved (--apply to move them)"
    exit 0
fi
if [[ ${FORCE} != true ]] && command -v squeue >/dev/null; then
    n=$(squeue -h -u "${USER}" 2>/dev/null | wc -l)
    (( n == 0 )) || die "${n} of your SLURM jobs are queued or running; wait for them (or --force)"
fi

for r in "${runs[@]}"; do
    mv "${PROJECT_DIR}/${r}" "${CALLS_DIR}/${r}"
done

# links into moved folders: absolute targets under PROJECT_DIR/<moved run>/... -> relative to the new place
fixed=0 broken=()
while IFS= read -r -d '' link; do
    t=$(readlink "${link}")
    new=""
    for r in "${runs[@]}"; do
        case $t in
            "${PROJECT_DIR}/${r}"|"${PROJECT_DIR}/${r}"/*) new="${CALLS_DIR}/${r}${t#"${PROJECT_DIR}/${r}"}" ;;
        esac
    done
    if [[ -n ${new} ]]; then
        rel_link "${new}" "${link}"
        fixed=$((fixed + 1))
    fi
    [[ -e ${link} ]] || broken+=("${link#"${CALLS_DIR}/"} -> $(readlink "${link}")")
done < <(find "${CALLS_DIR}" -type l -print0)

note "Moved ${#runs[@]} folders; rewrote ${fixed} links"
if (( ${#broken[@]} )); then
    note "WARNING: links that do not resolve (they pointed elsewhere before the move too, or at removed runs):"
    printf '    %s\n' "${broken[@]}" >&2
fi
