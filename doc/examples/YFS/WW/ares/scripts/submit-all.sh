#! /bin/bash
# Submit a whole SET: one array plus one chained merge per card.
#
#   ./scripts/submit-all.sh cheap            # everything except L6
#   ./scripts/submit-all.sh l6
#   ./scripts/submit-all.sh all
#   ./scripts/submit-all.sh l6 --clean       # wipe the set first
#   ./scripts/submit-all.sh l6 5             # 5 seeds per card
#
# Without --clean this refuses any card that already has yodas (see submit.sh):
# topping up is not the default here because several of the comparisons are
# differences between cards run on a common seed set, and extending one card
# without the others quietly breaks the cancellation.
#
# Prints the job count against MaxSubmitJobs=1000 BEFORE submitting anything.
# Chained merges count towards that limit, so the real cost of a card is
# NSEED+1 jobs.
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

CLEAN=0; PASS=()
args=()
for a in "$@"; do
    case "${a}" in
        --clean)    CLEAN=1 ;;
        --continue) PASS+=( --continue ) ;;
        --force)    PASS+=( --force ) ;;
        -*) echo "ERROR: unknown flag '${a}'" >&2; exit 1 ;;
        *)  args+=( "${a}" ) ;;
    esac
done
set -- ${args[@]+"${args[@]}"}

SET=${1:?usage: submit-all.sh <SET> [NSEED] [--clean] [--continue|--force]}
NSEED_OVERRIDE=${2:-}
resolve_set "${SET}"

# Cost preview. Everything here is one core per task, so core-hours are
# tasks x wall time, and the wall time below is the REQUEST, not a prediction.
total=0
echo "### set '${SET}': ${#cards[@]} card(s)"
for c in "${cards[@]}"; do
    n=${NSEED_OVERRIDE:-$(nseed_for "${c}")}
    total=$(( total + n + 1 ))
done
queued=$(squeue -u "${USER}" -h -t pending,running -r 2>/dev/null | wc -l | tr -d ' ')
echo "### ${total} SLURM jobs (arrays + chained merges), ${queued} already queued"
if [[ $(( total + queued )) -gt 1000 ]]; then
    echo "ERROR: that would exceed MaxSubmitJobs=1000 (pending+running, merges included)." >&2
    echo "       Submit in stages: 'cheap' first, then 'l6'." >&2
    exit 1
fi

# Pre-flight the whole set BEFORE submitting anything. submit.sh refuses a card
# that already has yodas, and finding that out card 40 of 81 leaves half a set
# queued and half not -- the state that is hardest to reason about later.
if [[ ${CLEAN} -eq 0 ]]; then
    dirty=()
    for c in "${cards[@]}"; do
        n=$(ls -1 "yodas/${c}"/*.yoda.gz 2>/dev/null | wc -l | tr -d ' ' || true)
        [[ ${n} -gt 0 ]] && dirty+=( "${c}(${n})" )
    done
    if [[ ${#dirty[@]} -gt 0 && " ${PASS[*]-} " != *--continue* && " ${PASS[*]-} " != *--force* ]]; then
        echo "ERROR: ${#dirty[@]} card(s) in '${SET}' already have yodas:" >&2
        printf '         %s\n' "${dirty[@]}" >&2
        echo "       Those seeds would be regenerated identically and overwritten." >&2
        echo "       --clean for a fresh run, --continue to add seeds, --force to redo." >&2
        echo "       (A smoke test leaves its four cards populated; clearing them with" >&2
        echo "        ./scripts/reset.sh smoke --yes is the normal thing to do next.)" >&2
        exit 1
    fi
fi

if [[ ${CLEAN} -eq 1 ]]; then
    echo "### --clean: clearing '${SET}' first"
    ./scripts/reset.sh "${SET}" --yes
    echo
fi

for c in "${cards[@]}"; do
    ./scripts/submit.sh "${c}" ${NSEED_OVERRIDE:+"${NSEED_OVERRIDE}"} ${PASS[@]+"${PASS[@]}"}
done

echo
echo "queued. watch with: squeue -u \$USER"
echo "when the merges are done:  ./scripts/link-mw.sh  then fetch-from-ares.sh"
