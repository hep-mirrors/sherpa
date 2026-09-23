#! /bin/bash
# One-screen state of a SET: seeds published, merge present, queue.
#
#   ./scripts/status.sh all
#   ./scripts/status.sh l6
#
# Exists so that monitoring is one short ssh call. A long chain of remote
# commands gets SIGHUP'd part way through, and the half that ran looks exactly
# like the half that did not.
set -uo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh >/dev/null
source ./scripts/card_sets.sh

SET=${1:-all}
resolve_set "${SET}" || exit 1

printf "%-26s %5s %6s %7s %6s  %s\n" card grids seeds merged fail state
tot_seeds=0; tot_done=0
for c in "${cards[@]}"; do
    grids=no
    { [[ -f "int-res/${c}.zip" ]] || [[ -n $(ls -A "int-res/${c}" 2>/dev/null || true) ]]; } && grids=yes
    n=$(ls -1 "yodas/${c}"/*.yoda.gz 2>/dev/null | wc -l | tr -d ' ')
    m="-"
    [[ -f "$(merged_dir "${c}")/${c}.yoda.gz" ]] && m="yes"
    # A task that ended without reaching the OK line. Counted separately from
    # "no seed" because a failed seed still leaves a log.
    logs=$(ls -1 "logs/${c}"/[0-9]*.out 2>/dev/null | wc -l | tr -d ' ')
    fail=$(( logs - n )); [[ ${fail} -lt 0 ]] && fail=0
    st=$(squeue -u "${USER}" -h -o "%T" -n "yfsww-gen" 2>/dev/null | sort -u | tr '\n' ',' )
    printf "%-26s %5s %6s %7s %6s  %s\n" "${c}" "${grids}" "${n}" "${m}" "${fail}" "${st%,}"
    tot_seeds=$(( tot_seeds + n ))
    [[ ${m} == yes ]] && tot_done=$(( tot_done + 1 ))
done
echo
echo "${tot_seeds} seed(s) published, ${tot_done}/${#cards[@]} card(s) merged"
q=$(squeue -u "${USER}" -h -t pending,running -r 2>/dev/null | wc -l | tr -d ' ')
echo "${q} job(s) in the queue (cap is MaxSubmitJobs=1000, chained merges included)"
