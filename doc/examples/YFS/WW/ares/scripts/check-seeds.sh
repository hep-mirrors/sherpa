#! /bin/bash
# Assert that every array task of a card really used a different RNG seed.
#
#   ./scripts/check-seeds.sh L6_nlo_161p0_mue
#   ./scripts/check-seeds.sh all
#
# The seed flag is -R. `-s` is accepted and does nothing, and the banner then
# reads 'Seed: 1' for every task: N identical samples, merged as if they were
# independent, giving a cross-section error smaller than the truth by sqrt(N).
# sherpa_array.sbatch asserts the banner matches its own -R, but only per task;
# this is the across-the-array check, and it also catches the case where two
# tasks were handed the same seed by a bad --array range.
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

SET=${1:?usage: check-seeds.sh <SET|CARD>}
resolve_set "${SET}"

bad=0
for c in "${cards[@]}"; do
    logs=( $(ls -1 "logs/${c}"/[0-9]*.out 2>/dev/null || true) )
    if [[ ${#logs[@]} -eq 0 ]]; then
        echo "skip ${c}: no task logs"
        continue
    fi
    seen=""
    for l in "${logs[@]}"; do
        want=$(basename "${l}" .out)
        got=$(grep -m1 -E '^Seed: ' "${l}" | awk '{print $2}' || true)
        if [[ -z ${got} ]]; then
            echo "  ${c} seed ${want}: NO 'Seed:' LINE (job may still be starting)" >&2
            bad=1; continue
        fi
        if [[ ${got} != "${want}" ]]; then
            echo "  ${c} seed ${want}: banner says 'Seed: ${got}'" >&2
            bad=1
        fi
        seen="${seen} ${got}"
    done
    uniq=$(echo "${seen}" | tr ' ' '\n' | grep -c . || true)
    dist=$(echo "${seen}" | tr ' ' '\n' | grep . | sort -u | wc -l | tr -d ' ')
    if [[ ${uniq} -ne ${dist} ]]; then
        echo "  ${c}: ${uniq} tasks but only ${dist} DISTINCT seeds" >&2
        bad=1
    else
        echo "  ok ${c}: ${dist} distinct seeds over ${uniq} tasks"
    fi
done
[[ ${bad} -eq 0 ]] || { echo "SEED CHECK FAILED" >&2; exit 1; }
echo "seed check passed"
