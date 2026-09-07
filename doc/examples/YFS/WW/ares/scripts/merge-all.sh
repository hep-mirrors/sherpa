#! /bin/bash
# Re-merge every card in a SET that has yodas. Normally unnecessary -- submit.sh
# chains a merge onto every array -- but needed after a manual top-up, after a
# merge job was killed, or when seeds were added with --continue.
#
#   ./scripts/merge-all.sh all
#   ./scripts/merge-all.sh l6
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

SET=${1:?usage: merge-all.sh <SET>   (aliases: $(set_names))}
resolve_set "${SET}"

any=0
for c in "${cards[@]}"; do
    n=$(ls -1 "yodas/${c}"/*.yoda.gz 2>/dev/null | wc -l | tr -d ' ' || true)
    if [[ ${n} -eq 0 ]]; then
        echo "skip ${c}: no yodas"
        continue
    fi
    mkdir -p "logs/${c}"
    jid=$(sbatch --parsable -A "${ACCOUNT}" -p "${PARTITION}" -x "${EXCLUDE_NODES}" --export=ALL,CARD="${c}" \
            -o "logs/${c}/merge_%j.out" -e "logs/${c}/merge_%j.err" \
            scripts/merge.sbatch)
    echo "merge ${jid}: ${c} (${n} seeds) -> $(merged_dir "${c}")/${c}.yoda.gz"
    any=1
done
[[ ${any} -eq 0 ]] && echo "nothing to merge in set '${SET}'"
exit 0
