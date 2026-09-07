#! /bin/bash
# Integrate every card in a SET once, into int-res/<CARD> and proc/<CARD>.
# Nothing can be submitted until this has finished for the cards involved.
#
#   ./scripts/integrate-all.sh all
#   ./scripts/integrate-all.sh l6            # only the expensive rung
#   ./scripts/integrate-all.sh rung:L6_nlo,energy:161p0
#
# Core count comes from intcores_for() in job_config.sh: 16 for L6_nlo, which
# integrates the BVR-matched process through OpenLoops, 4 for everything else.
# MPI is worth it here and nowhere else -- integration cannot be split by seed.
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

SET=${1:?usage: integrate-all.sh <SET>   (aliases: $(set_names))}
resolve_set "${SET}"

mkdir -p logs
echo "### ${#cards[@]} card(s) in set '${SET}'"
for c in "${cards[@]}"; do
    n=$(intcores_for "${c}")
    jid=$(sbatch --parsable -A "${ACCOUNT}" -p "${PARTITION}" -x "${EXCLUDE_NODES}" --ntasks="${n}" --export=ALL,CARD="${c}" \
            -o "logs/integrate_${c}_%j.out" -e "logs/integrate_${c}_%j.err" \
            scripts/integrate.sbatch)
    echo "integrate ${jid}: ${c} on ${n} core(s)"
done
echo
echo "wait for these, then: ./scripts/submit-all.sh ${SET}"
