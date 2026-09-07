#! /bin/bash
# Generate the run cards on Ares and lay out the campaign directory.
# Idempotent; run it after every push-to-ares.sh.
#
#   ssh ares 'cd /net/afscra/people/plgaprice/yfs-ww && ./scripts/prepare.sh'
#
# The cards are not rsynced: cards/make-cards.py is the source and the 81 yaml
# files are its deterministic output. Generating them here rather than shipping
# them means the cluster can never be running a card set that no longer matches
# the generator -- the failure mode being a hand-edited card on the cluster that
# nothing off it reproduces, i.e. a sample with no provenance.
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

[[ -f cards/make-cards.py ]] || {
    echo "ERROR: no cards/make-cards.py -- run push-to-ares.sh first" >&2; exit 1; }

echo "### generating run cards"
python3 cards/make-cards.py

mkdir -p logs yodas int-res proc merged work analysis

n=$(ls -1 cards/generated/*.yaml 2>/dev/null | wc -l | tr -d ' ')

# Cross-check the generator against the card matrix job_config.sh filters on.
# If make-cards.py grows a rung, an energy or a channel and job_config.sh does
# not, resolve_set silently stops matching the new family: integrate-all and
# submit-all skip it, reset.sh does not clean it, and nothing errors.
nr=$(echo "${RUNGS_ALL}"    | wc -w | tr -d ' ')
ne=$(echo "${ENERGIES_ALL}" | wc -w | tr -d ' ')
nc=$(echo "${CHANNELS_ALL}" | wc -w | tr -d ' ')
nm=$(echo "${MASSES_ALL}"   | wc -w | tr -d ' ')
expected=$(( (nr - 1) * ne * nc + ne * nm ))     # MW is the rung with its own count
if [[ ${n} -ne ${expected} ]]; then
    echo "ERROR: make-cards.py wrote ${n} cards, job_config.sh's matrix implies ${expected}." >&2
    echo "       ((${nr}-1) rungs x ${ne} energies x ${nc} channels) + (${ne} x ${nm} m_W probes)" >&2
    echo "       Bring RUNGS_ALL/ENERGIES_ALL/CHANNELS_ALL/MASSES_ALL back in step" >&2
    echo "       with cards/make-cards.py before running anything." >&2
    exit 1
fi

# Every card name must parse, or it lands in the wrong merged directory.
resolve_set all >/dev/null
echo "### ${n} cards, all parsed, matrix agrees with job_config.sh"

# Same assertion submit.sh makes, made once here for the whole set rather than
# 81 times at submit.
bad=0
for f in cards/generated/*.yaml; do
    grep -qE '^EVENT_GENERATION_MODE:[[:space:]]*Weighted[[:space:]]*$' "${f}" \
        || { echo "ERROR: $(basename "${f}") is not EVENT_GENERATION_MODE: Weighted" >&2; bad=1; }
done
[[ ${bad} -eq 0 ]] || exit 1
echo "### every card is weighted"

echo
echo "next: ./scripts/build-analysis.sh"
