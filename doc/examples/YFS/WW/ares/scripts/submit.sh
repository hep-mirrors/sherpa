#! /bin/bash
# Submit the event-generation array for ONE card, with its merge chained on.
#
#   ./scripts/submit.sh L6_nlo_161p0_mue              # per-rung default seeds
#   ./scripts/submit.sh L6_nlo_161p0_mue 20           # 20 seeds
#   EVT=50000 ./scripts/submit.sh L6_nlo_161p0_mue 20
#   ./scripts/submit.sh L6_nlo_161p0_mue 10 --continue   # 10 MORE seeds
#
# SEEDS START AT 1 AND ARE THE SAME FOR EVERY CARD, deliberately, and that is a
# change from the pion campaign's continue-from-the-highest-seed default. Two of
# the comparisons here are differences between cards -- the three m_W points of
# dsigma/dm_W, and rung-to-rung at fixed energy and channel -- and running them
# on the same seed set makes part of the MC fluctuation common and cancel.
# Topping up one card without the others silently breaks that alignment, so a
# top-up has to be asked for explicitly and should be applied to a whole family.
#
# The seed is both the RNG seed and the output filename, so re-running a seed
# regenerates the identical events AND overwrites the file. Existing seeds are
# therefore a hard error unless --continue (append a fresh range) or --force
# (redo the same range) is given.
#
# Integration must have run first: the array tasks read int-res/<CARD> and
# proc/<CARD>. This refuses rather than let N tasks each redo the integration.
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

CONTINUE=0; FORCE=0
args=()
for a in "$@"; do
    case "${a}" in
        --continue) CONTINUE=1 ;;
        --force)    FORCE=1 ;;
        -*) echo "ERROR: unknown flag '${a}'" >&2; exit 1 ;;
        *)  args+=( "${a}" ) ;;
    esac
done
set -- ${args[@]+"${args[@]}"}

CARD=${1:?usage: submit.sh <CARD> [NSEED] [--continue|--force]}
card_parse "${CARD}" >/dev/null
NSEED=${2:-$(nseed_for "${CARD}")}
EVT=${EVT:-$(evt_for "${CARD}")}
WALL=${WALL:-$(walltime_for "${CARD}")}

CARDFILE="cards/generated/${CARD}.yaml"
[[ -f ${CARDFILE} ]] || { echo "ERROR: no such runcard: ${CARDFILE}" >&2; exit 1; }

# The card, not the command line, decides the generation mode -- but a card that
# lost the line would run PartiallyUnweighted at an unweighting efficiency of
# 0.07-0.7% and produce a few hundred events where we asked for 200 000, which
# looks like a slow node rather than a mistake. Assert instead of overriding.
grep -qE '^EVENT_GENERATION_MODE:[[:space:]]*Weighted[[:space:]]*$' "${CARDFILE}" || {
    echo "ERROR: ${CARDFILE} does not set 'EVENT_GENERATION_MODE: Weighted'." >&2
    echo "       Regenerate the cards with cards/make-cards.py." >&2; exit 1; }
grep -q 'yfs-ww' "${CARDFILE}" || {
    echo "ERROR: ${CARDFILE} does not request the yfs-ww Rivet analysis." >&2; exit 1; }

if [[ ! -f "int-res/${CARD}.zip" ]] \
   && [[ -z $(ls -A "int-res/${CARD}" 2>/dev/null || true) ]]; then
    echo "ERROR: no integration grids for ${CARD}" >&2
    echo "       (looked for int-res/${CARD}.zip and a non-empty int-res/${CARD}/)" >&2
    echo "       integrate first:" >&2
    echo "       sbatch --ntasks=$(intcores_for "${CARD}") --export=ALL,CARD=${CARD} scripts/integrate.sbatch" >&2
    exit 1
fi

mkdir -p "logs/${CARD}" "yodas/${CARD}" "$(merged_dir "${CARD}")"

existing=$(ls -1 "yodas/${CARD}"/*.yoda.gz 2>/dev/null | wc -l | tr -d ' ' || true)
if [[ ${CONTINUE} -eq 1 ]]; then
    last=$(ls -1 "yodas/${CARD}"/*.yoda.gz 2>/dev/null \
             | sed -E 's|.*/([0-9]+)\.yoda\.gz|\1|' | sort -n | tail -1 || true)
    FIRST=$(( ${last:-0} + 1 ))
elif [[ ${existing} -gt 0 && ${FORCE} -ne 1 ]]; then
    echo "ERROR: yodas/${CARD} already holds ${existing} seed(s)." >&2
    echo "       Seeds 1-${NSEED} would be REGENERATED IDENTICALLY and overwritten." >&2
    echo "       --continue to add a fresh range, --force to redo, or" >&2
    echo "       ./scripts/reset.sh ${CARD} --yes for a clean run." >&2
    exit 1
else
    FIRST=1
fi
LAST=$(( FIRST + NSEED - 1 ))

maxidx=$(scontrol show config 2>/dev/null | awk '/^MaxArraySize/ {print $3}')
if [[ -n ${maxidx} && ${LAST} -ge ${maxidx} ]]; then
    echo "ERROR: seed range ${FIRST}-${LAST} exceeds MaxArraySize=${maxidx}." >&2
    exit 1
fi

# MaxSubmitJobs is 1000 for pending+running and CHAINED MERGES COUNT. Each card
# costs NSEED+1.
queued=$(squeue -u "${USER}" -h -t pending,running -r 2>/dev/null | wc -l | tr -d ' ')
if [[ $(( queued + NSEED + 1 )) -gt 1000 ]]; then
    echo "ERROR: ${queued} job(s) already queued; +${NSEED} array tasks +1 merge" >&2
    echo "       would exceed MaxSubmitJobs=1000. Wait, or submit a smaller set." >&2
    exit 1
fi

# Log file per SEED, not per job id: the seed already identifies the sample, and
# %A_%a would put a resubmitted seed in a second file so that logs/<CARD>/7.*
# stops meaning "everything about seed 7". The array script reads its own .out
# back to assert the Seed: banner, so the name has to be predictable.
jid=$(sbatch --parsable -A "${ACCOUNT}" -p "${PARTITION}" -x "${EXCLUDE_NODES}" \
        --export=ALL,CARD="${CARD}",EVT="${EVT}" \
        --array="${FIRST}"-"${LAST}" \
        --time="${WALL}" \
        -o "logs/${CARD}/%a.out" \
        -e "logs/${CARD}/%a.err" \
        scripts/sherpa_array.sbatch)
echo "array ${jid}: ${CARD}, seeds ${FIRST}-${LAST} (${NSEED} x ${EVT} ev, t=${WALL})"

# afterany, not afterok: afterok needs EVERY array task to exit 0, so one bad
# seed out of ten would cancel the merge. merge.sbatch refuses an empty set and
# is happy with a partial one -- nine seeds is unbiased, just less precise.
mid=$(sbatch --parsable -A "${ACCOUNT}" -p "${PARTITION}" -x "${EXCLUDE_NODES}" \
        --dependency=afterany:"${jid}" \
        --export=ALL,CARD="${CARD}" \
        -o "logs/${CARD}/merge_%j.out" \
        -e "logs/${CARD}/merge_%j.err" \
        scripts/merge.sbatch)
echo "merge ${mid}: after ${jid} -> $(merged_dir "${CARD}")/${CARD}.yoda.gz"
