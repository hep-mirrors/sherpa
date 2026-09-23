#!/bin/bash
#SBATCH -J yfsww-int
#SBATCH -p plgrid
#SBATCH -A plgsherpalumi-cpu
#SBATCH -t 24:00:00
#SBATCH -x ac[0770-0780],ac0714
#SBATCH -o logs/integrate_%j.out
#SBATCH -e logs/integrate_%j.err
#
# Build the Comix process setup and the YFS/OpenLoops integration grids for ONE
# card, into int-res/<CARD> (+ proc/<CARD>/Process). Every later array task
# reads both. Run once per card before submitting its array, otherwise all N
# tasks redo the identical integration -- and for L6 that integration is
# several minutes of OpenLoops setup each.
#
# MPI is used here and ONLY here. Integration cannot be split by seed, so it is
# the one stage where extra ranks buy something; for event generation each rank
# would pay the same ~60 s initialisation an array task pays, so MPI there costs
# core-hours for nothing.
#
#   sbatch --ntasks=16 --export=ALL,CARD=L6_nlo_161p0_mue scripts/integrate.sbatch
#
# Integration also prints sigma with its own error before any event is
# generated. That number is a free cross-check of every array task that follows,
# and for the cheap rungs it is already at the precision the budget needs.
set -euo pipefail

source ./scripts/job_config.sh
source ./scripts/card_sets.sh

CARD=${CARD:?set CARD=...}
card_parse "${CARD}" >/dev/null
NTASKS=${SLURM_NTASKS:-1}
TOP=${SLURM_SUBMIT_DIR:-${PWD}}

assert_sherpa_has_ladder_weights

CARDFILE="${TOP}/cards/generated/${CARD}.yaml"
RES="${TOP}/int-res/${CARD}"
CPP="${TOP}/proc/${CARD}"
# integrate-all.sh and smoke.sh override -o to logs/integrate_<CARD>_<jobid>.out
# so the log says which card it belongs to. The header default has no card in
# it. Both are accepted, but a log this cannot find is a hard error: the whole
# point of the file is the assertions read off it, and skipping them quietly
# would leave a run that was never checked looking like one that passed.
LOG=""
for cand in "${TOP}/logs/integrate_${CARD}_${SLURM_JOB_ID:-0}.out" \
            "${TOP}/logs/integrate_${SLURM_JOB_ID:-0}.out"; do
    [[ -f ${cand} ]] && { LOG=${cand}; break; }
done
mkdir -p "${RES}" "${CPP}"

# OL_PREFIX only for the rung that actually calls OpenLoops. Passing it to a
# card that does not use it would show up in Sherpa's "settings that have not
# been used" warning, and check_unused_settings would then have to whitelist
# it -- which is exactly the kind of blanket exemption that later hides a real
# dropped setting.
OLARG=()
[[ $(card_rung "${CARD}") == L6_nlo ]] && OLARG=( "OL_PREFIX:${OL_PREFIX}" )

echo "### integrating ${CARD} on ${NTASKS} core(s)"
echo "### card ${CARDFILE}"
if [[ ${NTASKS} -gt 1 ]]; then
    mpirun -n "${NTASKS}" "${sherpa}" -f "${CARDFILE}" -e 0 \
        "RESULT_DIRECTORY:${RES}" "SHERPA_CPP_PATH:${CPP}" \
        ${OLARG[@]+"${OLARG[@]}"} BATCH_MODE:2
else
    "${sherpa}" -f "${CARDFILE}" -e 0 \
        "RESULT_DIRECTORY:${RES}" "SHERPA_CPP_PATH:${CPP}" \
        ${OLARG[@]+"${OLARG[@]}"} BATCH_MODE:2
fi

# Sherpa writes the result directory out as int-res/<CARD>.zip and leaves
# int-res/<CARD>/ empty, so a `-d` test passes on a bare mkdir. Accept either
# form but insist one of them is non-trivial.
if [[ ! -f "${RES}.zip" ]] && [[ -z $(ls -A "${RES}" 2>/dev/null || true) ]]; then
    echo "ERROR: integration produced neither ${RES}.zip nor a non-empty ${RES}/" >&2
    exit 1
fi

if [[ -z ${LOG} ]]; then
    echo "ERROR: cannot find this job's log, so the setting and cross-section" >&2
    echo "       checks did not run. Submit via scripts/integrate-all.sh." >&2
    exit 1
fi
check_unused_settings "${LOG}"

# Stamp the grids with the exact generator that made them. sherpa_array.sbatch
# refuses to run against grids from a different build.
sherpa_fingerprint > "${TOP}/int-res/${CARD}.provenance"
echo "### generator stamp:"
sed 's/^/    /' "${TOP}/int-res/${CARD}.provenance"

echo "### integrated ${CARD}"
echo "### integration cross section -- a free cross-check of every array task"
echo "### that follows, and for the cheap rungs already at the target precision:"
grep -E "^[0-9].*: .* pb \+-" "${LOG}" | tail -5 || true
