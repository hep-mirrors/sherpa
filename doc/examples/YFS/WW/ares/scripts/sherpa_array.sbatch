#!/bin/bash
#SBATCH -J yfsww-gen
#SBATCH -p plgrid
#SBATCH -A plgsherpalumi-cpu
#SBATCH -t 24:00:00
#SBATCH --ntasks=1
#SBATCH -x ac[0770-0780],ac0714
#
# ONE core deliberately. Sherpa here is a single-threaded process with no Rivet
# subprocess (the analysis runs in-process), so `-c 2` would idle a core and
# double the charge for nothing.
#
# Sherpa event generation, one seed per array task. Submit via
# scripts/submit.sh <CARD> [NSEED], which sets CARD, --array, --time and the
# per-card log paths.
set -euo pipefail

source ./scripts/job_config.sh
source ./scripts/card_sets.sh

CARD=${CARD:?CARD not set - submit with scripts/submit.sh}
card_parse "${CARD}" >/dev/null
SEED=${SLURM_ARRAY_TASK_ID:?no SLURM_ARRAY_TASK_ID}
EVT=${EVT:-$(evt_for "${CARD}")}

TOP=${SLURM_SUBMIT_DIR:-${PWD}}
CARDFILE="${TOP}/cards/generated/${CARD}.yaml"
OUT="${TOP}/yodas/${CARD}"
WORK="${TOP}/work/${CARD}/${SEED}"
RES="${TOP}/int-res/${CARD}"
CPP="${TOP}/proc/${CARD}"
LOG="${TOP}/logs/${CARD}/${SEED}.out"
mkdir -p "${OUT}" "${WORK}"

assert_sherpa_has_ladder_weights

# Refuse to run without grids. Test for the zip as well as the directory:
# Sherpa writes int-res/<CARD>.zip and leaves int-res/<CARD>/ empty, so a plain
# `-d` test passes on a bare mkdir.
if [[ ! -f "${RES}.zip" ]] && [[ -z $(ls -A "${RES}" 2>/dev/null || true) ]]; then
    echo "ERROR: no integration grids for ${CARD} -- run scripts/integrate.sbatch first." >&2
    exit 1
fi

# The grids must come from the binary we are about to run.
check_provenance "${TOP}/int-res/${CARD}.provenance"

OLARG=()
[[ $(card_rung "${CARD}") == L6_nlo ]] && OLARG=( "OL_PREFIX:${OL_PREFIX}" )

# Run from the per-seed work directory. Sherpa writes Settings_Report/ into
# SHERPA_RUN_PATH, which is the cwd, and several hundred concurrent tasks
# sharing one cwd all write the same file. Everything Sherpa needs is passed as
# an absolute path, so a private cwd costs one directory and removes the race.
#   SHERPA_CPP_PATH -> proc/<CARD>, so the Comix process setup is PER CARD.
#     Left at the default it would be one Process/Comix.zip in the shared run
#     directory shared by all 81 cards, which are different final states.
#   RESULT_DIRECTORY -> int-res/<CARD>, the grids integrate.sbatch built.
#   -R ${SEED} is essential: no runcard sets RANDOM_SEED, so without it every
#     array task regenerates the identical sample. -s is accepted and does
#     nothing; check_seed_line below asserts the banner agrees with -R.
#   No -w: EVENT_GENERATION_MODE lives in the card and submit.sh asserts it
#     says Weighted, so there is exactly one place that decides it.
cd "${WORK}"
echo "### ${CARD} seed ${SEED}: ${EVT} weighted events"
t0=${SECONDS}
"${sherpa}" \
    -f "${CARDFILE}" \
    -e "${EVT}" \
    -R "${SEED}" \
    -A "${WORK}/${SEED}" \
    "RESULT_DIRECTORY:${RES}" \
    "SHERPA_CPP_PATH:${CPP}" \
    ${OLARG[@]+"${OLARG[@]}"} \
    BATCH_MODE:2
dt=$(( SECONDS - t0 ))
# Spaces around the ternary are load-bearing: Ares' /usr/bin/awk parses the
# `/t:0` in `t>0?e/t:0` as the start of a regex literal and dies with a syntax
# error, printing nothing. The rate then comes out empty, which is silent here
# and only surfaces much later as smoke.sh's sizing table dividing by zero.
echo "### generation took ${dt} s ($(awk -v e="${EVT}" -v t="${dt}" 'BEGIN{printf "%.1f", t > 0 ? e / t : 0}') ev/s)"

# Log assertions. Each of these was a silent wrong result, not a crash.
check_seed_line "${LOG}" "${SEED}"
check_no_weightname_error "${LOG}"
check_unused_settings "${LOG}"

[[ -s "${WORK}/${SEED}.yoda.gz" ]] || { echo "FAIL no-yoda for seed ${SEED}" >&2; exit 1; }

# Per-seed column check BEFORE publishing. A seed whose YFS.* set is short is
# not merely useless -- rivet-merge would drop the mismatched objects from the
# merge while still normalising by the full sumW, so it would drag the survivors
# down by this seed's share.
python3 "${TOP}/scripts/check-yoda.py" --require "$(required_columns "${CARD}")" \
        "${WORK}/${SEED}.yoda.gz"

# Atomic publish: Rivet rewrites its -o file periodically, so a partial file is
# indistinguishable from a finished one and merge.sbatch's glob would take it.
# Rename is atomic within $SCRATCH.
mv -f "${WORK}/${SEED}.yoda.gz" "${OUT}/${SEED}.yoda.gz" \
    || { echo "FAIL publish-yoda" >&2; exit 1; }

# Work dir only survives a failure, where it is the forensic record (the
# Settings_Report in particular). On success it is 4 files x 400 tasks of inode
# pressure on a filesystem already at 68% of its 1 M inode quota.
cd "${TOP}"
rm -rf "${WORK}"
echo "### OK ${CARD} seed ${SEED}: ${EVT} events -> ${OUT}/${SEED}.yoda.gz"
