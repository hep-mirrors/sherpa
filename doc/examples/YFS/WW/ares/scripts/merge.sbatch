#!/bin/bash
#SBATCH -J yfsww-merge
#SBATCH -p plgrid
#SBATCH -A plgsherpalumi-cpu
#SBATCH -t 02:00:00
#SBATCH -c 1
#SBATCH -x ac[0770-0780],ac0714
#
# rivet-merge, never yodamerge: yfsww::finalize() scales every histogram by
# crossSection()/sumOfWeights(), so averaging the FINALISED objects is wrong.
# rivet-merge sums the /RAW/ objects and re-runs finalize(), which is why it
# needs RIVET_ANALYSIS_PATH (job_config.sh sets it).
#
#   sbatch --export=ALL,CARD=L6_nlo_161p0_mue scripts/merge.sbatch
#
# Output goes to one directory per (rung, energy, channel), which is the layout
# budget.py reads:  merged/<energy>/<channel>/<rung>/<CARD>.yoda.gz
# The m_W probes go to merged/<energy>/_mw/<mass>/ and are linked in as
# <energy>/<channel>/mw by scripts/link-mw.sh.
set -euo pipefail

source ./scripts/job_config.sh
source ./scripts/card_sets.sh

CARD=${CARD:?CARD not set - submit with scripts/submit.sh}
card_parse "${CARD}" >/dev/null
TOP=${SLURM_SUBMIT_DIR:-${PWD}}

SRC="${TOP}/yodas/${CARD}"
DEST="${TOP}/$(merged_dir "${CARD}")"
mkdir -p "${DEST}"

# `|| true`: with `set -o pipefail`, a failing glob makes the whole pipeline
# non-zero and `set -e` kills the script with no output at all.
files=( $(ls -1 "${SRC}"/*.yoda.gz 2>/dev/null || true) )
if [[ ${#files[@]} -eq 0 ]]; then
    echo "ERROR: no yoda files in ${SRC}" >&2
    exit 1
fi
echo "### merging ${#files[@]} seed(s) for ${CARD}"

# ONE CARD ONLY. Never glob across cards: yfsww books photon-E-log and every
# energy-scaled histogram against sqrtS(), so the three scan points have
# DIFFERENT bin edges. rivet-merge silently drops objects whose binning does
# not match across inputs while still normalising by the total sumW, so a
# cross-energy merge would leave the survivors low by the fraction of seeds
# from the other energy -- with nothing in the output to say so.
#
# No -j: this Rivet's rivet-merge is out of step with its compiled bindings and
# the parallel path dies in getRawAOpaths, but only when len(files) > 2*jobs,
# so a small test merge passes and the real one does not.
#
# --ignore-beams because the cards set RIVET: IGNORE_BEAMS: 1, so the beam
# records are not consistent across files and -e otherwise aborts with
#   Equivalent merging requires matching beams across input files!
OUT="${DEST}/${CARD}.yoda.gz"
rivet-merge -e --assume-reentrant --ignore-beams -o "${OUT}" "${files[@]}"

# The merged file gets the same assertions as each seed did, plus the one that
# only a merge can fail: acceptance must still be exactly 1. A merge that
# dropped mismatched objects while keeping the full sumW normalisation shows up
# here as an acceptance below 1 by the dropped fraction.
python3 "${TOP}/scripts/check-yoda.py" --require "$(required_columns "${CARD}")" "${OUT}"

echo "### wrote ${OUT} from ${#files[@]} seed(s)"
