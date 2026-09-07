#! /bin/bash
# Push this campaign to Ares: cards/make-cards.py, budget.py, the Rivet
# analysis source and these scripts. Then generate the cards there.
#
# RSYNC IS THE ONLY TRANSPORT. None of this campaign is in Sherpa's git -- not
# the scripts, not the cards, not budget.py, not the analysis -- so there is no
# checkout on the cluster that could be assumed to already have them. Everything
# the jobs need is copied by this script and nothing else.
#
#   ./scripts/push-to-ares.sh
#
# The 81 run cards are NOT pushed. cards/make-cards.py is deterministic and
# takes no arguments, so the cluster generates its own (scripts/prepare.sh does
# it). Shipping the output instead would let a hand-edited card live on the
# cluster that nothing here reproduces -- a sample with no provenance.
#
# Excludes *.so and *.dylib: the macOS build of the Rivet plugin must never
# shadow the one scripts/build-analysis.sh builds on Ares. Rivet would fail to
# load a Mach-O object, but only when an analysis is looked up -- i.e. after the
# array has queued.
#
# Deliberately does NOT touch Sherpa. The worktree under plggyfsteam is shared
# and is updated and rebuilt by hand; nothing here may race that.
#
# --delete is scoped to the individual directories listed, so nothing this
# script does can reach yodas/, int-res/, proc/, merged/ or logs/.
set -euo pipefail

cd "$(dirname "$0")/.."
: "${ARES_DIR:=/net/afscra/people/plgaprice/yfs-ww}"

MAKECARDS=../cards/make-cards.py
BUDGET=../budget.py
ANALYSIS_SRC=${ANALYSIS_SRC:-$HOME/Documents/research/rivet-analysis/yfs-ww.cc}

for f in "${MAKECARDS}" "${BUDGET}" "${ANALYSIS_SRC}"; do
    [[ -f ${f} ]] || { echo "ERROR: missing ${f}" >&2; exit 1; }
done

echo "### pushing to ares:${ARES_DIR}"
ssh ares "mkdir -p ${ARES_DIR}/cards ${ARES_DIR}/analysis ${ARES_DIR}/scripts ${ARES_DIR}/logs"

rsync -az --exclude='*~' "${MAKECARDS}" ares:"${ARES_DIR}/cards/"
rsync -az --exclude='*~' "${BUDGET}"    ares:"${ARES_DIR}/"
# NO --delete here. The source is a single file, so --delete buys nothing, and
# the destination directory holds things this script must not touch: the
# RivetYFSWW.so built on Ares and build-analysis.sh's lockfile. The *.so exclude
# already stops the macOS Mach-O object being copied up; --delete would be a
# second, load-bearing reason not to remove that exclude, which is exactly the
# kind of hidden coupling that bites later.
rsync -az --exclude='*.so' --exclude='*.dylib' --exclude='*~' \
    "${ANALYSIS_SRC}" ares:"${ARES_DIR}/analysis/"
rsync -az --delete --exclude='*~' --exclude='__pycache__/' \
    scripts/ ares:"${ARES_DIR}/scripts/"

ssh ares "chmod +x ${ARES_DIR}/scripts/*.sh ${ARES_DIR}/scripts/*.py ${ARES_DIR}/budget.py ${ARES_DIR}/cards/make-cards.py"

echo "pushed."
echo "next, in separate short ssh calls (long chains get SIGHUP'd mid-way):"
echo "  ssh ares 'cd ${ARES_DIR} && ./scripts/prepare.sh'"
echo "  ssh ares 'cd ${ARES_DIR} && ./scripts/build-analysis.sh'"
