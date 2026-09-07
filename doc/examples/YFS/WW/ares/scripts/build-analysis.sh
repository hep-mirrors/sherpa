#! /bin/bash
# Build the yfs-ww Rivet plugin ON ARES.
#
#   ssh ares 'cd /net/afscra/people/plgaprice/yfs-ww && ./scripts/build-analysis.sh'
#
# The macOS build is a Mach-O object and useless here; push-to-ares.sh does not
# copy it, and this refuses to run if one turned up anyway. A stale plugin does
# not announce itself -- Rivet simply does not find the analysis, and the run
# fails after the array has already queued.
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh

cd analysis
[[ -f yfs-ww.cc ]] || { echo "ERROR: no yfs-ww.cc here -- run push-to-ares.sh" >&2; exit 1; }

# Refuse to run twice at once. Two concurrent builds share one output file: the
# second one's `rm -f` deletes the .so the first is still validating against,
# and the first then reports every card as naming an unknown analysis. Observed,
# not hypothetical. Both runs can still exit having printed a verdict, and the
# verdict of the loser is garbage.
exec 9>.build-analysis.lock
if ! flock -n 9; then
    echo "ERROR: another build-analysis.sh is running in this directory." >&2
    echo "       Two builds share one RivetYFSWW.so and the loser's checks are" >&2
    echo "       meaningless. Wait for it to finish." >&2
    exit 1
fi
if file ./*.so 2>/dev/null | grep -qi 'Mach-O'; then
    echo "ERROR: a Mach-O .so is present. Remove it; it will shadow this build." >&2
    exit 1
fi

rm -f RivetYFSWW.so
rivet-build RivetYFSWW.so yfs-ww.cc

[[ -f RivetYFSWW.so ]] || { echo "ERROR: rivet-build produced nothing" >&2; exit 1; }
[[ RivetYFSWW.so -nt yfs-ww.cc ]] || {
    echo "ERROR: RivetYFSWW.so is not newer than yfs-ww.cc -- the build did not run" >&2
    exit 1; }

echo "### ldd"
ldd RivetYFSWW.so | grep -iE 'rivet|yoda|hepmc|not found' || true
if ldd RivetYFSWW.so | grep -q 'not found'; then
    echo "ERROR: unresolved libraries in RivetYFSWW.so" >&2
    exit 1
fi

# job_config.sh points RIVET_ANALYSIS_PATH at this directory, so a successful
# build must make every analysis the runcards name visible. Anything else and
# the runs die at initialisation, one array task at a time.
#
# The list is taken ONCE. `rivet --list-analyses` loads every plugin in the
# path and takes seconds; calling it per card turned this script into a
# multi-minute job for no reason.
# awk '{print $1}' is load-bearing: `rivet --list-analyses` pads each name with
# trailing spaces, so a `grep -x yfs-ww` straight off it never matches and the
# check would report a perfectly good plugin as missing.
known=$(rivet --list-analyses 2>/dev/null | awk '{print $1}')
have() { echo "${known}" | grep -qx "$1"; }

for a in ${RIVET_ANALYSES}; do
    name=${a%%:*}
    have "${name}" || { echo "ERROR: '${name}' not visible to Rivet after the build" >&2; exit 1; }
    echo "ok: ${name}"
done

# Every card must ask for an analysis that now exists. A card naming one that is
# not here fails at initialisation, one array task at a time, after the
# integration has already been paid for.
missing=0
n=0
for f in ../cards/generated/*.yaml; do
    n=$(( n + 1 ))
    for a in $(awk '/^ *ANALYSES:/{g=1;next} g&&/^ *- /{print $2;next} g&&!/^ *- /{g=0}' "${f}"); do
        have "${a}" || { echo "ERROR: $(basename "${f}") asks for unknown analysis '${a}'" >&2; missing=1; }
    done
done
[[ ${missing} -eq 0 ]] || exit 1
echo "all ${n} cards name only analyses Rivet can see"
