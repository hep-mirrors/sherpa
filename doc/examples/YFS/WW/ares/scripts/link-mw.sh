#! /bin/bash
# Make the m_W scan visible to budget.py from every channel directory.
#
# budget.py wants <rundir>/mw/<mass> with <rundir> = merged/<energy>/<channel>,
# but the m_W probes are mue-only (make-cards.py generates nine MW cards, three
# masses at each of three energies, all in the mue channel). merge.sbatch puts
# them in merged/<energy>/_mw/<mass>, and this links that in as
# merged/<energy>/<channel>/mw for all six channels.
#
# Reusing one scan across channels is legitimate: budget.py only uses it for
# dln(sigma)/dm_W, which is a logarithmic derivative, so the leptonic branching
# fractions that distinguish the channels cancel out of it.
#
# Relative symlinks, so the tree survives rsync to the local machine unchanged.
#
#   ./scripts/link-mw.sh            (safe to re-run; run it locally too)
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh

n=0
for e in ${ENERGIES_ALL}; do
    src="merged/${e}/_mw"
    [[ -d ${src} ]] || { echo "skip ${e}: no ${src} yet"; continue; }
    for ch in ${CHANNELS_ALL}; do
        dst="merged/${e}/${ch}"
        [[ -d ${dst} ]] || continue
        # Replace only a symlink. A real directory called mw/ here would be
        # someone's data and silently clobbering it is exactly the kind of
        # thing that turns up three weeks later as an unexplained number.
        if [[ -e ${dst}/mw && ! -L ${dst}/mw ]]; then
            echo "ERROR: ${dst}/mw exists and is not a symlink -- refusing" >&2
            exit 1
        fi
        ln -sfn ../_mw "${dst}/mw"
        n=$(( n + 1 ))
    done
done
echo "linked ${n} channel director$([[ ${n} == 1 ]] && echo y || echo ies) to their _mw scan"
