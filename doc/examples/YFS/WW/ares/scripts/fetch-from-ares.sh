#! /bin/bash
# Pull the merged yodas back from Ares for local plotting and budget.py.
#
#   ./scripts/fetch-from-ares.sh                 # merged/ only
#   ./scripts/fetch-from-ares.sh --logs          # also the job logs
#
# Lands in results/ next to this directory, laid out exactly as budget.py wants:
#
#   results/<energy>/<channel>/{L0_born,L1_isr,L4_master,L6_nlo}/<CARD>.yoda.gz
#   results/<energy>/<channel>/mw -> ../_mw        (three m_W points)
#
# so   python3 ../budget.py results/161p0/mue   is the whole analysis step.
#
# NOT --delete: a fetch must never be able to remove local results. Re-running
# it after more seeds have merged simply refreshes the files that changed.
set -euo pipefail

cd "$(dirname "$0")/.."
: "${ARES_DIR:=/net/afscra/people/plgaprice/yfs-ww}"

LOGS=0
for a in "$@"; do case "${a}" in --logs) LOGS=1 ;; *) echo "ERROR: unknown flag '${a}'" >&2; exit 1 ;; esac; done

mkdir -p results
# -l keeps the mw symlinks as symlinks; they are relative, so they resolve here
# exactly as they do on Ares.
rsync -azl --progress ares:"${ARES_DIR}/merged/" results/

if [[ ${LOGS} -eq 1 ]]; then
    mkdir -p results-logs
    rsync -az --progress ares:"${ARES_DIR}/logs/" results-logs/
fi

# Belt and braces: if the merges on Ares ran before link-mw.sh did, the symlinks
# are not there to copy. Recreate them locally from the same rule.
if [[ -x ./scripts/link-mw.sh ]]; then
    ( cd results/.. >/dev/null
      MERGEDROOT=results ENERGIES_ALL="157p5 161p0 162p5" \
      CHANNELS_ALL="mue emu mutau taumu etau taue" \
      bash -c '
        for e in ${ENERGIES_ALL}; do
          [[ -d ${MERGEDROOT}/${e}/_mw ]] || continue
          for ch in ${CHANNELS_ALL}; do
            d=${MERGEDROOT}/${e}/${ch}
            [[ -d ${d} ]] || continue
            [[ -e ${d}/mw && ! -L ${d}/mw ]] && continue
            ln -sfn ../_mw ${d}/mw
          done
        done' )
fi

echo
echo "fetched into $(pwd)/results"
found=0
for d in results/*/*/; do
    [[ -d ${d} ]] || continue
    case "${d}" in */_mw/) continue ;; esac
    echo "  python3 ../budget.py ${d%/}"
    found=1
done
[[ ${found} -eq 1 ]] || echo "  (nothing merged yet)"
