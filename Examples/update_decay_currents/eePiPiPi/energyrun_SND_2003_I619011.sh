#!/usr/bin/env bash
# Runs Sherpa.eepipipi.yaml at every sqrt(s) point of SND_2003_I619011.
# Each energy gets its own sub-directory, so Results/Process from one
# energy are never reused for another.
set -euo pipefail

STRIDE=10   # run only every STRIDE-th energy bin (1 = all of them)
EVENTS=""   # events per energy; empty = the runcard's EVENTS
TAG=""      # output goes to <script dir>/<analysis>_<TAG>; empty = no suffix

ana=SND_2003_I619011
here=$(cd "$(dirname "$0")" && pwd)
card=$here/Sherpa.eepipipi.yaml
out=$here/${ana}${TAG:+_$TAG}
info=$(rivet-config --datadir)/$ana.info

# sqrt(s) [GeV] from the Energies: list of the analysis' .info file
energies=$(awk '/^Energies:/{f=1} f{print} f&&/\]/{exit}' "$info" |
           sed -e 's/Energies://' -e 's/[][,]/ /g')

energies=$(printf '%s\n' $energies | awk -v s="$STRIDE" '(NR-1)%s==0')
[ -n "$EVENTS" ] || EVENTS=$(awk '/^EVENTS:/{print $2}' "$card")

for rs in $energies; do
  eb=$(awk -v r="$rs" 'BEGIN{printf "%.10g", r/2}')
  dir=$out/$rs
  mkdir -p "$dir"
  sed -e "s/^EVENTS:.*/EVENTS: $EVENTS/" \
      -e "s/^BEAM_ENERGIES:.*/BEAM_ENERGIES: [$eb, $eb]/" \
      -e "s/^ANALYSIS_OUTPUT:.*/ANALYSIS_OUTPUT: ${ana}_$rs/" \
      "$card" > "$dir/Sherpa.yaml"
  echo "sqrt(s) = $rs GeV, E_beam = $eb GeV"
  (cd "$dir" && Sherpa > sherpa.log 2>&1) || echo "  failed, see $dir/sherpa.log"
done

# Different energies fill different bins, so merge them as separate
# processes (no -e).
rivet-merge -o "$out/$ana.yoda.gz" "$out"/*/${ana}_*.yoda.gz
