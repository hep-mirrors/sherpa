#!/usr/bin/env bash
# Runs Sherpa.eepipipi.yaml at the centre of every sqrt(s) bin of
# BESIII_2019_I1773081. The analysis fills sqrtS() into its reference
# binning and has no beam-energy check, so any point inside a bin works.
# Each energy gets its own sub-directory, so Results/Process from one
# energy are never reused for another.
set -euo pipefail

STRIDE=10   # run only every STRIDE-th energy bin (1 = all of them)
EVENTS=""   # events per energy; empty = the runcard's EVENTS
TAG=""      # output goes to <script dir>/<analysis>_<TAG>; empty = no suffix

ana=BESIII_2019_I1773081
hists="d01-x01-y01"   # reference histograms the analysis books
here=$(cd "$(dirname "$0")" && pwd)
card=$here/Sherpa.eepipipi.yaml
out=$here/${ana}${TAG:+_$TAG}
ref=$(rivet-config --datadir)/$ana.yoda.gz

# bin centres [GeV] of the booked reference histograms
energies=$(for h in $hists; do
  gzip -dc "$ref" | awk -v h="/REF/$ana/$h" '
    /^BEGIN/ {f=($3==h)}
    f && /^Edges\(A1\)/ {sub(/.*\[/,""); sub(/\].*/,""); n=split($0,e,",");
                         for (i=1;i<n;++i) printf "%.6g\n", (e[i]+e[i+1])/2}'
done)

energies=$(printf '%s\n' $energies | awk -v s="$STRIDE" '(NR-1)%s==0')
[ -n "$EVENTS" ] || EVENTS=$(awk '/^EVENTS:/{print $2}' "$card")

for rs in $energies; do
  eb=$(awk -v r="$rs" 'BEGIN{printf "%.10g", r/2}')
  dir=$out/$rs
  mkdir -p "$dir"
  sed -e "s/^EVENTS:.*/EVENTS: $EVENTS/" \
      -e "s/^BEAM_ENERGIES:.*/BEAM_ENERGIES: [$eb, $eb]/" \
      -e "s/^ANALYSIS_OUTPUT:.*/ANALYSIS_OUTPUT: ${ana}_$rs/" \
      -e "s/^  --analyses:.*/  --analyses: [MC_XS, $ana]/" \
      "$card" > "$dir/Sherpa.yaml"
  echo "sqrt(s) = $rs GeV, E_beam = $eb GeV"
  (cd "$dir" && Sherpa > sherpa.log 2>&1) || echo "  failed, see $dir/sherpa.log"
done

# Different energies fill different bins, so merge them as separate
# processes (no -e).
rivet-merge -o "$out/$ana.yoda.gz" "$out"/*/${ana}_*.yoda.gz
