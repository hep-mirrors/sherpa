#!/usr/bin/env bash
# Runs Sherpa.ep.yaml once for each beam energy of JEFFERSON_2004_I643262
# (one HEPData table per energy). The beams are set head-on in the e p
# centre-of-mass frame with the same sqrt(s) as the fixed-target beam; the
# analysis boosts back to the proton rest frame and picks the table from
# sqrt(s). Each energy gets its own sub-directory.
set -euo pipefail

EVENTS=""   # events per energy; empty = the runcard's EVENTS
TAG=""      # output goes to <script dir>/<analysis>_<TAG>; empty = no suffix

ana=JEFFERSON_2004_I643262
here=$(cd "$(dirname "$0")" && pwd)
anadir=${ANADIR:-$HOME/code_analysis/MyRivetAnalyses}   # where the .cc/.so/.yoda live
card=${1:-$here/Sherpa.ep.yaml}   # runcard from any directory: ./energyrun_...sh path/to/card.yaml
card=$(cd "$(dirname "$card")" && pwd)/$(basename "$card")   # absolute, since Sherpa runs in sub-dirs
out=$here/${ana}${TAG:+_$TAG}

# local plugin and reference data
export RIVET_ANALYSIS_PATH=$anadir${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}
export RIVET_DATA_PATH=$anadir${RIVET_DATA_PATH:+:$RIVET_DATA_PATH}

# lab beam energies [GeV] of Tables 1-7 (same as _ebeam in the .cc)
energies="1.148 1.882 2.235 3.114 4.104 4.413 5.494"

[ -n "$EVENTS" ] || EVENTS=$(awk '/^EVENTS:/{print $2}' "$card")

n=0
for en in $energies; do
  n=$((n+1))
  # e and p energies in the centre-of-mass frame, s = mp^2 + me^2 + 2 E mp
  read ee ep < <(awk -v E="$en" 'BEGIN{mp=0.938272; me=0.000511; s=mp*mp+me*me+2*E*mp; rs=sqrt(s);
                                      printf "%.10g %.10g\n", (s+me*me-mp*mp)/(2*rs), (s+mp*mp-me*me)/(2*rs)}')
  # angle windows of table n -> centre and elastic Q^2 range [GeV^2], so every event lands in its window
  windows=$(awk -v h="/REF/$ana/d0$n-x01-y01" -v E="$en" '
    function q2(t) {c = 1 - cos(t*atan2(0,-1)/180); return 2*E*E*c/(1 + E/0.938272*c)}
    /^BEGIN/ {f=($3==h)}
    f && /^Edges\(A1\)/ {sub(/.*\[/,""); sub(/\].*/,""); m=split($0,e,",");
                         for (i=1;i<m;i+=2) printf "%.4g %.6g %.6g\n", (e[i]+e[i+1])/2, q2(e[i]), q2(e[i+1])}' "$anadir/$ana.yoda")
  while read th q2lo q2hi; do
    dir=$out/$en/$th
    mkdir -p "$dir"
    sed -e "s/^EVENTS:.*/EVENTS: $EVENTS/" \
        -e "s/^BEAMS:.*/BEAMS: [11, 2212]/" \
        -e "s/^BEAM_ENERGIES:.*/BEAM_ENERGIES: [$ee, $ep]/" \
        -e "s/^- \[Q2,.*/- [Q2, 11, 11, $q2lo, $q2hi]/" \
        -e "s/^ANALYSIS_OUTPUT:.*/ANALYSIS_OUTPUT: ${ana}_${en}_$th/" \
        -e "s/^  --analyses:.*/  --analyses: [MC_XS, $ana]/" \
        "$card" > "$dir/Sherpa.yaml"
    echo "E_beam = $en GeV, theta = $th deg, $q2lo < Q2 < $q2hi GeV^2"
    (cd "$dir" && ${SHERPA:-Sherpa} > sherpa.log 2>&1) || echo "  failed, see $dir/sherpa.log"
  done <<< "$windows"
done

# Each run fills only its own angle window, so merge them as separate
# processes (no -e).
rivet-merge -o "$out/$ana.yoda.gz" "$out"/*/*/${ana}_*.yoda.gz
