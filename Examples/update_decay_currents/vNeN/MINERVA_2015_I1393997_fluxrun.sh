#!/usr/bin/env bash
# Runs each given runcard (default: the nu_e n and nu_ebar p cards next to this
# script) at the centre of every flux bin of MINERVA_2015_I1393997 (d01-x01-y01).
# Each card keeps its own BEAMS and PROCESSES, so nu and nubar need separate
# cards; only the neutrino energy, EVENTS, output name and analyses are
# changed. The analysis weights each event by the nu_e or nu_ebar flux
# of its bin, so the runs are summed with rivet-merge.
#   ./fluxrun.sh [card.yaml ...]
set -euo pipefail

EVENTS=""     # events per run; empty = each runcard's EVENTS
TAG=""        # output goes to <script dir>/<analysis>_<TAG>; empty = no suffix

ana=MINERVA_2015_I1393997
here=$(cd "$(dirname "$0")" && pwd)
anadir=${ANADIR:-$HOME/code_analysis/MyRivetAnalyses}   # where the .cc/.so/.yoda live
[ $# -gt 0 ] || set -- "$here/Sherpa.vNeN.yaml" "$here/Sherpa.vNeN_anti.yaml"
out=$here/${ana}${TAG:+_$TAG}
ref=$anadir/$ana.yoda

# local plugin and reference data
export RIVET_ANALYSIS_PATH=$anadir${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}
export RIVET_DATA_PATH=$anadir${RIVET_DATA_PATH:+:$RIVET_DATA_PATH}

# flux bin centres [GeV]
energies=$(awk -v h="/REF/$ana/d01-x01-y01" '
  /^BEGIN/ {f=($3==h)}
  f && /^Edges\(A1\)/ {sub(/.*\[/,""); sub(/\].*/,""); n=split($0,e,",");
                       for (i=1;i<n;++i) printf "%.6g\n", (e[i]+e[i+1])/2}' "$ref")

files=()   # outputs of this call only, so old runs in $out are never merged
for card in "$@"; do
  card=$(cd "$(dirname "$card")" && pwd)/$(basename "$card")   # absolute, since Sherpa runs in sub-dirs
  name=$(basename "$card" .yaml)
  nev=${EVENTS:-$(awk '/^EVENTS:/{print $2}' "$card")}
  for en in $energies; do
    dir=$out/$name/$en
    mkdir -p "$dir"
    sed -e "s/^EVENTS:.*/EVENTS: $nev/" \
        -e "s/^BEAM_ENERGIES: *\[[^,]*,/BEAM_ENERGIES: [$en,/" \
        -e "s/^ANALYSIS_OUTPUT:.*/ANALYSIS_OUTPUT: ${ana}_${name}_$en/" \
        -e "s/^  --analyses:.*/  --analyses: [MC_XS, $ana]/" \
        "$card" > "$dir/Sherpa.yaml"
    echo "$name, E_nu = $en GeV"
    (cd "$dir" && Sherpa > sherpa.log 2>&1) || echo "  failed, see $dir/sherpa.log"
    [ -f "$dir/${ana}_${name}_$en.yoda.gz" ] && files+=("$dir/${ana}_${name}_$en.yoda.gz")
  done
done

# Each run carries its own flux weight, so sum them as separate
# processes (no -e).
rivet-merge -o "$out/$ana.yoda.gz" "${files[@]}"
