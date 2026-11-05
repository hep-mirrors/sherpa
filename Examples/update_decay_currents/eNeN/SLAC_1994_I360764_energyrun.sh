#!/usr/bin/env bash
# Runs a runcard once per kinematic point of SLAC_1994_I360764:
#  - data points: each (Q^2, theta) of Tables 1-4, E0 from elastic kinematics
#  - sampling points: each Q^2 of Table 6 at the EPS values below, for the
#    Rosenbluth fits that give the Table 6 form factors
# Each run gets a Q2 selector matching its +-DW theta window, so every event
# lands in its window. Beams are head-on in the e p centre-of-mass frame;
# the analysis boosts back to the proton rest frame.
set -euo pipefail

EVENTS=""             # events per run; empty = the runcard's EVENTS
TAG=""                # output goes to <script dir>/<analysis>_<TAG>; empty = no suffix
DW=0.1                # theta half-window [deg] of sampling points, = _dw in the .cc
EPS="0.1 0.3 0.5"     # eps of sampling points, = _eps in the .cc; empty = data points only

ana=SLAC_1994_I360764
here=$(cd "$(dirname "$0")" && pwd)
anadir=${ANADIR:-$HOME/code_analysis/MyRivetAnalyses}   # where the .cc/.so/.yoda live
card=${1:-$here/Sherpa.ep.yaml}   # runcard from any directory: ./energyrun_...sh path/to/card.yaml
card=$(cd "$(dirname "$card")" && pwd)/$(basename "$card")   # absolute, since Sherpa runs in sub-dirs
out=$here/${ana}${TAG:+_$TAG}

# local plugin and reference data
export RIVET_ANALYSIS_PATH=$anadir${RIVET_ANALYSIS_PATH:+:$RIVET_ANALYSIS_PATH}
export RIVET_DATA_PATH=$anadir${RIVET_DATA_PATH:+:$RIVET_DATA_PATH}

# one line per run: tag, lab E0, CM E_e, CM E_p, Q^2 range of the theta window [GeV, GeV^2]
points=$(awk -v DW="$DW" -v EPS="$EPS" '
  function rad(t) {return t*atan2(0,-1)/180}
  function e0(Q, t,  s) {s=sin(rad(t)/2)^2; return Q/(4*M) + sqrt((Q/(4*M))^2 + Q/(4*s))}
  function q2(E, t,  c) {c=1-cos(rad(t)); return 2*E*E*c/(1+E/M*c)}
  function pt(tag, Q, tlo, thi,  E, s, rs) {E=e0(Q,(tlo+thi)/2); s=M*M+me*me+2*E*M; rs=sqrt(s)
    printf "%s %.10g %.10g %.10g %.6g %.6g\n", tag, E, (s+me*me-M*M)/(2*rs), (s+M*M-me*me)/(2*rs), q2(E,tlo), q2(E,thi)}
  BEGIN {M=0.938272; me=0.000511}
  /^BEGIN/ {split($3,p,"/"); h=p[4]}
  /^Edges\(A1\)/ {sub(/.*\[/,""); sub(/\].*/,""); n[h]=split($0,x,","); for (i=1;i<=n[h];++i) e[h,i]=x[i]+0}
  END {
    for (j=1;j<=4;++j) {h=sprintf("d%02d-x01-y01",j); Q=(e["d05-x01-y01",2*j-1]+e["d05-x01-y01",2*j])/2
      for (i=1;i<n[h];i+=2) pt(sprintf("d%02d_th%.4g",j,(e[h,i]+e[h,i+1])/2), Q, e[h,i], e[h,i+1])}
    m=split(EPS,eps," "); h="d06-x01-y01"
    for (i=1;i<n[h];i+=2) {Q=(e[h,i]+e[h,i+1])/2; tau=Q/(4*M*M)
      for (k=1;k<=m;++k) {t=2*atan2(sqrt((1/eps[k]-1)/(2*(1+tau))),1)*180/atan2(0,-1)
        pt(sprintf("q2_%.4g_eps%s",Q,eps[k]), Q, t-DW, t+DW)}}
  }' "$anadir/$ana.yoda")

[ -n "$EVENTS" ] || EVENTS=$(awk '/^EVENTS:/{print $2}' "$card")

while read tag en ee ep q2lo q2hi; do
  dir=$out/$tag
  mkdir -p "$dir"
  sed -e "s/^EVENTS:.*/EVENTS: $EVENTS/" \
      -e "s/^BEAMS:.*/BEAMS: [11, 2212]/" \
      -e "s/^BEAM_ENERGIES:.*/BEAM_ENERGIES: [$ee, $ep]/" \
      -e "s/^- \[Q2,.*/- [Q2, 11, 11, $q2lo, $q2hi]/" \
      -e "s/^ANALYSIS_OUTPUT:.*/ANALYSIS_OUTPUT: ${ana}_$tag/" \
      -e "s/^  --analyses:.*/  --analyses: [MC_XS, $ana]/" \
      "$card" > "$dir/Sherpa.yaml"
  echo "$tag: E0 = $en GeV, $q2lo < Q2 < $q2hi GeV^2"
  (cd "$dir" && ${SHERPA:-Sherpa} > sherpa.log 2>&1) || echo "  failed, see $dir/sherpa.log"
done <<< "$points"

# Each run fills only its own window, so merge them as separate
# processes (no -e). The form factors are refitted in finalize.
rivet-merge -o "$out/$ana.yoda.gz" "$out"/*/${ana}_*.yoda.gz
