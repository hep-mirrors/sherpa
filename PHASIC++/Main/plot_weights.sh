#!/bin/bash

test -z "$1" && exit 1;

rdb=$1;
n=6;
for i in $(unzip -l $rdb '*WD*/2_*' | grep WD | awk '{print $NF}'); do
  if ! test -z "$2"; then
    if ! $(echo $i | awk '{ split($0,a,"/"); if (match(a[3],"'$2'")>0) exit 0; exit 1; }'); then continue; fi
  fi
  (( ++n )); if (( n ==5 )); then (( ++n )); fi
  mkdir -p $(dirname $i);
  unzip -p $rdb $i > $i; sed -e '1 d' -i $i;
  if test -z "$plotcmd"; then plotcmd="plot ";
  else plotcmd=$plotcmd", "; fi;
  t=$(echo $i | sed -e 's|.*WD_._.__.*/||g;s|__QCD\(.*\)| \1|g;s|BVI|S|g;s|j (RS|(H|g;s|[0-9]_[0-9]__||g;s|_| |g');
  plotcmd=$plotcmd"'$i' u (10**\$1):2 t '$t' lc $n lt $n";
done;

gnuplot <<EOF
set key left top
set term postscript color;
set output 'weights_plot.ps';
set logscale xy;
set format y "%2.0te{%L}"
set format x "%2.0te{%L}"
#set xrange [1e-10:1e-3];
#set yrange [1e1:1e4];
set ylabel "# Points"
set xlabel "Weight [1/GeV^2]"
$plotcmd
EOF
