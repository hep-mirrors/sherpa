#!/bin/bash

test -z "$1" && exit 1;

rdb=$1;
n=0;
for i in $(sqlite3 $rdb "Select file from path where file like '%Statistics.dat'"); do
  (( ++n ));
  mkdir -p $(dirname $i);
  sqlite3 $rdb "select content from path where file = '"$i"'" > $i;
  if test -z "$plotcmd"; then plotcmd="plot ";
  else plotcmd=$plotcmd", "; fi;
  t=$(echo $i | sed -e 's|.*MC_._.__||g;s|/Statistics.dat||g;s|__QCD.*||g');
  fitcmd=$fitcmd"f$n(x) = a$n; fit f$n(x) '$i' using 5:(abs(\$2)):(2*\$3) via a$n;"
  plotcmd=$plotcmd"'$i' u 5:(abs(\$2)):(2*\$3) w yerr t '$t' lc $n lt $n"
  plotcmd=$plotcmd", f$n(x) t sprintf(\"%3.4g pb\",a$n) lc $n lt $n";
done;

gnuplot <<EOF
set term postscript color;
set output 'stats_plot.ps';
set logscale xy;
set xlabel "number of points"
set ylabel "cross section [pb]"
$fitcmd
$plotcmd
EOF
