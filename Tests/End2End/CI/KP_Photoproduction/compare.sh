#!/bin/sh
# usage: compare.sh <log with process mapping> <log without process mapping>
# Every integrated cross section must be finite and agree between the two runs
# within five combined standard deviations.
set -eu

# process, cross section and error of the last integration of each process
# (the log file is appended to on repeated runs)
xs() {
  sed -n 's/^\([^ ]*\) : \([^ ]*\) pb +- ( \([^ ]*\) pb.*/\1 \2 \3/p' "$1" |
    awk '{ last[$1] = $0 } END { for (p in last) print last[p] }' | sort
}

xs "$1" > mapped.xs
xs "$2" > unmapped.xs
[ -s mapped.xs ] || { echo "no cross sections in $1"; exit 1; }

join mapped.xs unmapped.xs | awk -v n="$(wc -l < mapped.xs)" '
  { a = $2 + 0; b = $4 + 0; s = sqrt($3 * $3 + $5 * $5)
    ok = ($2 ~ /^[-+0-9.eE]+$/ && $4 ~ /^[-+0-9.eE]+$/ && (a - b) ^ 2 <= 25 * s * s)
    printf "%s %s: %s +- %s vs %s +- %s pb\n", (ok ? "ok  " : "FAIL"), $1, $2, $3, $4, $5
    if (!ok) bad = 1; m++ }
  END { if (m != n || m == 0) { print "process lists differ"; bad = 1 }
        exit bad }'
