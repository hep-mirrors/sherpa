#! /bin/bash
# Run budget.py over every (energy, channel) that has merged output.
#
#   ./scripts/budget-all.sh                 # on Ares, over merged/
#   ./scripts/budget-all.sh results         # locally, after fetch-from-ares.sh
#
# Read the cross-section table before reading any shape. budget.py marks any
# rung whose own error exceeds 5% with (!) -- a rung whose error is comparable
# to the shift it is meant to measure is not a measurement, and the shapes from
# that sample are not either.
set -euo pipefail

cd "$(dirname "$0")/.."
ROOT=${1:-merged}
[[ -d ${ROOT} ]] || { echo "ERROR: no such directory: ${ROOT}" >&2; exit 1; }

found=0
for d in "${ROOT}"/*/*/; do
    [[ -d ${d} ]] || continue
    case "${d}" in */_mw/) continue ;; esac
    python3 budget.py "${d%/}"
    found=1
done
[[ ${found} -eq 1 ]] || echo "nothing merged under ${ROOT} yet"
