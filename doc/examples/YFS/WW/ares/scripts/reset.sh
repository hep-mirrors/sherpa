#! /bin/bash
# Clear everything belonging to a card SET so the next submission starts from
# nothing instead of topping up what is there.
#
#   ./scripts/reset.sh l6                    # DRY RUN: list what would go
#   ./scripts/reset.sh l6 --yes              # actually delete
#   ./scripts/reset.sh l6 --yes --grids      # also drop int-res/ and proc/
#   ./scripts/reset.sh L6_nlo_161p0_mue --yes
#
# WHY this exists rather than `rm -rf yodas/*`: merge.sbatch globs a card's whole
# directory, so leftovers from an earlier generator or analysis version are not
# an error -- they are silently averaged into the new merge. After any change to
# the YFS code, the cards or the Rivet analysis that is exactly wrong, and
# nothing in the output says so.
#
# --grids is separate because a clean EVENT run does not need a re-integration,
# and re-integrating 81 cards is not free. Use it after a change to the matrix
# element or the process definition; skip it after a change to the Rivet
# analysis or to the number of events.
#
# THE CARD LIST COMES FROM card_sets.sh, the same file submit-all.sh and
# integrate-all.sh use. A reset carrying its own stale copy does not error --
# it leaves one family's yodas in place and the next merge averages two physics
# versions.
#
# DRY RUN IS THE DEFAULT, and deletion is not reversible here: `rip` is a
# local-only tool, it does not exist on Ares, and `rip ... 2>/dev/null` there is
# a silent no-op that looks like it worked. This uses plain rm.
set -uo pipefail

cd "$(dirname "$0")/.."
# NOT `|| true`: without job_config.sh there is no merged_dir(), and reset
# would then quietly leave every merged yoda in place.
source ./scripts/job_config.sh >/dev/null
source ./scripts/card_sets.sh

YES=0; GRIDS=0
positional=()
for a in "$@"; do
    case "${a}" in
        --yes|-y)  YES=1 ;;
        --grids)   GRIDS=1 ;;
        --help|-h) sed -n '2,30p' "$0"; exit 0 ;;
        -*) echo "ERROR: unknown flag '${a}'" >&2; exit 1 ;;
        *)  positional+=( "${a}" ) ;;
    esac
done
SET=${positional[0]:-}
[[ -n ${SET} ]] || { echo "usage: reset.sh <SET|CARD> [--yes] [--grids]" >&2
                     echo "       aliases: $(set_names)" >&2; exit 1; }
resolve_set "${SET}" || exit 1

# Refuse while anything is queued. A task that publishes its yoda after the
# directory has been cleared leaves a single-seed sample that looks complete,
# and its chained merge fires on afterany regardless -- merged/ would then hold
# one seed's statistics with nothing to indicate it.
queued=$(squeue -u "${USER}" -h -t pending,running -r 2>/dev/null | wc -l | tr -d ' ')
if [[ ${queued} -gt 0 ]]; then
    echo "WARNING: ${queued} job(s) still queued or running for ${USER}." >&2
    echo "         Clearing now races them. Cancel (scancel -u ${USER}) or wait." >&2
    if [[ ${YES} -eq 1 ]]; then
        echo "         REFUSING to delete with jobs in flight." >&2
        exit 1
    fi
fi

targets=()
add() { [[ -e $1 || -L $1 ]] && targets+=( "$1" ); return 0; }
for c in "${cards[@]}"; do
    add "yodas/${c}"
    add "logs/${c}"
    add "work/${c}"
    add "$(merged_dir "${c}")/${c}.yoda.gz"
    if [[ ${GRIDS} -eq 1 ]]; then
        add "int-res/${c}"
        add "int-res/${c}.zip"
        add "proc/${c}"
    fi
done

echo "set '${SET}': ${#cards[@]} card(s)"
if [[ ${GRIDS} -eq 1 ]]; then
    echo "  --grids: int-res/ and proc/ WILL be removed (re-integration needed)"
else
    echo "  integration grids and process setup kept (--grids drops them too)"
fi
echo

if [[ ${#targets[@]} -eq 0 ]]; then
    echo "nothing to remove -- already clean."
    exit 0
fi

total=0
for t in "${targets[@]}"; do
    if [[ -d ${t} ]]; then
        n=$(find "${t}" -type f 2>/dev/null | wc -l | tr -d ' ')
        printf "  %-56s dir   %6s files  %6s\n" "${t}" "${n}" "$(du -sh "${t}" 2>/dev/null | cut -f1)"
        total=$(( total + n ))
    else
        printf "  %-56s file                 %6s\n" "${t}" "$(du -h "${t}" 2>/dev/null | cut -f1)"
        total=$(( total + 1 ))
    fi
done
echo
echo "  ${#targets[@]} path(s), ${total} file(s)"

if [[ ${YES} -ne 1 ]]; then
    echo
    echo "DRY RUN -- nothing deleted. Re-run with --yes to go ahead."
    exit 0
fi

echo
for t in "${targets[@]}"; do rm -rf -- "${t}" || echo "  FAILED to remove ${t}" >&2; done

left=0
for t in "${targets[@]}"; do [[ -e ${t} ]] && left=$(( left + 1 )); done
if [[ ${left} -gt 0 ]]; then
    echo "ERROR: ${left} path(s) survived deletion -- check permissions." >&2
    exit 1
fi
echo "removed ${#targets[@]} path(s). '${SET}' is clean."
