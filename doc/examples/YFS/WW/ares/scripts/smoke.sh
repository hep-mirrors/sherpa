#! /bin/bash
# Smoke test: one short task per rung at 161.0 GeV, mue -- the reference point.
#
#   ./scripts/smoke.sh integrate      # 4 integration jobs
#   ./scripts/smoke.sh run            # 4 one-seed arrays + their merges
#   ./scripts/smoke.sh report         # verdict + measured rates + sizing
#   ./scripts/smoke.sh clean          # clear the four cards for production
#
# Three separate calls on purpose. A long ssh chain gets SIGHUP'd part way
# through and the half that ran looks like the half that did not; submit in one
# short call, poll in another.
#
# THIS RUNS THE PRODUCTION PATH, not a parallel copy of it: the same
# integrate.sbatch, the same sherpa_array.sbatch, the same merge.sbatch, only
# with EVT turned down. A smoke test that exercises different code proves
# nothing about the thing that will run for two days. The consequence is that
# it leaves four production cards populated with a tiny sample --
# `smoke.sh clean` removes them, and submit-all.sh refuses to run over them
# until it has been.
#
# What it is checking, all of which fail silently otherwise:
#   * Seed:            the banner agrees with -R, per task and across tasks
#   * no weight name   'GenCrossSection::set_xsec: no weight with given name',
#                      which used to kill the NLO cards at ~486 events
#   * unused settings  the whole YFS block being dropped by a binary that
#                      predates Ladder_Weights
#   * YFS.* columns    NoCoulomb/NoIFI on L4, plus LO/NLO on L6
#   * acceptance       xs histogram / _XSEC, which is 1 by construction here
#   * rate             events/s per rung, which is what sizes the arrays
set -euo pipefail

cd "$(dirname "$0")/.."
source ./scripts/job_config.sh
source ./scripts/card_sets.sh

SMOKE_EVT_L0_BORN=${SMOKE_EVT_L0_BORN:-20000}
SMOKE_EVT_L1_ISR=${SMOKE_EVT_L1_ISR:-20000}
SMOKE_EVT_L4_MASTER=${SMOKE_EVT_L4_MASTER:-20000}
SMOKE_EVT_L5_LO=${SMOKE_EVT_L5_LO:-20000}
# L6 at ~5.5 ev/s on a laptop core: 2000 events is minutes, and already enough
# to see the NLO columns and to time the rate.
SMOKE_EVT_L6_NLO=${SMOKE_EVT_L6_NLO:-2000}

smoke_evt() {
    case "$(card_rung "$1")" in
        L0_born)   echo "${SMOKE_EVT_L0_BORN}" ;;
        L1_isr)    echo "${SMOKE_EVT_L1_ISR}" ;;
        L4_master) echo "${SMOKE_EVT_L4_MASTER}" ;;
        L5_lo)     echo "${SMOKE_EVT_L5_LO}" ;;
        L6_nlo)    echo "${SMOKE_EVT_L6_NLO}" ;;
        *)         echo 20000 ;;
    esac
}

resolve_set smoke

case "${1:?usage: smoke.sh <integrate|run|report|clean>}" in

integrate)
    assert_sherpa_has_ladder_weights
    mkdir -p logs
    for c in "${cards[@]}"; do
        n=$(intcores_for "${c}")
        jid=$(sbatch --parsable -A "${ACCOUNT}" -p "${PARTITION}" -x "${EXCLUDE_NODES}" --ntasks="${n}" --export=ALL,CARD="${c}" \
                -o "logs/integrate_${c}_%j.out" -e "logs/integrate_${c}_%j.err" \
                scripts/integrate.sbatch)
        echo "integrate ${jid}: ${c} on ${n} core(s)"
    done
    echo
    echo "poll:  squeue -u \$USER"
    echo "then:  ./scripts/smoke.sh run"
    ;;

run)
    for c in "${cards[@]}"; do
        EVT=$(smoke_evt "${c}") WALL=02:00:00 ./scripts/submit.sh "${c}" 1
    done
    echo
    echo "poll:  squeue -u \$USER"
    echo "then:  ./scripts/smoke.sh report"
    ;;

report)
    rc=0
    echo "==================================================================="
    echo " smoke report -- 161.0 GeV, mue"
    echo "==================================================================="
    printf "%-20s %10s %10s %10s  %s\n" rung events "wall[s]" "ev/s" "sigma"
    for c in "${cards[@]}"; do
        log="logs/${c}/1.out"
        if [[ ! -f ${log} ]]; then
            printf "%-20s %10s\n" "$(card_rung "${c}")" "NO LOG"
            rc=1; continue
        fi
        evt=$(smoke_evt "${c}")
        dt=$(grep -m1 '^### generation took' "${log}" | awk '{print $4}' || true)
        rate=$(grep -m1 '^### generation took' "${log}" | sed -E 's/.*\(([0-9.]+) ev\/s\).*/\1/' || true)
        sig=$(grep -m1 -E '^    sigma = ' "${log}" | sed -E 's/^ *sigma = //' || true)
        printf "%-20s %10s %10s %10s  %s\n" \
               "$(card_rung "${c}")" "${evt}" "${dt:-?}" "${rate:-?}" "${sig:-?}"
    done

    echo
    echo "--- assertions ---"
    for c in "${cards[@]}"; do
        log="logs/${c}/1.out"
        [[ -f ${log} ]] || continue
        ok=1
        grep -qE '^### seed ok: 1$' "${log}" || { echo "  ${c}: SEED CHECK DID NOT PASS"; ok=0; rc=1; }
        grep -q 'no weight with given name' "${log}" && { echo "  ${c}: GenCrossSection weight-name error PRESENT"; ok=0; rc=1; }
        # Prints the whole unused list for review and fails only on the keys
        # whose loss would not show up any other way (see UNUSED_WATCHLIST).
        check_unused_settings "${log}" || { ok=0; rc=1; }
        need=$(required_columns "${c}")
        if [[ -n ${need} ]]; then
            for w in ${need}; do
                grep -qE "^    columns .*${w}( |$)" "${log}" \
                    || { echo "  ${c}: MISSING column ${w}"; ok=0; rc=1; }
            done
        fi
        grep -qE '^    acceptance = .* = 1\.0000' "${log}" \
            || { echo "  ${c}: acceptance is not 1 (see log)"; ok=0; rc=1; }
        grep -q '^### OK ' "${log}" || { echo "  ${c}: task did not reach the OK line"; ok=0; rc=1; }
        [[ ${ok} -eq 1 ]] && echo "  ${c}: all checks passed"
    done

    echo
    echo "--- dipole set (YFS: Dump_Dipoles), L4 ---"
    sed -n '/YFS dipole set/,/form factor =/p' "logs/L4_master_161p0_mue/1.out" 2>/dev/null | sed 's/^/  /' || true

    echo
    echo "--- sizing from the measured rates ---"
    echo "  Target 0.5% on sigma per configuration. For L0/L1/L4/MW that target is"
    echo "  met almost immediately, so their production event counts are set by the"
    echo "  statistics the HISTOGRAMS need -- particularly the high bins of"
    echo "  photon-E-log -- not by sigma. Only L6 is sigma-limited."
    for c in "${cards[@]}"; do
        log="logs/${c}/1.out"
        [[ -f ${log} ]] || continue
        rate=$(grep -m1 '^### generation took' "${log}" | sed -E 's/.*\(([0-9.]+) ev\/s\).*/\1/' || true)
        rel=$(grep -m1 -E '^    sigma = ' "${log}" | sed -E 's/.*\(([0-9.]+)%\).*/\1/' || true)
        evt=$(smoke_evt "${c}")
        [[ -n ${rate:-} && -n ${rel:-} ]] || continue
        # LC_ALL=C: awk otherwise formats with the node's decimal separator, and
        # a table of core-hours written with commas is both unreadable and
        # unparseable by anything downstream.
        prod=$(evt_for "${c}"); nsd=$(nseed_for "${c}")
        LC_ALL=C awk -v r="${rate}" -v e="${evt}" -v p="${rel}" -v c="$(card_rung "${c}")" \
                     -v pe="${prod}" -v ns="${nsd}" 'BEGIN{
            # The relative error scales as 1/sqrt(N) for a weighted sample too,
            # provided the weight spread is unchanged -- which it is: same card,
            # same integration grid, only more events.
            need = e * (p/0.5)^2;
            printf "  %-11s %8.4g ev/s | %5.2f%% at %7d ev | 0.5%% needs %9.4g ev = %7.2f core-h\n",
                   c, r, p, e, need, need/r/3600.;
            printf "  %-11s production default %d seeds x %d ev = %.4g ev = %.2f core-h/config",
                   "", ns, pe, ns*pe, ns*pe/r/3600.;
            if (ns*pe < need)
                printf "   <-- BELOW the 0.5%% target, raise EVT_%s\n", toupper(c);
            else
                printf "  (%.2fx the 0.5%% target)\n", ns*pe/need;
        }'
    done
    echo
    echo "  Rates measured with the node to ourselves are ~40% optimistic:"
    echo "  a wave of single-core tasks packs ~22 to a 48-core node and"
    echo "  memory-bandwidth contention costs roughly that much. Size with headroom."
    exit ${rc}
    ;;

clean)
    ./scripts/reset.sh smoke "${@:2}"
    ;;

*) echo "usage: smoke.sh <integrate|run|report|clean>" >&2; exit 1 ;;
esac
