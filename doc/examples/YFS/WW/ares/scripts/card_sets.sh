#! /bin/bash
# SET spec -> card list. Sourced by submit-all.sh, integrate-all.sh, merge-all.sh
# and reset.sh; not executed directly. Requires job_config.sh first.
#
#   resolve_set <SPEC>    fills the array `cards`
#   set_names             prints the accepted aliases, for usage messages
#
# THE CARD LIST IS NEVER WRITTEN DOWN HERE. It is enumerated from
# cards/generated/*.yaml, which scripts/prepare.sh produces by running
# cards/make-cards.py on the cluster. Five scripts need the same mapping, and a
# hand-maintained copy in any of them drifts silently: a reset.sh that has
# fallen one family behind does not error, it leaves that family's yodas in
# place and the next merge averages two generator versions. Adding a rung to
# make-cards.py is enough for everything here to pick it up -- prepare.sh
# additionally checks that job_config.sh's matrix still accounts for every card
# make-cards.py wrote, so a new family cannot be silently unmatched by the
# filters either.
#
# SPEC grammar -- comma-separated filters, ANDed:
#   rung:L6_nlo        energy:161p0        channel:mue        mass:80p379
# or one of the aliases below, or a literal card name.
#
#   ./scripts/submit-all.sh all
#   ./scripts/submit-all.sh l6
#   ./scripts/submit-all.sh rung:L6_nlo,energy:161p0
#   ./scripts/submit-all.sh L6_nlo_161p0_mue

SET_ALIASES="all cheap l0 l1 l4 l5 l6 mw smoke e157p5 e161p0 e162p5"
set_names() { echo "${SET_ALIASES}"; }

_carddir() { echo "${CARDDIR:-cards/generated}"; }

# Every card that exists, in a stable order.
all_cards() {
    local f
    for f in "$(_carddir)"/*.yaml; do
        [[ -e ${f} ]] || {
            echo "ERROR: no cards in $(_carddir)." >&2
            echo "       The cards are generated, not shipped: ./scripts/prepare.sh" >&2
            return 1; }
        basename "${f}" .yaml
    done | sort
}

_alias_to_spec() {
    case "$1" in
        all)     echo "" ;;                       # no filter
        l0)      echo "rung:L0_born" ;;
        l1)      echo "rung:L1_isr" ;;
        l4)      echo "rung:L4_master" ;;
        l5)      echo "rung:L5_lo" ;;
        l6)      echo "rung:L6_nlo" ;;
        mw)      echo "rung:MW" ;;
        # Everything except the expensive NLO rung. This is the set that can be
        # run to completion in under an hour of wall time.
        cheap)   echo "@not:rung:L6_nlo" ;;
        # One card per rung at the reference point, for scripts/smoke.sbatch.
        smoke)   echo "@list:L0_born_161p0_mue L1_isr_161p0_mue L4_master_161p0_mue L5_lo_161p0_mue L6_nlo_161p0_mue" ;;
        e157p5)  echo "energy:157p5" ;;
        e161p0)  echo "energy:161p0" ;;
        e162p5)  echo "energy:162p5" ;;
        *)       return 1 ;;
    esac
    return 0
}

resolve_set() {
    local spec=$1 expanded c f key val keep all
    cards=()
    [[ -n ${spec} ]] || { echo "ERROR: empty SET spec. One of: $(set_names)" >&2; return 1; }

    if expanded=$(_alias_to_spec "${spec}"); then spec=${expanded}; fi

    # Explicit list (the smoke alias).
    if [[ ${spec} == @list:* ]]; then
        for c in ${spec#@list:}; do cards+=( "${c}" ); done
    else
        all=$(all_cards) || return 1
        for c in ${all}; do
            keep=1
            if [[ ${spec} == @not:* ]]; then
                # single negated filter
                key=${spec#@not:}; val=${key#*:}; key=${key%%:*}
                key=${key#@not:}
                case "${key}" in
                    rung)    [[ $(card_rung "${c}")    == "${val}" ]] && keep=0 ;;
                    energy)  [[ $(card_energy "${c}")  == "${val}" ]] && keep=0 ;;
                    channel) [[ $(card_channel "${c}") == "${val}" ]] && keep=0 ;;
                    mass)    [[ $(card_mass "${c}")    == "${val}" ]] && keep=0 ;;
                    *) echo "ERROR: unknown filter key '${key}'" >&2; return 1 ;;
                esac
            elif [[ -n ${spec} ]]; then
                local IFS=,
                for f in ${spec}; do
                    key=${f%%:*}; val=${f#*:}
                    case "${key}" in
                        rung)    [[ $(card_rung "${c}")    == "${val}" ]] || keep=0 ;;
                        energy)  [[ $(card_energy "${c}")  == "${val}" ]] || keep=0 ;;
                        channel) [[ $(card_channel "${c}") == "${val}" ]] || keep=0 ;;
                        mass)    [[ $(card_mass "${c}")    == "${val}" ]] || keep=0 ;;
                        *)
                            # Not a filter -- treat the whole spec as a literal
                            # card name and let the existence check below decide.
                            cards=( "${spec}" ); keep=-1 ;;
                    esac
                done
                [[ ${keep} == -1 ]] && break
            fi
            [[ ${keep} == 1 ]] && cards+=( "${c}" )
        done
    fi

    if [[ ${#cards[@]} -eq 0 ]]; then
        echo "ERROR: SET spec '${1}' matched no cards." >&2
        echo "       aliases: $(set_names)" >&2
        echo "       filters: rung:<${RUNGS_ALL// /|}> energy:<${ENERGIES_ALL// /|}> channel:<...> mass:<...>" >&2
        return 1
    fi

    # Every resolved name must be a real card and must parse. A typo that
    # silently resolves to nothing, or to a card in the wrong family, is the
    # failure mode this guards.
    for c in "${cards[@]}"; do
        [[ -f "$(_carddir)/${c}.yaml" ]] || {
            echo "ERROR: no such card: $(_carddir)/${c}.yaml" >&2; return 1; }
        card_parse "${c}" >/dev/null || return 1
    done
    return 0
}
