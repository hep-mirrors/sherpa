#! /bin/bash
# Shared environment for the FCC-ee threshold YFS WW correction-budget jobs on
# Ares. Sourced by every *.sbatch and every driver script here; never executed
# directly.
#
# Campaign layout on Ares (${ARES_DIR}):
#   cards/make-cards.py           card generator, rsynced from the workstation
#   cards/generated/<CARD>.yaml   the 81 cards, GENERATED HERE by scripts/prepare.sh
#   analysis/yfs-ww.cc            Rivet source, BUILT ON ARES into RivetYFSWW.so
#   proc/<CARD>/Process/          per-card Comix process setup (see SHERPA_CPP_PATH)
#   int-res/<CARD>[.zip]          per-card integration grids
#   yodas/<CARD>/<seed>.yoda.gz   one file per array task, published atomically
#   logs/<CARD>/<seed>.{out,err}  live logs
#   merged/<E>/<CH>/<RUNG>/       one merged yoda per config -- budget.py reads this
#   merged/<E>/_mw/<MASS>/        the m_W scan, symlinked in as <E>/<CH>/mw

TOPDIR=${SLURM_SUBMIT_DIR:-${PWD}}

# ---------------------------------------------------------------------------
# Allocation
# ---------------------------------------------------------------------------
# plgsherpalumi-cpu runs to 2027-02-15, 72 h job limit. As of 2026-08-30 it had
# consumed 417 675 of 500 000 h, i.e. ~82 000 h left -- comfortable for this
# campaign (~400 h) but not for a careless resubmission of everything.
# plgpions-cpu is exhausted; plgmuone*-cpu no longer exist.
ACCOUNT=plgsherpalumi-cpu
export ACCOUNT

# Every sbatch carries this. plgrid for production AND builds: plgrid-now and
# plgrid-testing have MaxSubmitJobsPU=1, so anything queued there serialises.
PARTITION=plgrid
EXCLUDE_NODES='ac[0770-0780],ac0714'
export PARTITION EXCLUDE_NODES

ARES_DIR=/net/afscra/people/plgaprice/yfs-ww
export ARES_DIR

# ---------------------------------------------------------------------------
# Toolchain
# ---------------------------------------------------------------------------
# ~/.bashrc brings in the module function and rivetenv (Rivet 4 + HepMC3);
# ~/.load-modules adds libzip, python, openmpi, cmake -- all needed at runtime.
# Both reference unset variables, so relax `set -u` while sourcing (callers run
# with `set -euo pipefail`).
# YFSWW_NO_MODULES=1 skips this. It exists ONLY so the driver scripts can be
# exercised off-cluster against stub sbatch/squeue; never set it in a job.
if [[ ${YFSWW_NO_MODULES:-0} != 1 ]]; then
    _had_u=0; [[ $- == *u* ]] && _had_u=1
    set +u
    # Silenced: the module set prints ~200 lines of load/unload/reload chatter,
    # which at 450 array tasks is most of what logs/ would contain and buries
    # the lines the assertions actually read. Discarding it is only safe
    # because the toolchain is asserted explicitly below -- a module that
    # failed to load then produces a named error instead of a missing command
    # 40 minutes into a job.
    source ~/.bashrc      >/dev/null 2>&1
    source ~/.load-modules >/dev/null 2>&1
    [[ ${_had_u} -eq 1 ]] && set -u
    unset _had_u

    _missing=""
    for _t in python3 rivet rivet-build rivet-merge mpirun sbatch; do
        command -v "${_t}" >/dev/null 2>&1 || _missing="${_missing} ${_t}"
    done
    if [[ -n ${_missing} ]]; then
        echo "ERROR: toolchain incomplete after loading modules:${_missing}" >&2
        echo "       Re-run with YFSWW_NO_MODULES=0 and inspect" >&2
        echo "       'source ~/.bashrc; source ~/.load-modules' by hand." >&2
        return 1 2>/dev/null || exit 1
    fi
    unset _missing _t
fi

# Sherpa. The branch is ap-port-yfs-nlo. THE BUILD IS NOT OWNED BY THESE
# SCRIPTS: the worktree under plggyfsteam is shared, and the source there is
# updated and rebuilt by hand. Nothing in this directory pushes source or
# starts a build; everything here only asserts that the binary it is about to
# use is new enough.
#
# The cards need the YFS settings Ladder_Weights, Dump_Dipoles and
# NLO_Weight_Breakdown. assert_sherpa_has_ladder_weights below is the guard
# against running against a binary that predates them: without it Sherpa
# accepts the unknown YFS keys, ignores them, emits no YFS.* columns, and
# budget.py silently falls back to the nominal for every rung -- reporting
# L2 = L3 = L4 in a table that looks entirely reasonable.
export SHERPA_BRANCH=ap-port-yfs-nlo
export SHERPA_WORKTREE=/net/pr2/projects/plgrid/plggyfsteam/software/sherpa-worktree/${SHERPA_BRANCH}
# CMAKE_INSTALL_PREFIX for this build IS the build directory, so `--target
# install` is what publishes bin/Sherpa. A rebuild without the install step is
# a silent no-op. build-sherpa.sbatch asserts on that.
export sherpa=${SHERPA_WORKTREE}/build/bin/Sherpa

# OpenLoops. The L6 cards use Loop_Generator/Real_Generator: OpenLoops on
# 11 -11 -> 4 leptons, which needs libopenloops_eellll_ew_lt.so. That is in the
# stock group OpenLoops (verified 2026-08-30), which is also what this Sherpa
# was configured against (OPENLOOPS_PREFIX in build/CMakeCache.txt). OL-sqed is
# the scalar-QED build used by the pion campaign and is NOT wanted here.
export OL_PREFIX=/net/pr2/projects/plgrid/plggyfsteam/software/OpenLoops
export LD_LIBRARY_PATH=${OL_PREFIX}/lib:${LD_LIBRARY_PATH:-}

# Sherpa concatenates the SHERPA_CPP_PATH environment variable and the
# SHERPA_CPP_PATH setting (Run_Parameter.C:143 seeds the variable from the
# environment, line 256 appends the setting to it). We pass the setting on the
# command line, so an inherited environment value would produce a spliced
# nonsense path. Clear it.
unset SHERPA_CPP_PATH

# One analysis, from one source file.
export RIVET_ANALYSES="yfs-ww"
export RIVET_ANALYSIS_PATH=${TOPDIR}/analysis

# ---------------------------------------------------------------------------
# Card naming: <RUNG>_<ENERGY>_<CHANNEL>, or MW_<ENERGY>_<MASS>_mue
# ---------------------------------------------------------------------------
# Everything downstream keys off the card name, so parsing it wrong would
# quietly put a sample in the wrong merged directory. card_parse validates by
# reassembling the name and comparing.
RUNGS_ALL="L0_born L1_isr L4_master L5_lo L6_nlo MW"
ENERGIES_ALL="157p5 161p0 162p5"
CHANNELS_ALL="mue emu mutau taumu etau taue"
MASSES_ALL="80p279 80p379 80p479 80p4335"
export RUNGS_ALL ENERGIES_ALL CHANNELS_ALL MASSES_ALL

card_rung()    { case "$1" in MW_*) echo MW ;; *) echo "$1" | cut -d_ -f1-2 ;; esac; }
card_energy()  { case "$1" in MW_*) echo "$1" | cut -d_ -f2 ;; *) echo "$1" | cut -d_ -f3 ;; esac; }
card_channel() { echo "$1" | cut -d_ -f4 ; }
card_mass()    { case "$1" in MW_*) echo "$1" | cut -d_ -f3 ;; *) echo "" ;; esac; }

# Reassemble and compare. A card name this cannot round-trip is a hard error
# rather than something that lands in merged/<garbage>/.
card_parse() {
    local c=$1 r e ch m rebuilt
    r=$(card_rung "$c"); e=$(card_energy "$c"); ch=$(card_channel "$c"); m=$(card_mass "$c")
    if [[ ${r} == MW ]]; then rebuilt="MW_${e}_${m}_${ch}"; else rebuilt="${r}_${e}_${ch}"; fi
    if [[ ${rebuilt} != "${c}" ]]; then
        echo "ERROR: cannot parse card name '${c}' (round-tripped to '${rebuilt}')" >&2
        return 1
    fi
    echo "${RUNGS_ALL}"    | tr ' ' '\n' | grep -qx "${r}"  || { echo "ERROR: unknown rung '${r}' in '${c}'" >&2; return 1; }
    echo "${ENERGIES_ALL}" | tr ' ' '\n' | grep -qx "${e}"  || { echo "ERROR: unknown energy '${e}' in '${c}'" >&2; return 1; }
    echo "${CHANNELS_ALL}" | tr ' ' '\n' | grep -qx "${ch}" || { echo "ERROR: unknown channel '${ch}' in '${c}'" >&2; return 1; }
    if [[ ${r} == MW ]]; then
        echo "${MASSES_ALL}" | tr ' ' '\n' | grep -qx "${m}" || { echo "ERROR: unknown m_W '${m}' in '${c}'" >&2; return 1; }
    fi
    return 0
}

# Where merge.sbatch puts the merged yoda, relative to ${ARES_DIR}. budget.py
# wants <rundir>/{L0_born,L1_isr,L4_master,L6_nlo} and <rundir>/mw/<mass>, with
# <rundir> = merged/<energy>/<channel>. The m_W scan is mue-only, so it is
# merged once per energy into _mw/ and symlinked in as mw/ for every channel
# (scripts/link-mw.sh) -- dln(sigma)/dm_W is the same for all six.
merged_dir() {
    local c=$1 r e ch m
    r=$(card_rung "$c"); e=$(card_energy "$c"); ch=$(card_channel "$c"); m=$(card_mass "$c")
    if [[ ${r} == MW ]]; then echo "merged/${e}/_mw/${m}"; else echo "merged/${e}/${ch}/${r}"; fi
}

# ---------------------------------------------------------------------------
# Sizing, per rung
# ---------------------------------------------------------------------------
# EVT_* is events PER SEED, NSEED_* is array tasks per card, so total events per
# configuration = EVT * NSEED. WALL_* is the array task walltime request.
#
# PROVISIONAL until the Ares smoke test lands: the numbers below are laptop
# rates (L6 ~5.5 ev/s; L0/L1/L4 ~8600 ev/s) with a factor ~2 slowdown assumed
# for a packed Ares node. Measured single-node rates are ~40% optimistic once
# ~22 tasks share a 48-core node's memory bandwidth, so do not shave these.
#
# Override at submit time: EVT=... NSEED=... ./scripts/submit.sh <CARD>
EVT_L0_BORN=${EVT_L0_BORN:-2000000}    ; NSEED_L0_BORN=${NSEED_L0_BORN:-4}  ; WALL_L0_BORN=${WALL_L0_BORN:-08:00:00}
EVT_L1_ISR=${EVT_L1_ISR:-2000000}      ; NSEED_L1_ISR=${NSEED_L1_ISR:-4}    ; WALL_L1_ISR=${WALL_L1_ISR:-08:00:00}
EVT_L4_MASTER=${EVT_L4_MASTER:-2000000}; NSEED_L4_MASTER=${NSEED_L4_MASTER:-4}; WALL_L4_MASTER=${WALL_L4_MASTER:-08:00:00}
# L5_lo is L4 with BETA: 0 -- the LO baseline the whole ladder is now measured
# against, so it needs at least the statistics of the rungs compared to it.
# Same cost as L4: turning the beta expansion off changes a weight, not the
# amount of work per event.
EVT_L5_LO=${EVT_L5_LO:-2000000}        ; NSEED_L5_LO=${NSEED_L5_LO:-4}      ; WALL_L5_LO=${WALL_L5_LO:-08:00:00}
EVT_L6_NLO=${EVT_L6_NLO:-20000}        ; NSEED_L6_NLO=${NSEED_L6_NLO:-10}   ; WALL_L6_NLO=${WALL_L6_NLO:-24:00:00}
EVT_MW=${EVT_MW:-2000000}              ; NSEED_MW=${NSEED_MW:-4}            ; WALL_MW=${WALL_MW:-08:00:00}

# Cores for the MPI integration. Integration cannot be split by seed, so this is
# the one place MPI earns its keep. L6 integrates the BVR-matched process with
# OpenLoops and is the only slow one.
INTCORES_L6_NLO=${INTCORES_L6_NLO:-16}
INTCORES_OTHER=${INTCORES_OTHER:-4}

_rungvar() { # _rungvar <PREFIX> <CARD> -> value of <PREFIX>_<RUNG uppercased>
    local r; r=$(card_rung "$2" | tr '[:lower:]' '[:upper:]')
    eval "echo \"\${$1_${r}}\""
}
evt_for()      { _rungvar EVT "$1"; }
nseed_for()    { _rungvar NSEED "$1"; }
walltime_for() { _rungvar WALL "$1"; }
intcores_for() { [[ $(card_rung "$1") == L6_nlo ]] && echo "${INTCORES_L6_NLO}" || echo "${INTCORES_OTHER}"; }

# ---------------------------------------------------------------------------
# What each rung MUST emit as YFS.* named weights
# ---------------------------------------------------------------------------
# Derived from YFS_Handler::BuildNamedWeights / BuildLadderWeights, not guessed:
#   NoCoulomb  needs COULOMB: true
#   NoIFI      needs IFI_Sub: 1, FULL_FORM >= 1, TChannel 0, Fixed_Order != NLO
#   LO / NLO   need a non-Born NLO type with NLO available
# L0 has no YFS block at all and L1 has neither COULOMB nor IFI_Sub, so both
# legitimately emit nothing -- budget.py only reads their nominal.
# budget.py's L2 and L3 rungs are YFS.NoIFI and YFS.NoCoulomb on the L4 sample:
# if those go missing the budget does not fail, it silently reports L2 = L3 = L4.
# That is the whole reason this list exists.
required_columns() {
    case "$(card_rung "$1")" in
        # L5_lo carries COULOMB: true, IFI_Sub: 1 and FULL_FORM: 1 exactly as L4
        # does -- BETA: 0 changes the beta expansion, not the ladder factors --
        # so it must emit the same two columns. If it silently did not, budget.py
        # would build the ladder on an L5 whose L2/L3 rungs collapse onto it.
        L4_master|L5_lo) echo "YFS.NoCoulomb YFS.NoIFI" ;;
        L6_nlo)          echo "YFS.LO YFS.NLO YFS.NoCoulomb YFS.NoIFI" ;;
        *)               echo "" ;;
    esac
}

# ---------------------------------------------------------------------------
# Log assertions
# ---------------------------------------------------------------------------
# Sherpa accepts an unrecognised setting and then drops it, printing only a
# WARNING with the top-level keys involved. If the binary predates
# YFS: Ladder_Weights the whole YFS block shows up here and the run produces a
# plausible cross section with no ladder columns. Fail on the keys that matter;
# anything else is reported but tolerated.
# Only the keys whose loss would be INVISIBLE are fatal:
#   YFS                     -> no ladder columns, budget reports L2 = L3 = L4
#   PARTICLE_DATA           -> massless leptons, or a tau that decays after all
#   BEAMS / BEAM_ENERGIES   -> a perfectly good run at the wrong energy
#   EVENT_GENERATION_MODE   -> unweighted, 456 events instead of 200 000, which
#                              reads as a slow node rather than as a mistake
# PROCESSES and RIVET are deliberately NOT here: losing either aborts the run or
# produces no yoda, so they are already loud, and they are the ones most likely
# to be reported as partially-unused for a benign reason. A watchlist that fires
# on all 81 cards would be turned off within a day, and then the YFS case goes
# with it.
UNUSED_WATCHLIST="YFS PARTICLE_DATA BEAMS BEAM_ENERGIES EVENT_GENERATION_MODE"

check_unused_settings() { # check_unused_settings <logfile>
    local log=$1 bad=0 k
    grep -q "have not been used" "${log}" || return 0
    echo "### Sherpa reported UNUSED settings:"
    sed -n '/have not been used or include subsettings/,/Settings Report/p' "${log}" | sed 's/^/    /'
    for k in ${UNUSED_WATCHLIST}; do
        if sed -n '/have not been used or include subsettings/,/Settings Report/p' "${log}" \
             | grep -qE "^[[:space:]]*-[[:space:]].*:.*(^|[[:space:]])${k}([[:space:]]|$)"; then
            echo "ERROR: '${k}' was defined but NOT USED by Sherpa." >&2
            bad=1
        fi
    done
    [[ ${bad} -eq 0 ]] || {
        echo "ERROR: a setting this study varies was silently dropped. Refusing." >&2
        return 1
    }
    return 0
}

# ---------------------------------------------------------------------------
# Generator provenance
# ---------------------------------------------------------------------------
# The Sherpa worktree is SHARED and rebuilt by hand, so the binary can change
# between the integration of a card and the array that reads its grids. Grids
# from build A feeding events from build B is not an error anywhere: the run
# succeeds and the cross section looks fine. This is the stamp that makes it an
# error, and it doubles as the provenance record a paper needs -- "a sample
# whose provenance is unknown is not usable".
sherpa_fingerprint() {
    local lib=${SHERPA_WORKTREE}/build/lib64/SHERPA-MC/libYFSMain.so
    echo "sherpa_sha  $(sha256sum "${sherpa}" 2>/dev/null | awk '{print $1}')"
    echo "yfslib_sha  $(sha256sum "${lib}"    2>/dev/null | awk '{print $1}')"
    echo "git_commit  $(cd "${SHERPA_WORKTREE}" && git rev-parse HEAD 2>/dev/null)"
    echo "git_dirty   $(cd "${SHERPA_WORKTREE}" && git status --porcelain 2>/dev/null | grep -c . )"
}

# Compare against the stamp written at integration time. Only the two hashes
# are enforced: git_commit/git_dirty are recorded for the paper but a
# whitespace-only edit that does not change the binary must not invalidate
# perfectly good grids.
check_provenance() { # check_provenance <stampfile>
    local stamp=$1 k
    if [[ ! -f ${stamp} ]]; then
        echo "ERROR: no generator stamp at ${stamp}." >&2
        echo "       These grids predate provenance recording. Re-integrate the card" >&2
        echo "       so the events can be tied to a specific Sherpa build." >&2
        return 1
    fi
    for k in sherpa_sha yfslib_sha; do
        local was now
        was=$(awk -v k="${k}" '$1==k{print $2}' "${stamp}")
        now=$(sherpa_fingerprint | awk -v k="${k}" '$1==k{print $2}')
        if [[ ${was} != "${now}" ]]; then
            echo "ERROR: Sherpa has been rebuilt since these grids were integrated." >&2
            echo "       ${k} at integration: ${was:-<none>}" >&2
            echo "       ${k} now:            ${now:-<none>}" >&2
            echo "       Generating events now would mix two generator versions in one" >&2
            echo "       merge. Re-integrate this card (reset.sh <CARD> --yes --grids)." >&2
            return 1
        fi
    done
    return 0
}

# The seed flag is -R. -s is accepted and does nothing, and the banner then
# prints 'Seed: 1' for every task, i.e. N identical samples merged as if
# independent. Assert the banner agrees with what we asked for.
check_seed_line() { # check_seed_line <logfile> <expected seed>
    local log=$1 want=$2 got
    got=$(grep -m1 -E '^Seed: ' "${log}" | awk '{print $2}')
    if [[ -z ${got} ]]; then
        echo "ERROR: no 'Seed:' line in ${log} -- cannot confirm the RNG seed." >&2
        return 1
    fi
    if [[ ${got} != "${want}" ]]; then
        echo "ERROR: Sherpa reports 'Seed: ${got}' but we asked for -R ${want}." >&2
        return 1
    fi
    echo "### seed ok: ${got}"
    return 0
}

# The HepMC3 weight-name binding is fixed by the first event, so a weight name
# that appears later kills the run. That bug used to stop the NLO cards at ~486
# events; it is fixed, and this asserts it stays fixed.
check_no_weightname_error() { # check_no_weightname_error <logfile>
    if grep -q "no weight with given name" "$1"; then
        echo "ERROR: 'GenCrossSection::set_xsec: no weight with given name' in $1" >&2
        echo "       The YFS named-weight set is not constant across events." >&2
        return 1
    fi
    return 0
}

assert_sherpa_has_ladder_weights() {
    if ! grep -q "Ladder_Weights" "${SHERPA_WORKTREE}/YFS/Main/YFS_Base.C" 2>/dev/null; then
        echo "ERROR: ${SHERPA_WORKTREE} has no YFS: Ladder_Weights." >&2
        echo "       The worktree does not carry the YFS changes these cards need." >&2
        echo "       Update and rebuild it by hand before running this campaign." >&2
        return 1
    fi
    [[ -x ${sherpa} ]] || { echo "ERROR: no Sherpa binary at ${sherpa}" >&2; return 1; }

    # The decisive check: does the INSTALLED library actually contain the
    # setting? Grepping the source only proves someone edited a file. The YFS
    # code links into libYFSMain.so, NOT into bin/Sherpa -- a strings check on
    # the executable finds nothing even for a perfectly good build, so look in
    # the right place. grep -a rather than `strings` so this does not depend on
    # binutils being among the loaded modules.
    local yfslib=${SHERPA_WORKTREE}/build/lib64/SHERPA-MC/libYFSMain.so
    if [[ ! -f ${yfslib} ]]; then
        echo "ERROR: no installed ${yfslib}" >&2
        echo "       CMAKE_INSTALL_PREFIX for this build IS the build directory," >&2
        echo "       so the 'install' target is what publishes it." >&2
        return 1
    fi
    if ! grep -aq "Ladder_Weights" "${yfslib}"; then
        echo "ERROR: ${yfslib} does not contain 'Ladder_Weights'." >&2
        echo "       The source has it but the installed library does not, so the" >&2
        echo "       build did not run or ran without the install step. Sherpa would" >&2
        echo "       accept the YFS keys, drop them, and emit no ladder columns." >&2
        return 1
    fi
    # The install step is what publishes bin/Sherpa (CMAKE_INSTALL_PREFIX is the
    # build dir). A source file newer than the binary means the last build was
    # skipped or the install target was not run.
    if [[ ${SHERPA_WORKTREE}/YFS/Main/YFS_Handler.C -nt ${sherpa} ]]; then
        echo "ERROR: YFS/Main/YFS_Handler.C is NEWER than ${sherpa}." >&2
        echo "       The binary does not contain the current source -- either a build" >&2
        echo "       is still in flight, or one ran without the install step (this" >&2
        echo "       build's CMAKE_INSTALL_PREFIX IS the build directory, so 'install'" >&2
        echo "       is what publishes bin/Sherpa and skipping it is a silent no-op)." >&2
        return 1
    fi
    return 0
}
