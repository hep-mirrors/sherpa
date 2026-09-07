# YFS corrections to leptonic WW at the FCC-ee threshold

A correction budget for `e+e- -> W+W- -> 4 leptons` across the FCC-ee W-mass
threshold scan: what each ingredient of the YFS treatment is worth in the total
cross section, and what that is worth in `m_W`.

FCC-ee aims to measure `m_W` to <= 0.5 MeV from the threshold cross section,
with a theory target of ~0.01% on `sigma(WW)`. That is a factor ~10 beyond LEP2,
and beyond what KandY (YFSWW3+KoralW) or RACOONWW deliver.

## The ladder

Six rungs in four runs. `YFS: Ladder_Weights: 1` emits two of them as named
weights on the L4 sample, so their ratios are exactly correlated and carry
almost no MC error of their own.

| rung | what it adds | where it comes from |
|---|---|---|
| L0 | Born | `L0_born` run (no `YFS:` block) |
| L1 | + ISR | `L1_isr` run (`MODE: ISR`) |
| L2 | + FSR | L4 x `YFS.NoCoulomb` x `YFS.NoIFI` |
| L3 | + IFI | L4 x `YFS.NoCoulomb` |
| L4 | + Coulomb | L4 nominal |
| L6 | + NLO matching | `L6_nlo` run (`NLO_Part: BVR`) |

Generate the cards with `cards/make-cards.py`, analyse with `budget.py`.

```
cd cards && python3 make-cards.py
# ... run cards/generated/*.yaml, one directory per rung ...
python3 budget.py <rundir>
```

`budget.py` expects `<rundir>/{L0_born,L1_isr,L4_master,L6_nlo}` and, for the
`m_W` conversion, `<rundir>/mw/<mass>` holding the three-point scan.

Cards carry one **common `SELECTORS` block** -- identical on every rung and
every channel, because every rung has to see identical phase space or the cross
sections are not comparable. It is a regulator, not the acceptance: charged
leptons at `10-170 deg` and `E > 5 GeV`, well outside the Rivet acceptance
(`18.2 deg`, `10 GeV`), which is where the acceptance still lives.

Run cut-free, the four channels containing an electron pick up non-resonant
single-`W`/t-channel-photon diagrams whose final-state `e` goes collinear to the
beam. Measured on 2M-event samples that region put `sigma(mue)` at
`0.0735 pb` against `0.0517 pb` for `mutau` -- **40% of the "WW" cross section
was not WW** -- and threw rare enormous weights: one seed in four read 1.2-2.5%
error where its siblings read 0.07%, with a central value up to 2.4% high.
`rivet-merge` combines seeds by `sumW`, not by inverse variance, so that seed
dragged the merged number with it, worth ~16 MeV in `m_W`. The electron-free
channels never had the problem (four seeds agreeing to four decimals at 0.07%).

## Event generation mode: run weighted

`EVENT_GENERATION_MODE` defaults to `PartiallyUnweighted`, and the unweighting
efficiency on this process is terrible:

| rung | unweighting efficiency |
|---|---|
| L0 born | 0.274% |
| L1 isr | 0.689% |
| L4 master | 0.243% |
| L6 nlo | **0.0658%** |

L6 discards ~1500 trial points per accepted event, which is why it manages only
~456 accepted events/day on one core. But this study wants cross sections and
weighted histograms, and neither needs unweighted events. The same L6 run's
**integration** had already delivered

```
2_4__e-__e+__ve__vmub__mu-__e+__EW(BVR) : 0.0605106 pb +- 3.68 %
```

in 568 s, against 0.0600 pb for L4 -- i.e. the NLO matching is worth about
+0.9%, from a number that was sitting in the integration output while event
generation ground away at a 14-day ETA.

So set `EVENT_GENERATION_MODE: Weighted` for production. The figure of merit to
check when doing so is the achieved relative error on sigma per unit wall time,
not events/day: weighted events carry a weight spread, so their effective
statistics per event is below one, and YFS weights do have tails (the
`CheckStability` skips). Measure it, do not assume it.

## Running these locally

Two things that cost time to rediscover:

- **The seed flag is `-R`, not `-s`.** `-s` is silently accepted and does
  nothing; the banner prints `Seed: 1` and every run is identical. Check the
  `Seed:` line in the log before believing that two runs are independent. A
  common seed is what you want for the correlated comparisons (scheme A vs B,
  the m_W derivative); it is not what you want when accumulating statistics.
- **Background jobs must outlive the launching shell.** `cmd &` inside a shell
  that then exits gets the job killed part-way with no error in the log -- it
  simply stops mid-event-generation and writes no summary.

For statistics, run several seeds in parallel and let `budget.py` combine them:
point it at a rung directory holding one subdirectory per seed and it takes the
inverse-variance mean per column.

## dsigma/dm_W

Taken from the generator, not the literature. Three ISR-only runs at
`m_W = 80.279 / 80.379 / 80.479`, central-differenced. `EW_SCHEME` defaults to
`Gmu`, in which `m_W` is an input and `sin^2(theta_W)` and `alpha` follow from
it, so overriding `PARTICLE_DATA: 24: Mass` is a consistent SM variation rather
than a bare propagator shift.

## Rivet analysis

`yfs-ww` (in `~/Documents/research/rivet-analysis/yfs-ww.cc`). Flavour-blind:
one negative and one positive charged lepton of any flavour, so all six
different-flavour channels fill the same histograms.

The observable it exists for is **`photon-E-log`**: the inclusive photon energy
spectrum, log-binned with `Gamma_W = 2.085 GeV` forced onto a bin edge.
arXiv:1906.09071 shows the photon distribution "transmutes" across
`E_gamma ~ Gamma_W` -- below it the radiation field sees the four decay
fermions, above it a single charged W -- and that no EEX scheme reproduces the
transition.

Every observable is also filled for the subsample with total photon energy below
`OMEGA` (default 2 GeV, roughly both the FCC-ee photon resolution and
`Gamma_W`). The same paper warns the non-factorisable interferences grow
strongly under exactly such mild cuts, so an inclusive-only budget would
understate them. (The "factors of 2-5" figure quoted around this study is our
own reading of that paper and has not been checked line-by-line against its
body text -- treat it as an order of magnitude, not a citation.)

### Cut families

Each core observable is filled up to six ways, from one production:

| suffix | what it selects |
|---|---|
| *(none)* | inclusive |
| `-fid` | both dressed leptons in the acceptance: `LEPCOS` (default `\|cos θ_l\| < 0.95`), `LEPE` (default `E_l > 10` GeV) |
| `-cut` | inclusive, total photon energy below `OMEGA` |
| `-fid-cut` | both of the above |
| `-sel` | inclusive, passing the missing-ET selection |
| `-fid-sel` | acceptance and the missing-ET selection |

`-sel` is the LEP2 idiom for this final state: a cut on `p_T^miss/E_beam`
(`PTMISS_EBEAM`, default 0.05 -- OPAL's acoplanar-dilepton selection,
hep-ex/9909052), plus an optional acollinearity veto (`ACOL_MAX`, **off** by
default). It is a ratio and not a GeV number so one value covers all three scan
energies. The acollinearity default is off deliberately: LEP experiments ranged
from ALEPH's 20 deg through DELPHI's 50 deg to L3's 165 deg, and OPAL's
simplified phase-space definition applies none, so any default here would be an
arbitrary choice dressed up as a convention.

The photon veto and the missing-ET selection are **not** crossed with each
other -- they answer different questions (scheme sensitivity vs. experimental
acceptance) and the joint subsample is too small to say anything at these
statistics.

**These are signal-only acceptances.** The production is `WW -> 4 leptons` and
nothing else: no gamma-gamma, no Bhabha-like `Zee`, no `ZZ`, no single-`W`, no
`tau tau`. In a real analysis the missing-ET cut is what suppresses exactly
those; here it only measures how much signal a realistic cut keeps. Quote
`-sel` numbers as a fiducial convention for comparison with an experimental
analysis, never as evidence the backgrounds are under control.

Do **not** set these as Rivet analysis options in the run cards
(`- yfs-ww:PTMISS_EBEAM=0.05`). Rivet accepts the syntax but folds the option
string into the analysis name, so every histogram path becomes
`/yfs-ww:PTMISS_EBEAM=0.05/...` and `scripts/check-yoda.py` -- which reads
`/yfs-ww/xs` and treats its absence as "the analysis did not run" -- fails
every card. Change the defaults in `yfs-ww.cc` instead.

Run the tau channels with `PARTICLE_DATA: 15: {Stable: 1}` (the generated cards
do) so the histograms measure WW radiation and not tau decay.

## Diagnostics

- `YFS: Dump_Dipoles: 1` prints the dipole set once per scheme state: flavours,
  charges, `ChargeNorm`, which pair radiates, the soft cutoff each form factor
  uses, and each dipole's own term `Y` in `FormFactorSum()`. Which pairs the
  construction picked is not obvious from the process definition and it
  determines the whole radiation pattern.
- The run reports how many events fell inside the leading-pole window and how
  many production recoils were rejected.

For the flat scheme on `e+e- -> mu- nubar_mu e+ nu_e` the dump shows **one** FF
dipole, `(mu-, e+)` -- a pair spanning the two W's. Each W has exactly one
charged decay product, so there is no within-W dipole at all. That is the flat
scheme working as designed, not a bug, but it is worth seeing.

## WW_Scheme: flat vs pole

`YFS: WW_Scheme` selects how a `W+W- -> 4f` final state is exponentiated.

- `flat` (default): the four external fermions are the emitters. Correct for
  photons softer than `Gamma_W`, where the W lives too briefly to be resolved.
- `pole`: the pole expansion -- a production stage with the W's as charged
  emitters plus a decay stage per W. The YFSWW3 picture, correct for photons
  harder than `Gamma_W`.

Neither is right across `E_gamma ~ Gamma_W`; the difference between them is the
scheme uncertainty.

**`pole` is implemented but not finished, and is not part of the production
budget.** What works: the dipole set is built correctly (production stage
`(e-,e+)`, `(W-,W+)` and four `(beam,W)`; decay stage `(W,l)` per W), each stage
is separately charge conserving so the photon mass cancels stage by stage, the
Coulomb singularity is subtracted per eq (9) of hep-ph/9606429, the decay stage
uses eq (12) of hep-ph/0302065 in closed form (the generic parametrisation is
singular there -- a two-body decay has `Q^2 = 0` identically), and the
production recoil is carried down to the four fermions. The production-stage
interference terms cancel to 0.4% of their individual size, which is the
`Gamma_W/M_W` suppression falling out of the construction.

What does not: the FSR crude/exact weight bookkeeping assumes the radiating
dipole's legs are event particles, and in the pole scheme they are W's. The
kinematic recoil is done; the matching Jacobian and crude normalisation are not.
At 161 GeV that leaves `pole` a factor ~2.7 below `flat`, which is far too large
to be a scheme difference.

`WW_FORM`/`BVV_WW` is the **superseded** predecessor and should not be used: it
omits the initial-initial pair, so its photon-mass logarithm does not cancel and
the factor came out ~3300x. It is retired from the live path.

## References

- S. Jadach, W. Placzek, M. Skrzypek, *QED exponentiation for quasi-stable
  charged particles: the e-e+ -> W-W+ process*, arXiv:1906.09071.
- S. Jadach, W. Placzek, M. Skrzypek, B.F.L. Ward, *Gauge invariant YFS
  exponentiation of (un)stable W+W- production*, hep-ph/9606429 -- Appendix A
  has the numerically stable closed forms for `ReB(s)`, `ReB(t)` and `Btilde`.
- W. Placzek, S. Jadach, *Multiphoton radiation in leptonic W-boson decays*,
  hep-ph/0302065 -- the decay stage.

Cut conventions:

- OPAL collaboration, *W+W- production and triple gauge boson couplings at LEP
  energies up to 183 GeV*, hep-ex/9909052 -- the acoplanar-dilepton selection
  the `-sel` family follows: `p_T^miss/E_beam > 0.05`, with a tighter 0.022
  variant under an acoplanarity-dependent scaling that is not reproduced here.
- ALEPH, DELPHI, L3, OPAL and the LEP EW Working Group, *Electroweak
  measurements in electron-positron collisions at W-boson-pair energies at
  LEP*, arXiv:1302.3415 -- the cross-experiment cut tables. Worth reading for
  the definitive per-experiment numbers; the acollinearity spread quoted above
  is from secondary sources and has not been checked against it directly.
- CEPC/FCC-ee threshold-scan study, arXiv:1812.09855 -- confirms 157-163 GeV as
  the operating window and sets an effective background level of ~0.3 pb, but
  for the semi-leptonic `2 jets + lepton + MET` topology rather than this one.
  Use it for the background list, not for cut values.
