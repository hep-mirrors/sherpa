"""Writes the per-experiment notebooks of examples/RemnantTests/SimpleLEP/XSecs.

    python3 make_notebooks.py        # (re)creates the *.ipynb listed in NOTEBOOKS

The notebooks only read generator output (see deuteron_common.py for the paths);
run them with Jupyter.  Edit the texts and code here or in the notebooks
themselves (re-running this script overwrites the notebooks).
"""
import json, os, sys

from paths import XSECS

def md(s):   return {"cell_type": "markdown", "metadata": {}, "source": s.strip("\n").splitlines(True)}
def code(s): return {"cell_type": "code", "metadata": {}, "execution_count": None, "outputs": [], "source": s.strip("\n").splitlines(True)}
def notebook(cells):
    return {"cells": cells, "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python"}}, "nbformat": 4, "nbformat_minor": 5}

SETUP = code(r'''
%matplotlib inline
import sys
sys.path.insert(0, "../../../../scripts")   # doc/scripts: deuteron_common.py, paths.py
from deuteron_common import *      # paths, readers, fits, figure style (see deuteron_common.py)
''')

GENERATORS = md(r'''
### The generator setups compared in this notebook

| name | what it is | output |
|---|---|---|
| **Sherpa s3_053** | Sherpa (branch hadronicRescattering = master + deuteron phase) with the setting found by the parameter scan in `../TuneScan` (round 3, point s3_053), scored against ALICE's 7 TeV $\pi$, $K$, $p$ and (anti)deuteron spectra; the same setting for all experiments; run cards `../../../../deuteron_tune_report/runcards/ScanTune_s3_053/` | `DeuteronFormationFromXSecs/scanTune_s3_053/` |
| **Sherpa s4_005** | the fine scan around s3_053 (`../TuneScan/scan_config_fine.yaml`, point s4_005): more semi-hard scatterings suppressed (`SIGMA_ND_NORM` 0.749), larger kick cap (`PT_MAX` 2.94), slightly fewer diquarks (`BARYON_FRACTION` 0.163); fewer protons above 3 GeV than s3_053; run cards `../../../../deuteron_tune_report/runcards/ScanTune_s4_005/` | `DeuteronFormationFromXSecs/scanTune_s4_005/` |
| **Sherpa s6_008**, **s6_010** | two further rounds around s4_005 (`../TuneScan/scan_config_leading.yaml`, `scan_config_kicks.yaml`): much larger transverse kicks in cluster decays (`KT_0` 2.51 / 3.18, `PT_MAX` 3.50 / 3.35; s3_053: 1.77, 2.60), `SIGMA_ND_NORM` 0.754 / 0.803, and the z shape of splittings with a perturbative quark (`ALPHA_H` 0.131, `BETA_H` 1.321, `GAMMA_H` 0.13 / 0.14); protons within 10 % of ALICE from 0.3 to 3 GeV; run cards `../../../../deuteron_tune_report/runcards/ScanTune_s6_008/`, `.../ScanTune_s6_010/`; at 7 TeV four seeds merged (400k events) | `DeuteronFormationFromXSecs/scanTune_s6_008/`, `.../scanTune_s6_010/` |
| **Pythia 8.318** | the Pythia 8 event generator with its default tune and its **own** implementation of the same deuteron model; how it is run is explained in `../PythiaReference/README.md` | `DeuteronFormationFromXSecs/pythia/` |

Older Sherpa setups (before tuning, the first tunes, master's defaults) stay available in `deuteron_common.ALL_MODELS`; add a name to `SHOWN` there to plot it again.

All three use the same deuteron model (Dal & Raklev, [arXiv:1504.07242](https://arxiv.org/abs/1504.07242)): every (anti)nucleon pair forms an (anti)deuteron with probability $P=\sigma(k)/\sigma_0$. The cross sections $\sigma(k)$ are fixed by measurements; $1/\sigma_0$ is the one free number. **All deuteron yields are proportional to $1/\sigma_0$**, so every generator's yields are first expressed "per unit of $1/\sigma_0$" and then $1/\sigma_0$ is fitted to each measurement. A good description means: after this fit the model follows the data in shape (small $\chi^2$), and $1/\sigma_0$ comes out similar for different experiments.
''')

NOTEBOOKS = {}

# =============================================================================
NOTEBOOKS["Model_Validation.ipynb"] = notebook([
md(r'''
# Model validation: the cross sections, the kinematics and the formation in the events

This notebook does not compare with any experiment. It checks that the Dal & Raklev model is implemented correctly in Sherpa, before any data comparison makes sense.

**The model in one paragraph.** After an event has been generated (hadronisation and hadron decays done), the `Deuteron_Formation` phase (`SHERPA/Single_Events/Deuteron_Formation.C`) looks at every pair of nucleons (or of antinucleons) in the event. For each pair it computes $k=|\vec p_1-\vec p_2|$ in the pair's rest frame, and for each of eight reactions $NN\to dX$ ($pn\to d\gamma$, $pp\to d\pi^+$, ...; the paper's Table I) the probability $P_i=\sigma_i(k)/\sigma_0$. $\sigma_i(k)$ are fits to measured cross sections (paper Eqs. (7)–(14), coded in `HADRON_RESCATTERING/XSecs/BaryonBaryon.C`). If the pair binds, the deuteron and the other products ($\gamma$, $\pi$, $\pi\pi$) are emitted isotropically in the pair's rest frame.

Sections: 1. cross sections; 2. final-state kinematics; 3. the nucleon pairs at the Z pole and which reactions make the deuterons; 4. deuterons actually formed vs. expected.
'''),
SETUP,
md(r'''
## 1. The cross sections $\sigma_i(k)$

**What the next cell shows.** The eight cross sections in $\mu$b versus $k$, as the event generation uses them. They come from Sherpa's self-test (`TEST: true` in the run card), which builds a nucleon pair of the right flavours at 400 values of $k$ and calls the same function the event generation calls, so thresholds and the conversion $k\to\eta=q/m_{\pi^+}$ for the $d\pi$ channels are included. Left: reactions of $pn$ pairs (paper's Fig. 4); right: $pp$ and $nn$ (Fig. 5); same axes as the paper.

**What to take from it.** The curves reproduce Figs. 4 and 5 of the paper: the $d\pi$ reactions dominate and peak at the $\Delta(1232)$ resonance near $k\approx1$ GeV; the $d\pi\pi$ reactions start at their thresholds ($k\approx1.04$–1.06 GeV). An independent re-implementation of the paper's formulas agrees with these numbers to about $10^{-6}$.
'''),
code(r'''
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
groups = [["pn_dgamma", "pn_dpi0", "pn_dpi0pi0", "pn_dpipluspiminus"],
          ["pp_dpiplus", "pp_dpipluspi0", "nn_dpiminus", "nn_dpiminuspi0"]]
for ax, group, title in zip(axes, groups, [r"cf. Fig. 4: $np\to dX$", r"cf. Fig. 5: $pp\to dX$, $nn\to dX$"]):
    for ch in group:
        x, y, _, _ = load_hist("test_sigma_" + ch, DATA_DIR)
        ax.plot(x, y, label=LABELS[ch])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 4.); ax.set_ylim(1., 5e3)
    ax.set_xlabel(r"$k$ [GeV]"); ax.set_ylabel(r"$\sigma$ [$\mu$b]"); ax.set_title(title); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
'''),
md(r'''
## 2. Final-state kinematics

**What the next cell shows.**
* Panels 1–3: for the three-body reactions ($d\pi\pi$) at $k=3$ GeV, the deuteron momentum $|p_d^*|$ in the pair's rest frame as the generator samples it ("sampled", uniform in the Dalitz plot as the paper prescribes, Sec. III D) and as flat three-body phase space gives it, computed independently ("phase space").
* Panel 4: the direction of the deuteron in the pair's rest frame, $\cos\theta^*$, for all eight reactions, normalised so that a flat distribution is 0.5.

**What to take from it.** "sampled" and "phase space" agree, and the $\cos\theta^*$ distributions are flat: the deuteron kinematics follow the paper.
'''),
code(r'''
fig, axes = plt.subplots(1, 4, figsize=(18, 3.8))
for ax, ch in zip(axes, ["pn_dpi0pi0", "pn_dpipluspiminus", "nn_dpiminuspi0"]):
    for name, lab in [("test_pd_cm_" + ch, "sampled"), ("test_pd_cm_" + ch + "_phasespace", "phase space")]:
        x, y, e, _ = load_hist(name, DATA_DIR); ax.plot(x, y, label=lab)
    ax.set_title(LABELS[ch]); ax.set_xlabel(r"$|p_d^*|$ [GeV]"); ax.legend(fontsize=8)
for ch in CHANNELS:
    x, y, e, _ = load_hist("test_costheta_d_" + ch, DATA_DIR); axes[3].plot(x, y, lw=0.8)
axes[3].set_ylim(0.4, 0.6); axes[3].set_xlabel(r"$\cos\theta_d^*$"); axes[3].set_title("isotropy, all reactions")
plt.tight_layout(); plt.show()
'''),
md(r'''
## 3. Which nucleon pairs make the deuterons (Z pole)

**What the next cell shows.** The paper's Fig. 6 for Sherpa in $e^+e^-\to Z\to$ hadrons: antinucleon pairs per event as a function of their momentum difference $k$ (rest frame of the pair, y axis) and their total momentum $|\vec p_1+\vec p_2|$ (lab, x axis), for $\bar p\bar n$, $\bar p\bar p$ and $\bar n\bar n$ pairs. The three maps share one colour scale (divided by the same maximum), as in the paper.

**What to take from it.** Most pairs have $k$ of 1–3 GeV, few have small $k$. There are about twice as many $\bar p\bar n$ pairs as $\bar p\bar p$ or $\bar n\bar n$ pairs (printed), as the paper states.
'''),
code(r'''
def load_map(name):
    d = np.loadtxt(os.path.join(DATA_DIR, "pairmap_" + name + ".dat"))
    n = int(round(np.sqrt(len(d))))
    m = d[:, 2].reshape(n, n)                               # [k bin, |p1+p2| bin], 0.1 GeV bins
    return m.reshape(n//2, 2, n//2, 2).sum(axis=(1, 3))     # 0.2 GeV bins, less noise
maps = {n: load_map(n) for n in ["pbarnbar", "pbarpbar", "nbarnbar"]}
vmax = max(m.max() for m in maps.values())
fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
for ax, (name, lab) in zip(axes, [("pbarnbar", r"$\bar p\bar n$"), ("pbarpbar", r"$\bar p\bar p$"), ("nbarnbar", r"$\bar n\bar n$")]):
    im = ax.imshow(maps[name]/vmax, origin="lower", extent=[0, 5, 0, 5], vmin=0, vmax=1, cmap="YlGnBu_r", aspect="auto")
    ax.set_xlabel(r"$|\vec p_1+\vec p_2|$ [GeV]"); ax.set_ylabel(r"$k$ [GeV]"); ax.set_title(lab + " pairs, cf. Fig. 6")
    ax.grid(False); ax.minorticks_off(); fig.colorbar(im, ax=ax)
plt.tight_layout(); plt.show()
tot = {n: m.sum() for n, m in maps.items()}
print(f"pbar nbar / pbar pbar = {tot['pbarnbar']/tot['pbarpbar']:.2f},  pbar nbar / nbar nbar = {tot['pbarnbar']/tot['nbarnbar']:.2f}   (paper: 'roughly a factor two')")
'''),
md(r'''
**What the next cell shows.**
* Left: pairs per event per GeV of $k$ (nucleon and antinucleon pairs together): how many pairs there are.
* Right: expected (anti)deuterons per event per GeV of $k$ for each reaction: the same pairs, each weighted with its probability $P_i(k)$.
* Table: the share of each reaction in the antideuteron yield.

**What to take from it.** Deuterons come mostly from pairs near the $\Delta$ resonance ($k\approx1$ GeV, dotted line) through the $d\pi$ reactions (~80 %); $d\pi\pi$ gives ~20 %, $d\gamma$ ~1 %.
'''),
code(r'''
summary, header = read_summary(DATA_DIR)
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
for pair in ["pn", "pp", "nn"]:
    x, y, e, _ = load_hist("pairs_k_" + pair, DATA_DIR); axes[0].plot(x, y, label=pair + " (+ antinucleon pairs)")
axes[0].set_xlabel(r"$k$ [GeV]"); axes[0].set_ylabel("pairs per event / GeV"); axes[0].legend(fontsize=8)
for ch in CHANNELS:
    x, y, e, _ = load_hist("expected_k_" + ch, DATA_DIR); axes[1].plot(x, y, label=LABELS[ch])
axes[1].set_yscale("log"); axes[1].set_ylim(1e-8, None); axes[1].axvline(1.0, color="k", ls=":", lw=1)
axes[1].set_xlabel(r"$k$ [GeV]"); axes[1].set_ylabel(r"expected $d+\bar d$ per event / GeV"); axes[1].legend(fontsize=7)
plt.tight_layout(); plt.show()
tot = sum(summary.get("expected_dbar_" + ch, (0, 0))[0] for ch in CHANNELS)
for ch in CHANNELS:
    v = summary.get("expected_dbar_" + ch, (0., 0.))[0]; print(f"{ch:20s} {v:10.3e} dbar per event   {v/tot:5.1%}")
'''),
md(r'''
## 4. Deuterons formed in the events vs. the expected yield

**Two ways to count deuterons.**
* **Expected** (used for all comparisons with data): for every pair and reaction the code adds the probability $P_i$ (the deuteron's momentum is drawn 10 times from the reaction's kinematics, each with weight $P_i/10$). Divided by the number of events and the bin width, this is the expected number of deuterons per event and unit, without any randomness, so the curves are smooth (paper Sec. III G, in the limit of many samples).
* **Formed**: the paper's procedure (Sec. III F): a random number per reaction, the pair binds if it is below $P_i$, used nucleons are not used again, and the deuteron is added to the event.

With $1/\sigma_0=1.8$/b only ~$4\times10^{-5}$ deuterons are formed per event, a handful in 200k events. A second run with $1/\sigma_0=18$/b forms ten times more; divided by 10 it must agree with the expected yield (yields are proportional to $1/\sigma_0$ as long as $P\ll1$).

**What the next cell shows.** Deuterons (top) and antideuterons (bottom) per event per unit of $p$, $p_T$ and rapidity $y$ at the Z pole: expected (line), formed at 18/b divided by 10 (circles), formed at 1.8/b (squares). Shaded: the ALEPH momentum window.

**What to take from it.** The circles follow the line, and the totals (printed) agree: the formation in the events is consistent with the expected yields.
'''),
code(r'''
binning = {"p": np.arange(0, 5.01, 0.25), "pT": np.arange(0, 5.01, 0.2), "y": np.arange(-4.8, 4.81, 0.4)}   # whole multiples of the Sherpa bins
fig, axes = plt.subplots(2, 3, figsize=(16, 8.5))
for row, sp in enumerate(["d", "dbar"]):
    for col, (var, xl) in enumerate([("p", r"$p$ [GeV]"), ("pT", r"$p_T$ [GeV]"), ("y", r"$y$")]):
        ax = axes[row, col]
        x, y, e = rebin(f"expected_{sp}_{var}", binning[var], DATA_DIR)
        ax.plot(x, y, label=r"expected, $1/\sigma_0=1.8$/b"); ax.fill_between(x, y-e, y+e, alpha=0.3)
        x2, y2, e2 = rebin(f"formed_{sp}_{var}", binning[var], VAL_DIR, scale=0.1)
        ax.errorbar(x2, y2, yerr=e2, fmt="o", ms=3, label=r"formed, $1/\sigma_0=18$/b, $\times0.1$")
        x3, y3, e3 = rebin(f"formed_{sp}_{var}", binning[var], DATA_DIR)
        ax.errorbar(x3[y3 > 0], y3[y3 > 0], yerr=e3[y3 > 0], fmt="s", ms=3, alpha=0.6, label=r"formed, $1/\sigma_0=1.8$/b")
        if var == "p": ax.axvspan(0.62, 1.03, color="grey", alpha=0.15, label="ALEPH window")
        ax.set_yscale("log"); ax.set_ylim(3e-3*y.max(), 5*y.max())
        ax.set_xlabel(xl); ax.set_ylabel(("d" if sp == "d" else r"$\bar d$") + " per event / unit"); ax.legend(fontsize=7)
plt.tight_layout(); plt.show()
vsum, _ = read_summary(VAL_DIR)
for sp in ["d", "dbar"]:
    f, fe = vsum.get(f"formed_{sp}_total", (0, 0)); x, xe = summary[f"expected_{sp}_total"]
    print(f"{sp:5s} per event: formed (18/b)/10 = {f/10:.3e} +- {fe/10:.1e},   expected (1.8/b) = {x:.3e} +- {xe:.1e}")
'''),
])

# =============================================================================
NOTEBOOKS["LEP_ALEPH_OPAL.ipynb"] = notebook([
md(r'''
# LEP (CERN): antideuterons in Z decays, ALEPH and OPAL

## The experiment in plain words

**LEP** collided electrons and positrons at CERN (1989–2000). At $\sqrt s=91.2$ GeV the beams annihilate into a **Z boson**, which mostly decays into a quark and an antiquark. These turn into two back-to-back jets of hadrons ($\pi$, $K$, $p$, $n$, ...). Such an event is a **"hadronic Z decay"**; LEP's experiments recorded millions of them, and all rates here are **per hadronic Z decay**.

Antideuterons are extremely rare in these events (a few per million) and hard to identify, because they look like other particles except in a narrow momentum range.

**ALEPH** (Phys. Lett. B 639 (2006) 192, [hep-ex/0604023](https://arxiv.org/abs/hep-ex/0604023)) identified antideuterons through their energy loss in the gas of their tracking chamber, which works only for momenta $0.62<p<1.03$ GeV, and only in the angular range $|\cos\theta|<0.95$ ($\theta$ = angle to the beam). Counting the candidates in this **momentum window** and correcting for efficiency gives
$$\frac{\#\,\bar d\ \text{in the window}}{\#\,\text{hadronic Z decays}}=(5.9\pm1.8\pm0.5)\times10^{-6}\quad(\text{statistical}\pm\text{systematic error}).$$
Antideuterons are used, not deuterons: deuterons are also knocked out of the detector material by the collisions, which would fake a signal.

**OPAL** (Z. Phys. C 67 (1995) 203) searched in $0.35<p<1.1$ GeV and found **no candidate**. With their detection efficiency $\epsilon=0.234$ and $N=1.64\times10^6$ hadronic Z decays, a model predicting $n$ antideuterons per event in that range predicts $s=\epsilon N n$ candidates; seeing none when $s$ were expected is unlikely if $s$ is large. The paper (Eqs. (21), (22)) adds $s$ to the $\chi^2$.

**What is compared.** Two numbers per model: the antideuterons per event in ALEPH's window, and the expected number of OPAL candidates. From them the paper fits $1/\sigma_0$ (Sec. IV E, Tables VII, VIII, row "$\bar d$, LEP": Herwig++ 1.93/b, Pythia 8.186 1.80/b). Also ALEPH's **coalescence parameter $B_2$**, which compares the antideuteron rate with the square of the antiproton rate.
'''),
GENERATORS,
SETUP,
md(r'''
## 1. The model numbers in the two windows

**What the next cell does.** For each generator it reads from `summary.dat` the antideuterons per hadronic Z event in the ALEPH window ($0.62<p<1.03$ GeV, $|\cos\theta|<0.95$) and in the OPAL window ($0.35<p<1.1$ GeV), and divides by the $1/\sigma_0$ of the run, giving `a` and `b` = yields **per unit of $1/\sigma_0$**. For orientation it also prints the values at $1/\sigma_0=1.8$/b and the number of OPAL candidates this would predict.

**What to take from it.** How many antideuterons each generator makes at LEP for the same model parameter; differences come from the generators' nucleons.
'''),
code(r'''
A, sA = 5.9e-6, np.hypot(1.8e-6, 0.5e-6)     # ALEPH: dbar per hadronic Z decay in its window, total error
eps, Nev = 0.234, 1.64e6                      # OPAL: efficiency and number of hadronic Z decays
lep = {}
print(f"{'model':14s} {'ALEPH window':>14s} {'OPAL window':>14s}   at 1/sigma_0 = 1.8/b: ALEPH window, OPAL candidates   [events]")
for m in MODELS:
    if not available(m, "LEP"): continue
    a, ae = deuteron_count(m, "LEP", "dbar_ALEPH_window")      # per unit of 1/sigma_0
    b, be = deuteron_count(m, "LEP", "dbar_OPAL_window")
    lep[m] = (a, ae, b, be)
    print(f"{m:28s} {a:14.3e} {b:14.3e}   {1.8*a:.2e}, {eps*Nev*1.8*b:.2f}   [{nevents(read_summary(model_path(m, 'LEP'))[1]):.2g}]")
print("(first two columns: dbar per hadronic Z event per unit of 1/sigma_0 [b])")
'''),
md(r'''
## 2. The paper's fit of $1/\sigma_0$

**What the next cell does.** For each generator, with $x=1/\sigma_0$:
$$\chi^2(x)=\frac{(A-a\,x)^2}{\sigma_A^2}+\epsilon N\,b\,x .$$
The first term compares with ALEPH ($A\pm\sigma_A$), the second is the number of candidates OPAL should have seen (it saw none). Because everything is linear in $x$, the minimum is analytic: $x=(A-\epsilon N b\,\sigma_A^2/(2a))/a$, with error $\sigma_A/a$. ALEPH alone would give $x=A/a$.

**What to take from it.** The fitted $1/\sigma_0$ and $\chi^2$ per generator, next to the paper's values. OPAL's null result always pulls $1/\sigma_0$ below the ALEPH-only value. With only two numbers, any model can be fitted; what matters is whether the value agrees with the other experiments (BaBar, ALICE notebooks).
'''),
code(r'''
lep_fit = {}
print(f"{'model':14s} {'ALEPH alone':>14s} {'ALEPH+OPAL':>20s} {'chi2':>6s}")
for m, (a, ae, b, be) in lep.items():
    x_aleph = A/a
    x = max(0., (A - eps*Nev*b*sA**2/(2*a))/a); xe = sA/a
    chi2 = (A - a*x)**2/sA**2 + eps*Nev*b*x
    lep_fit[m] = (x, xe, chi2)
    print(f"{m:28s} {x_aleph:8.2f} /b     {x:8.2f} +- {xe:4.2f} /b {chi2:6.2f}")
print("paper (Tables VII, VIII, row 'dbar, LEP'): Herwig++ 1.93/b (chi2 3.53), Pythia 8.186 1.80/b (chi2 3.55)")
'''),
md(r'''
## 3. The antideuteron momentum spectra

**What the next cell shows.** The antideuteron momentum spectrum in hadronic Z decays (all angles) for each generator, at its own best-fit $1/\sigma_0$ from the previous cell, in 0.1 GeV bins. Shaded: ALEPH's window. There is no measured spectrum: ALEPH measured one number, and OPAL saw no candidate. The black point is ALEPH's rate spread over the window, $A/(0.41~{\rm GeV})$, divided by 0.95 to go from $|\cos\theta|<0.95$ to all angles (assuming isotropic antideuterons). Apart from it, the figure shows only where in momentum each generator puts the antideuterons. Sherpa's curves are smooth (expected yields); Pythia's are formed by random draws and fluctuate more.
'''),
code(r'''
edges = np.arange(0., 3.001, 0.1)
fig, ax = plt.subplots(figsize=(8, 5))
for m in lep:
    x, y, e, w = deuteron_hist(m, "LEP", "dbar", "p")
    c = 0.5*(edges[1:]+edges[:-1])
    v = np.array([(y[(x > l) & (x < h)]*w).sum()/(h-l) for l, h in zip(edges[:-1], edges[1:])])
    ax.plot(c, lep_fit[m][0]*v, label=f"{m}, $1/\\sigma_0$ = {lep_fit[m][0]:.2f} b$^{{-1}}$", **STYLE[m])
ax.axvspan(0.62, 1.03, color="grey", alpha=0.15, label="ALEPH window")
ax.errorbar(0.825, A/0.41/0.95, yerr=sA/0.41/0.95, xerr=0.205, fmt="o", color="black", ms=6, zorder=5,
            label=r"ALEPH, $(5.9\pm1.9)\times10^{-6}$ in the window")
ax.set_xlabel(r"$p$ [GeV]"); ax.set_ylabel(r"$\bar d$ per hadronic Z decay per GeV"); ax.set_xlim(0, 3); ax.set_ylim(0, None)
ax.set_title(r"$e^+e^-\to Z\to$ hadrons, $\sqrt{s}=91.2$ GeV"); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
'''),
md(r'''
## 4. ALEPH's coalescence parameter $B_2$

**What it is.** $B_2$ compares the antideuteron yield with the square of the antiproton yield at half the antideuteron momentum (a deuteron is made of two nucleons, each carrying about half its momentum):
$$B_2=\frac{E_d\,d^3N_{\bar d}/dp_d^3}{\left(E_p\,d^3N_{\bar p}/dp_p^3\right)^2},\qquad p_p=p_d/2 .$$
$E\,d^3N/dp^3$ is the **invariant yield** (number per event per unit of momentum-space volume, multiplied by the energy). ALEPH (Eqs. (4), (5)) computes it from the counts in the window: $E\,d^3N/dp^3\simeq N\,E/(p^2\,\Delta\Omega\,\Delta p)$ at the window centre, with $\Delta\Omega=2\pi\times2\times0.95$ the solid angle of $|\cos\theta|<0.95$. ALEPH finds $B_2=(3.3\pm1.0\pm0.8)\times10^{-3}$ GeV².

**What the next cell does.** The same calculation for each generator: antideuterons in the ALEPH window at the best-fit $1/\sigma_0$, antiprotons (not from weak decays) at half the momenta ($0.31<p<0.515$ GeV, `nucleon_pbar_ALEPH_halfwindow`).

**What to take from it.** $B_2$ tests the formation relative to the generator's own antiprotons; if a generator has the right $B_2$ but the wrong antideuteron rate, the problem is its antiproton rate.
'''),
code(r'''
dOmega = 2*np.pi*2*0.95
def inv_window(n, m, pmin, pmax):
    p0 = 0.5*(pmin+pmax); E = np.sqrt(p0**2 + m**2)
    return n*E/(p0**2*dOmega*(pmax-pmin))
for m, (a, ae, b, be) in lep.items():
    s = read_summary(model_path(m, "LEP"))[0]
    npb, npbe = s.get("nucleon_pbar_ALEPH_halfwindow", (0., 0.))
    nd = a*lep_fit[m][0]
    B2 = inv_window(nd, 1.875613, 0.62, 1.03)/inv_window(npb, 0.938272, 0.31, 0.515)**2
    print(f"{m:28s}: antiprotons per event in 0.31 < p < 0.515 GeV = {npb:.4f},  B2 = {B2*1e3:.2f}e-3 GeV^2 (at its best-fit 1/sigma_0)")
print("ALEPH: B2 = (3.3 +- 1.0 +- 0.8)e-3 GeV^2")
'''),
])

# =============================================================================
NOTEBOOKS["BaBar_ARGUS_CLEO.ipynb"] = notebook([
md(r'''
# BaBar (SLAC), with the ARGUS and CLEO limits: antideuterons at $\sqrt s\approx10.6$ GeV

## The experiment in plain words

**BaBar** recorded $e^+e^-$ collisions at $\sqrt s=10.58$ GeV, the energy of the $\Upsilon(4S)$ resonance, a "B factory". Besides the resonance, a large share of the collisions are **"continuum" events**, $e^+e^-\to q\bar q$ directly: like LEP's Z decays, but at a ninth of the energy. BaBar measured antideuterons both in decays of $\Upsilon(nS)$ resonances and in the continuum (Phys. Rev. D 89 (2014) 111102, [arXiv:1403.4409](https://arxiv.org/abs/1403.4409)). **We simulate only the continuum** ($e^+e^-\to q\bar q$); Sherpa cannot yet hadronise $\Upsilon\to ggg$.

**What BaBar measured for the continuum** (HEPData [ins1286317](https://www.hepdata.net/record/ins1286317)):
* Table 5: the number of antideuterons per hadronic continuum event, $(3.01\pm0.13^{+0.37}_{-0.31})\times10^{-6}$ (BaBar measured $0.35<p<2.25$ GeV and extrapolated to all momenta with a fit to its spectrum);
* Table 4: the cross section $\sigma(e^+e^-\to\bar dX)=9.63$ fb.
The momentum spectrum (the paper's Fig. 9) is only published as a figure, so we compare with the rate.

**ARGUS** (Phys. Lett. B 236 (1990) 102) and **CLEO** (Phys. Rev. D 75 (2007) 012009) at about 10 GeV found no continuum antideuterons and give upper limits: fewer than $1.7\times10^{-5}$ per hadronic event (ARGUS) and $\sigma(e^+e^-\to\bar dX)<31$ fb (CLEO), both at 90 % confidence.

**What is compared.** The rate per event, from which $1/\sigma_0$ is fixed (paper Tables VII, VIII, row "$\bar d$, BABAR": Herwig++ 0.63/b, Pythia 8.186 1.13/b), the shape of the momentum spectrum between the generators, and the two limits.
'''),
GENERATORS,
SETUP,
md(r'''
## 1. The rate, and $1/\sigma_0$ from it

**What the next cell does.** For each generator: the total antideuterons per hadronic event, per unit of $1/\sigma_0$ (`a`), from `summary.dat`. Matching BaBar's rate $R$ requires $1/\sigma_0=R/a$, with the error from BaBar's error.

**What to take from it.** The $1/\sigma_0$ each generator needs at BaBar, to be compared with its LEP value (LEP notebook) and with the paper's values. A model that describes nucleons well should need similar values at LEP and BaBar, as both are $e^+e^-\to q\bar q$.
'''),
code(r'''
row = read_hepdata_blocks(os.path.join(HEPDATA, "BaBar_continuum_ratio_HEPData_ins1286317_Table5.csv"))[0][0]
R, R_err = row[1], np.hypot(0.5*(row[2]-row[3]), 0.5*(row[4]-row[5]))
babar = {}
print(f"BaBar: dbar per hadronic event = ({R*1e6:.2f} +- {R_err*1e6:.2f})e-6")
for m in MODELS:
    if not available(m, "BaBar_10.58GeV"): continue
    a, ae = deuteron_count(m, "BaBar_10.58GeV", "dbar_total")
    babar[m] = (a, R/a, R_err/a)
    print(f"{m:28s}: {a:.3e} dbar per event per unit 1/sigma_0  ->  1/sigma_0 = {R/a:.2f} +- {R_err/a:.2f} /b")
print("paper (Tables VII, VIII, row 'dbar, BABAR'): Herwig++ 0.63/b, Pythia 8.186 1.13/b (fitted to the spectrum)")
'''),
md(r'''
## 2. The momentum spectrum

**What the next cell shows.** The antideuteron momentum spectrum $\frac{1}{N}\frac{dN_{\bar d}}{dp}$ (= $\frac{d\sigma(\bar dX)/dp}{\sigma_{\rm had}}$, the quantity of the paper's Fig. 9) for each generator, normalised to BaBar's rate (each at its $1/\sigma_0$ from the previous cell), over BaBar's measured range in 0.1 GeV bins.

**What to take from it.** The shape each generator predicts; BaBar's measured spectrum (paper's Fig. 9) peaks at about 0.9 GeV.
'''),
code(r'''
edges = np.arange(0.35, 2.2501, 0.1)
fig, ax = plt.subplots(figsize=(8, 5))
for m, (a, x, xe) in babar.items():
    xs, y, e, w = deuteron_hist(m, "BaBar_10.58GeV", "dbar", "p")
    c = 0.5*(edges[1:]+edges[:-1])
    v = np.array([(y[(xs > l) & (xs < h)]*w).sum()/(h-l) for l, h in zip(edges[:-1], edges[1:])])
    ax.plot(c, x*v, label=f"{m}, $1/\\sigma_0$ = {x:.2f} b$^{{-1}}$", **STYLE[m])
ax.set_xlabel(r"$p$ [GeV]"); ax.set_ylabel(r"$\frac{1}{N_{\rm had}}\frac{dN_{\bar d}}{dp}$ [GeV$^{-1}$]")
ax.set_xlim(0.3, 2.3); ax.set_ylim(0, None); ax.legend(fontsize=8)
ax.set_title(r"$e^+e^-\to q\bar q$ at $\sqrt{s}=10.58$ GeV, normalised to BaBar's rate (cf. Fig. 9)")
plt.tight_layout(); plt.show()
'''),
md(r'''
## 3. ARGUS and CLEO limits

**What the next cell does.** At the $1/\sigma_0$ fitted to BaBar and at the LEP-like value 1.8/b: the antideuteron rate per event (ARGUS limit $<1.7\times10^{-5}$) and the cross section $\sigma(e^+e^-\to\bar dX)=$ rate × $\sigma_{\rm had}$, where $\sigma_{\rm had}=9.63\,{\rm fb}/3.01\times10^{-6}=3.2$ nb follows from BaBar's two tables (CLEO limit $<31$ fb).
'''),
code(r'''
t4 = read_hepdata_blocks(os.path.join(HEPDATA, "BaBar_continuum_ratio_HEPData_ins1286317_Table5.csv"))[0][0]
sig_had = 9.63/t4[1]                                    # fb
for m, (a, x, xe) in babar.items():
    for lab, inv in [("BaBar fit", x), ("1.8/b", 1.8)]:
        print(f"{m:28s} at {lab:9s}: rate {a*inv:.2e} (ARGUS < 1.7e-5),  sigma = {a*inv*sig_had:.1f} fb (CLEO < 31 fb)")
'''),
])

# =============================================================================
NOTEBOOKS["ALICE_2018_pp.ipynb"] = notebook([
md(r'''
# ALICE (CERN LHC): deuterons in proton–proton collisions at 0.9, 2.76 and 7 TeV

## The experiment in plain words

**ALICE** is one of the four big experiments at the LHC. We use its measurement of deuterons and antideuterons in proton–proton collisions at $\sqrt s=0.9$, 2.76 and 7 TeV (Phys. Rev. C 97 (2018) 024615, [arXiv:1709.08522](https://arxiv.org/abs/1709.08522), data on HEPData [ins1625294](https://www.hepdata.net/record/ins1625294)). ALICE identifies (anti)deuterons from their energy loss in the gas of its time-projection chamber and their flight time to a time-of-flight detector.

**What is measured, and the words used.**
* **Minimum bias / inelastic events.** ALICE records essentially every collision in which the protons break up ("inelastic"); yields are **per inelastic event**. Some of these are **diffractive** (one or both protons only get excited, few particles are made). Sherpa (Amisic) and our Pythia runs generate only **non-diffractive** (ND) collisions; as the paper, we assume diffractive events make no deuterons and multiply the ND yields by $f_{\rm ND}=1-\sigma_{\rm SD}/\sigma_{\rm inel}-\sigma_{\rm DD}/\sigma_{\rm inel}\approx0.68$ (SD/DD: single/double diffraction, measured by ALICE, HEPData ins1181770).
* **Transverse momentum $p_T$**: the momentum component perpendicular to the beams. **Rapidity $y$**: a measure of the velocity along the beams ($y=0$: at right angles to the beams). ALICE measures at $|y|<0.5$, the central region.
* **Invariant yield** $E\frac{d^3N}{dp^3}=\frac{1}{2\pi p_T}\frac{d^2N}{dy\,dp_T}$: deuterons per event per unit of $p_T$ and $y$, divided by $2\pi p_T$. This is the quantity of the paper's Fig. 7, in GeV$^{-2}$.
* **Statistical and systematic errors** are added in quadrature for the fits, and only the experimental errors are used (as the paper).

**Fitting.** For each energy and species ($d$, $\bar d$), $1/\sigma_0$ is fitted so that the model matches the data best (least squares). $\chi^2$ measures how well the **shape** then matches (roughly, $\chi^2\approx$ number of points means agreement within errors). The paper used *preliminary* ALICE data normalised to triggered events; we use the *final* data per inelastic event, so its $\chi^2$ and $1/\sigma_0$ (Tables VII, VIII) are not one-to-one comparable.

Sections: 1. the data; 2. deuteron spectra (Fig. 7); 3. $1/\sigma_0$ vs energy; 4. integrated yields and $\bar d/d$; 5. $B_2$; 6. the protons (Fig. 8); 7. fine scan around s3_053 (proton tail, 7 TeV).
'''),
GENERATORS,
SETUP,
md(r'''
## 1. What ALICE measured

**What the next cell shows.** Only the data: the deuteron (left) and antideuteron (right) invariant yields at the three energies, versus $p_T$; horizontal bars = bin widths, vertical bars = errors.

**What to take from it.** The spectra fall steeply with $p_T$; there are more deuterons at higher energy (more particles per collision), and somewhat fewer antideuterons than deuterons (the colliding protons carry baryon number).
'''),
code(r'''
fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), sharey=True)
for col, (sp, block, lab) in enumerate([("d", 0, "d"), ("dbar", 1, r"\bar d")]):
    for (rs, sqrts, table), c in zip(ENERGIES, ["C0", "C1", "C2"]):
        lo, hi, val, err = alice_spectrum(table, block)
        axes[col].errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), fmt="o", ms=4, color=c, label=f"{rs} TeV")
    axes[col].set_yscale("log"); axes[col].set_xlabel(r"$p_T$ [GeV]"); axes[col].set_title(f"${lab}$, ALICE, $|y|<0.5$"); axes[col].legend()
axes[0].set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]")
plt.tight_layout(); plt.show()
'''),
md(r'''
## 2. Deuteron spectra: the paper's Fig. 7

**What the next cell does.** For each energy and species, each generator's invariant yield per unit of $1/\sigma_0$ (×$f_{\rm ND}$) is averaged over ALICE's bins, and $1/\sigma_0$ is fitted (`fit_scale`). The plot shows the data and each generator at its best fit as a line; the legend gives $1/\sigma_0$ and $\chi^2$/points.

**What to take from it.** Sherpa before is too steep at 2.76 and 7 TeV (too high at low $p_T$, too low at high $p_T$); tuned Sherpa follows the data, with $\chi^2$ comparable to the paper's generators. Pythia's curves are jagged at high $p_T$ because it only forms deuterons by random draws.
'''),
code(r'''
fits = {}
fig, axes = plt.subplots(3, 2, figsize=(13, 13))
for row, (rs, sqrts, table) in enumerate(ENERGIES):
    fnd, _ = f_nd(sqrts); setup = f"ALICE_{rs}TeV"
    for col, (sp, block, lab) in enumerate([("d", 0, "d"), ("dbar", 1, r"\bar d")]):
        ax = axes[row, col]
        lo, hi, val, err = alice_spectrum(table, block)
        ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE", zorder=5, **STYLE["data"])
        for m in MODELS:
            if not available(m, setup): continue
            h = deuteron_hist(m, setup, sp, "invyield_midy")
            mc, _ = in_bins(h, lo, hi, scale=fnd)
            x, xe, chi2 = fit_scale(val, err, mc)
            fits[(m, rs, sp)] = (x, xe, chi2, len(val))
            xs, ys, _, _ = h; ok = ys > 0
            ax.plot(xs[ok], x*fnd*ys[ok], label=f"{m}: $1/\\sigma_0$ = {x:.2f} b$^{{-1}}$, $\\chi^2$ = {chi2:.1f}/{len(val)}", **STYLE[m])
        ax.set_yscale("log"); ax.set_xlim(0, 4); ax.set_ylim(1e-8, 1e-3)
        ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]")
        ax.set_title(f"${lab}$, pp at $\\sqrt{{s}}$ = {rs} TeV, $|y|<0.5$"); ax.legend(fontsize=8, loc="lower left")
plt.tight_layout(); savefig(fig, "alice_fig7"); plt.show()
'''),
md(r'''
## 3. $1/\sigma_0$ and $\chi^2$ versus energy

**What the next cell shows.** The fitted $1/\sigma_0$ for each generator versus $\sqrt s$, together with the paper's values (Pythia 8.186, Herwig++; preliminary data), and the table of all fits next to the paper's $\chi^2$.

**What to take from it.** In a perfect description $1/\sigma_0$ would be the same at all energies (it is meant to be a property of the deuteron, not of the collision). The paper's Pythia falls by ×1.36 from 0.9 to 7 TeV, our Pythia 8.318 by ×1.35; Sherpa before by ×3, tuned Sherpa by ×2.1: the tuning removed part of the energy dependence, not all of it.
'''),
code(r'''
paper = {"Pythia 8.186 (paper)": {("0.9","d"): (3.58, 0.77), ("0.9","dbar"): (3.63, 5.35), ("2.76","d"): (2.93, 9.22),
                                  ("2.76","dbar"): (2.88, 14.1), ("7","d"): (2.63, 55.1), ("7","dbar"): (2.58, 12.7)},
         "Herwig++ (paper)":     {("0.9","d"): (5.85, 0.71), ("0.9","dbar"): (6.15, 5.88), ("2.76","d"): (4.03, 10.9),
                                  ("2.76","dbar"): (4.00, 25.3), ("7","d"): (3.35, 231.), ("7","dbar"): (3.40, 97.4)}}
E = np.array([0.9, 2.76, 7.])
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
for ax, sp, lab in zip(axes, ["d", "dbar"], ["d", r"\bar d"]):
    for m in MODELS:
        v = [fits.get((m, rs, sp), (np.nan, 0))[:2] for rs, _, _ in ENERGIES]
        ax.errorbar(E, [a for a, b in v], yerr=[b for a, b in v], marker=MARKER[m], color=STYLE[m]["color"], ls=STYLE[m]["ls"], label=m)
    for name, mk in zip(paper, ["v", "o"]):
        ax.plot(E, [paper[name][(rs, sp)][0] for rs, _, _ in ENERGIES], mk, mfc="none", color="0.3", ls=":", label=name)
    ax.set_xscale("log"); ax.set_xlabel(r"$\sqrt{s}$ [TeV]"); ax.set_title(f"${lab}$"); ax.legend(fontsize=8)
axes[0].set_ylabel(r"fitted $1/\sigma_0$ [b$^{-1}$]"); axes[0].set_ylim(0, None)
plt.tight_layout(); plt.show()
print(f"{'':10s} {'':5s}" + "".join(f"{m:>30s}" for m in MODELS) + "   paper chi2: Herwig++, Pythia 8.186")
for rs, _, _ in ENERGIES:
    for sp in ("d", "dbar"):
        line = f"{rs+' TeV':10s} {sp:5s}"
        for m in MODELS:
            x, xe, c2, n = fits.get((m, rs, sp), (np.nan, np.nan, np.nan, 0)); line += f"         {x:5.2f}/b chi2 {c2:6.1f}/{n:<2d}"
        print(line + f"   {paper['Herwig++ (paper)'][(rs, sp)][1]:6.2f}, {paper['Pythia 8.186 (paper)'][(rs, sp)][1]:6.2f}")
'''),
md(r'''
## 4. Integrated yields and the antideuteron-to-deuteron ratio

**What the next cell does.** $dN/dy$ at $|y|<0.5$, i.e. deuterons per event per unit rapidity, summed over ALICE's **measured** bins only: $\sum_i 2\pi p_{T,i}\,\times$ invariant yield $\times\,\Delta p_{T,i}$ (ALICE's published $dN/dy$ adds an extrapolation to unmeasured $p_T$, 41–52 %). For the generators: at the $1/\sigma_0$ fitted to each spectrum. Right: $\bar d/d$, which does not depend on $\sigma_0$ at all (it cancels).

**What to take from it.** At the fitted $1/\sigma_0$ the yields agree by construction; $\bar d/d$ is a parameter-free test of the generators' (anti)baryon production. ALICE: 0.93, 0.89, 0.84 at 0.9, 2.76, 7 TeV; the generators agree within about two standard deviations (Sherpa somewhat high at the lower energies, Pythia somewhat low).
'''),
code(r'''
def dndy(val, lo, hi):
    xc = 0.5*(lo+hi); return (2*np.pi*xc*val*(hi-lo)).sum()
rows = {}
for rs, sqrts, table in ENERGIES:
    fnd, _ = f_nd(sqrts); setup = f"ALICE_{rs}TeV"
    for sp, block in [("d", 0), ("dbar", 1)]:
        lo, hi, val, err = alice_spectrum(table, block)
        rows[("ALICE", rs, sp)] = (dndy(val, lo, hi), np.sqrt(((2*np.pi*0.5*(lo+hi)*err*(hi-lo))**2).sum()))
        for m in MODELS:
            if not available(m, setup): continue
            mc, _ = in_bins(deuteron_hist(m, setup, sp, "invyield_midy"), lo, hi, scale=fnd)
            rows[(m, rs, sp)] = (dndy(mc, lo, hi), 0.)        # per unit 1/sigma_0
fig, ax = plt.subplots(figsize=(7, 4.5))
d_ = lambda who, rs: rows[(who, rs, "dbar")][0]/rows[(who, rs, "d")][0]
ax.errorbar(E, [d_("ALICE", rs) for rs, _, _ in ENERGIES],
            yerr=[d_("ALICE", rs)*np.hypot(rows[("ALICE", rs, "dbar")][1]/rows[("ALICE", rs, "dbar")][0], rows[("ALICE", rs, "d")][1]/rows[("ALICE", rs, "d")][0]) for rs, _, _ in ENERGIES],
            label="ALICE", **STYLE["data"])
for m in MODELS:
    have = [i for i, (rs, _, _) in enumerate(ENERGIES) if (m, rs, "d") in rows]     # energies this setup was run at
    ax.plot(E[have]*(1+0.03*list(MODELS).index(m)), [d_(m, ENERGIES[i][0]) for i in have], marker=MARKER[m], color=STYLE[m]["color"], ls=STYLE[m]["ls"], label=m)
ax.axhline(1, color="0.5", lw=0.8); ax.set_xscale("log"); ax.set_ylim(0.5, 1.4)
ax.set_xlabel(r"$\sqrt{s}$ [TeV]"); ax.set_ylabel(r"$\bar d/d$ (measured $p_T$ range)"); ax.legend(fontsize=8)
ax.set_title(r"$\bar d/d$, independent of $\sigma_0$ (cf. Fig. 6 of arXiv:1709.08522)")
plt.tight_layout(); plt.show()
for rs, _, _ in ENERGIES:
    print(f"{rs:>4s} TeV: dbar/d ALICE {d_('ALICE', rs):.2f}  " + "  ".join(f"{m} {d_(m, rs):.2f}" for m in MODELS if (m, rs, "d") in rows))
'''),
md(r'''
## 5. The coalescence parameter $B_2$

**What it is.** $B_2(p_T/A)=\dfrac{E_d\,d^3N_d/dp_d^3\ \text{at}\ p_{T,d}=2\,p_T/A}{\left(E_p\,d^3N_p/dp_p^3\ \text{at}\ p_{T,p}=p_T/A\right)^2}$: the deuteron yield divided by the square of the proton yield at half the deuteron's $p_T$ ($A=2$ nucleons). Deuterons go with protons, antideuterons with antiprotons. ALICE's values: HEPData Tables 6–8 (paper's Fig. 8).

**What the next cell does.** The same ratio for each generator, with the deuteron yield at its best-fit $1/\sigma_0$ (Sec. 2) and the generator's own (anti)protons (not from weak decays). Both yields per inelastic event carry a factor $f_{\rm ND}$, so $B_2$ gets $1/f_{\rm ND}$.

**What to take from it.** $B_2$ separates the formation step from the nucleon spectrum: a generator with the wrong protons but the right $B_2$ forms deuterons correctly from the protons it has.
'''),
code(r'''
fig, axes = plt.subplots(1, 3, figsize=(16, 4.5), sharey=True)
for ax, (rs, sqrts, table) in zip(axes, ENERGIES):
    fnd, _ = f_nd(sqrts); setup = f"ALICE_{rs}TeV"
    blocks = read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp_deuterons_HEPData_ins1625294_Table{5+table}.csv"))
    for (sp, psp, lab, c), blk in zip([("d", "p", "d", "black"), ("dbar", "pbar", r"\bar d", "0.55")], blocks):
        r_ = [r for r in blk if r[3] is not None]
        lo, hi, val = (np.array([r[i] for r in r_]) for i in (1, 2, 3))
        err = np.hypot(0.5*(np.array([r[4] for r in r_])-np.array([r[5] for r in r_])), 0.5*(np.array([r[6] for r in r_])-np.array([r[7] for r in r_])))
        ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), fmt="o" if sp == "d" else "o", mfc=c if sp == "d" else "white",
                    color="black", ms=4, label=f"ALICE ${lab}$")
        if sp != "d": continue
        for m in MODELS:
            if not available(m, setup): continue
            dv, _ = in_bins(deuteron_hist(m, setup, sp, "invyield_midy"), 2*lo, 2*hi)
            pv = at_points(f"nucleon_{psp}_invyield_midy", model_path(m, setup), 0.5*(lo+hi))   # protons at the bin centre
            ax.plot(0.5*(lo+hi), fits[(m, rs, sp)][0]*dv/(fnd*pv**2), marker=MARKER[m], color=STYLE[m]["color"], ls=STYLE[m]["ls"], label=f"{m} $d$")
    ax.set_title(f"pp {rs} TeV"); ax.set_xlabel(r"$p_T/A$ [GeV]"); ax.legend(fontsize=7)
axes[0].set_ylabel(r"$B_2$ [GeV$^2$]"); axes[0].set_ylim(0, None)
plt.tight_layout(); plt.show()
'''),
md(r'''
## 6. The protons: the paper's Fig. 8

**What it is.** The deuterons are built from the nucleons, so their spectrum cannot be described better than the nucleons'. ALICE measured the spectrum of **primary** protons + antiprotons at 7 TeV (Eur. Phys. J. C 75 (2015) 226, [arXiv:1504.00024](https://arxiv.org/abs/1504.00024), HEPData ins1357424, Table 1): $\frac{1}{N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ at $|y|<0.5$, "primary" meaning not from weak decays (e.g. $\Lambda\to p\pi^-$), the same definition as for the nucleons that enter our deuteron formation.

**What the next cell shows.** That spectrum as an invariant yield (divided by $2\pi p_T$, as in the paper's Fig. 8), ALICE against the three generators (×$f_{\rm ND}$), with the ratio generator/ALICE below.

**What to take from it.** Sherpa before had about twice ALICE's protons at $p_T\approx0.3$–0.5 GeV and too few at 1–2.5 GeV. Tuned Sherpa: 1.20× ALICE at 0.3–0.5 GeV (Pythia 1.33), 0.76× at 1–2 GeV (Pythia 0.85), right at 2.5–3 GeV (1.03), but ~45 % high above 3 GeV; that excess was already ~35 % before any tuning, and neither the AHADIC settings nor Amisic's `PT_Min` remove it (scans in `../deuteron_branch_changes.md`). Above 3 GeV Sherpa's curve also fluctuates (30–90 protons per bin with 100k events). These protons hardly matter for the ALICE deuterons, which are made from nucleons of about half the deuteron $p_T$.
'''),
code(r'''
rows_ = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, "ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv"))[2] if r[3] is not None]
alo, ahi, aval = (np.array([r[i] for r in rows_]) for i in (1, 2, 3)); aerr = np.hypot([r[4] for r in rows_], [r[6] for r in rows_])
apt = 0.5*(alo+ahi); ainv, ainverr = aval/(2*np.pi*apt), aerr/(2*np.pi*apt)
f7, _ = f_nd(7000.)
fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
ax.errorbar(apt, ainv, yerr=ainverr, xerr=0.5*(ahi-alo), label=r"ALICE $p+\bar p$", zorder=5, **dict(STYLE["data"], ms=3))
for m in MODELS:
    if not available(m, "ALICE_7TeV"): continue
    tot = None
    for sp in ("p", "pbar"):
        x, y, e, w = load_hist(f"nucleon_{sp}_invyield_midy", model_path(m, "ALICE_7TeV")); tot = y if tot is None else tot + y
    tot = f7*tot; ok = tot > 0
    ax.plot(x[ok], tot[ok], label=m, **STYLE[m])
    keep = ok & (x > alo[0]) & (x < ahi[-1]); rx.plot(x[keep], tot[keep]/np.interp(x[keep], apt, ainv), **STYLE[m])
ax.set_yscale("log"); ax.set_xlim(0, 4); ax.set_ylim(3e-5, 0.6); ax.legend(fontsize=8)
ax.set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]")
ax.set_title(r"primary $p+\bar p$, pp at $\sqrt{s}$ = 7 TeV, $|y|<0.5$ (cf. Fig. 8)")
rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0, 2.5); rx.set_ylabel("MC / ALICE"); rx.set_xlabel(r"$p_T$ [GeV]")
plt.tight_layout(); savefig(fig, "alice_fig8"); plt.show()
'''),md(r'''
## 7. Three tuning rounds after s3_053: the proton tail above 3 GeV

**Why.** s3_053 describes the protons well up to 3 GeV but is too high above: 1.09× ALICE at 3–3.5 GeV, 1.66× at 3.5–4 GeV (0.9M events). The decays are not the cause: the fractions of protons made directly (≈70 %) and from $\Delta$ decays (≈27 %) are the same in Sherpa and Pythia at every $p_T$, and both kinds are too hard by the same amount. The excess is in how hard the primary baryons are made.

**What was run.** Three rounds of random points at 7 TeV (100k events each), scored on π, K, p, the proton tail (3–5 GeV), the deuterons and $B_2$ (`../TuneScan/`):

| round | config | parameters varied | best point | result |
|---|---|---|---|---|
| 4 | `scan_config_fine.yaml` | the eight s3_053 parameters, ±10–20 % | s4_005 | lower tail, but fewer protons at 1–3 GeV; nothing beats s3_053 overall |
| 5 | `scan_config_leading.yaml` | `ALPHA_H`, `BETA_H`, `GAMMA_H`, `KT_0`, `PT_MAX` around s4_005 | s5_003 | the `_H` set hardly moves the tail; larger `KT_0`, `PT_MAX` lift the 1–3 GeV protons |
| 6 | `scan_config_kicks.yaml` | `KT_0` 2.4–3.2, `PT_MAX` 3.1–4.0, `SIGMA_ND_NORM`, `GAMMA_H` around s5_003 | s6_008, s6_010 | protons within 10 % of ALICE from 0.3 to 3.5 GeV, tail at 3.5–4 GeV reduced |

| | `SIGMA_ND_NORM` | `PT_0(ref)` | `ETA_Q` | `PT_MAX` | `KT_0` | `ALPHA_L` | `DECAY_THRESHOLD` | `BARYON_FRACTION` | `ALPHA_H`, `BETA_H`, `GAMMA_H` |
|---|---|---|---|---|---|---|---|---|---|
| s3_053 | 0.664 | 1.567 | 0.606 | 2.599 | 1.768 | 2.776 | 0.022 | 0.171 | −0.6, 1.8, 0.024 (default) |
| s4_005 | 0.749 | 1.612 | 0.568 | 2.944 | 1.770 | 2.976 | 0.048 | 0.163 | −0.6, 1.8, 0.024 (default) |
| s6_008 | 0.754 | 1.612 | 0.568 | 3.496 | 2.512 | 2.976 | 0.048 | 0.163 | 0.131, 1.321, 0.13 |
| s6_010 | 0.803 | 1.612 | 0.568 | 3.348 | 3.176 | 2.976 | 0.048 | 0.163 | 0.131, 1.321, 0.14 |

**What moves the tail** (linear fits over each round, change of MC/ALICE per one standard deviation of the parameter). `SIGMA_ND_NORM` sets its level (−0.09 to −0.10 at 3–4 GeV) but moves the 1–2 GeV protons the same way. `PT_MAX` raises both the 1–3 GeV protons and the tail. `KT_0` is the one parameter that tilts the spectrum: it raises the 1–3 GeV protons without raising the tail. The `_H` set, the z shape of cluster splittings that contain a quark of the perturbative stage (`Cluster_Splitter.C`), changes the tail by at most 0.03–0.05.

s6_008 and s6_010 are four seeds merged at 7 TeV (400k events), s3_053 and s4_005 one run of 100k: their curves fluctuate above 3 GeV (tens of protons per bin).
'''),
code(r'''
COMPARE = ["Sherpa s3_053", "Sherpa s4_005", "Sherpa s6_008", "Sherpa s6_010", "Pythia 8.318"]
setup7 = "ALICE_7TeV"; f7, _ = f_nd(7000.)
def mc_sum(m, names):
    "f_ND x sum of the histograms `names` of model m at 7 TeV, or None if the model has no such histograms"
    tot = None
    for n in names:
        if not os.path.exists(os.path.join(model_path(m, setup7), n + ".dat")): return None
        x, y, e, w = load_hist(n, model_path(m, setup7)); tot = y if tot is None else tot + y
    return x, f7*tot
def ratio_ranges(m, names, apt, ainv, ranges):
    c = mc_sum(m, names)
    if c is None: return ["  -- "]*len(ranges)
    x, t = c; out = []
    for a, b in ranges:
        ma, mm = (apt > a) & (apt < b), (x > a - 1e-9) & (x < b - 1e-9)
        out.append(f"{t[mm].mean()/ainv[ma].mean():5.2f}")
    return out

# protons, as Sec. 6
fig, (ax, rx) = plt.subplots(2, 1, figsize=(8, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
ax.errorbar(apt, ainv, yerr=ainverr, xerr=0.5*(ahi-alo), label=r"ALICE $p+\bar p$", zorder=5, **dict(STYLE["data"], ms=3))
rx.fill_between(apt, 1-ainverr/ainv, 1+ainverr/ainv, color="0.85", step="mid")
for m in COMPARE:
    c = mc_sum(m, ("nucleon_p_invyield_midy", "nucleon_pbar_invyield_midy"))
    if c is None: continue
    x, tot = c; ok = tot > 0
    ax.plot(x[ok], tot[ok], label=m, **STYLE[m])
    keep = ok & (x > alo[0]) & (x < ahi[-1]); rx.plot(x[keep], tot[keep]/np.interp(x[keep], apt, ainv), **STYLE[m])
ax.set_yscale("log"); ax.set_xlim(0, 4); ax.set_ylim(3e-5, 0.6); ax.legend(fontsize=8)
ax.set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]")
ax.set_title(r"primary $p+\bar p$, pp at $\sqrt{s}$ = 7 TeV, $|y|<0.5$: tuning rounds")
rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0, 2.5); rx.set_ylabel("MC / ALICE"); rx.set_xlabel(r"$p_T$ [GeV]")
plt.tight_layout(); savefig(fig, "alice_fig8_fine"); plt.show()
PRANGES = [(0.3, 0.6), (0.6, 1), (1, 2), (2, 3), (3, 3.5), (3.5, 4)]
print("protons MC/ALICE   " + "".join(f"{a:.1f}-{b:<5.1f}" for a, b in PRANGES))
for m in COMPARE:
    print(f"{m:18s} " + "    ".join(ratio_ranges(m, ("nucleon_p_invyield_midy", "nucleon_pbar_invyield_midy"), apt, ainv, PRANGES)))
'''),
md(r'''
**Pions and kaons.** The same parameters set the meson spectra, so a setting that fixes the proton tail must not spoil them. ALICE's $\pi^++\pi^-$ and $K^++K^-$ from the same paper (HEPData ins1357424, Table 1), with the ratio MC/ALICE. For Pythia, note that with lifetimes above 100 fm/c kept stable the $\eta$ and $\eta'$ do not decay, so the pions from their decays (which ALICE counts) are missing, a few per cent at low $p_T$.
'''),
code(r'''
blocks_ = read_hepdata_blocks(os.path.join(HEPDATA, "ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv"))
fig, axes = plt.subplots(2, 2, figsize=(14, 7), sharex="col", gridspec_kw={"height_ratios": [3, 1]})
for col, (label, block, names, ranges) in enumerate([
        (r"$\pi^++\pi^-$", 0, ("meson_pi_plus_invyield_midy", "meson_pi_minus_invyield_midy"), [(0.1, 0.5), (0.5, 1), (1, 2), (2, 3)]),
        (r"$K^++K^-$",     1, ("meson_K_plus_invyield_midy", "meson_K_minus_invyield_midy"),   [(0.2, 0.5), (0.5, 1), (1, 2), (2, 3), (3, 4)])]):
    r_ = [r for r in blocks_[block] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in r_]) for i in (1, 2, 3)); err = np.hypot([r[4] for r in r_], [r[6] for r in r_])
    pt = 0.5*(lo+hi); inv, inverr = val/(2*np.pi*pt), err/(2*np.pi*pt)
    ax, rx = axes[0, col], axes[1, col]
    ax.errorbar(pt, inv, yerr=inverr, xerr=0.5*(hi-lo), label="ALICE", zorder=5, **dict(STYLE["data"], ms=3))
    rx.fill_between(pt, 1-inverr/inv, 1+inverr/inv, color="0.85", step="mid")
    print(f"{label} MC/ALICE".ljust(26) + "".join(f"{a:.1f}-{b:<5.1f}" for a, b in ranges))
    for m in COMPARE:
        c = mc_sum(m, names)
        if c is None: continue
        x, tot = c; ok = tot > 0
        ax.plot(x[ok], tot[ok], label=m, **STYLE[m])
        keep = ok & (x > lo[0]) & (x < hi[-1]); rx.plot(x[keep], tot[keep]/np.interp(x[keep], pt, inv), **STYLE[m])
        print(f"  {m:22s}  " + "    ".join(ratio_ranges(m, names, pt, inv, ranges)))
    ax.set_yscale("log"); ax.set_xlim(0, hi[-1]); ax.legend(fontsize=8); ax.set_title(f"{label}, pp 7 TeV, $|y|<0.5$")
    rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0, 2); rx.set_xlabel(r"$p_T$ [GeV]")
axes[0, 0].set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]"); axes[1, 0].set_ylabel("MC / ALICE")
plt.tight_layout(); plt.show()
'''),
md(r'''
**Deuterons and $B_2$ at 7 TeV.** As Secs. 2 and 5: $1/\sigma_0$ fitted to each spectrum, $B_2$ with the generator's own protons at that $1/\sigma_0$.
'''),
code(r'''
table7 = [t for rs, _, t in ENERGIES if rs == "7"][0]
fits7 = {}
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
for col, (sp, block, lab) in enumerate([("d", 0, "d"), ("dbar", 1, r"\bar d")]):
    ax = axes[col]
    lo, hi, val, err = alice_spectrum(table7, block)
    ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE", zorder=5, **STYLE["data"])
    for m in COMPARE:
        if not available(m, setup7): continue
        h = deuteron_hist(m, setup7, sp, "invyield_midy")
        mc, _ = in_bins(h, lo, hi, scale=f7)
        x, xe, chi2 = fit_scale(val, err, mc); fits7[(m, sp)] = (x, chi2, len(val), x*mc/val, 0.5*(lo+hi))
        xs, ys, _, _ = h; ok = ys > 0
        ax.plot(xs[ok], x*f7*ys[ok], label=f"{m}: $1/\\sigma_0$ = {x:.2f} b$^{{-1}}$, $\\chi^2$ = {chi2:.1f}/{len(val)}", **STYLE[m])
    ax.set_yscale("log"); ax.set_xlim(0, 3.2); ax.set_ylim(1e-7, 1e-3); ax.set_xlabel(r"$p_T$ [GeV]")
    ax.set_title(f"${lab}$, pp 7 TeV, $|y|<0.5$"); ax.legend(fontsize=7, loc="lower left")
axes[0].set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]")
ax = axes[2]
blk = read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp_deuterons_HEPData_ins1625294_Table{5+table7}.csv"))[0]
r_ = [r for r in blk if r[3] is not None]
lo, hi, val = (np.array([r[i] for r in r_]) for i in (1, 2, 3))
err = np.hypot(0.5*(np.array([r[4] for r in r_])-np.array([r[5] for r in r_])), 0.5*(np.array([r[6] for r in r_])-np.array([r[7] for r in r_])))
ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE $d$", **STYLE["data"])
for m in COMPARE:
    if (m, "d") not in fits7: continue
    dv, _ = in_bins(deuteron_hist(m, setup7, "d", "invyield_midy"), 2*lo, 2*hi)
    pv = at_points("nucleon_p_invyield_midy", model_path(m, setup7), 0.5*(lo+hi))   # protons at the bin centre
    ax.plot(0.5*(lo+hi), fits7[(m, "d")][0]*dv/(f7*pv**2), marker=MARKER[m], color=STYLE[m]["color"], ls=STYLE[m]["ls"], label=m)
ax.set_title(r"$B_2$ ($d$), pp 7 TeV"); ax.set_xlabel(r"$p_T/A$ [GeV]"); ax.set_ylabel(r"$B_2$ [GeV$^2$]"); ax.set_ylim(0, None); ax.legend(fontsize=8)
plt.tight_layout(); savefig(fig, "alice_7TeV_fine"); plt.show()
print(f"{'':18s} {'sp':5s} {'1/sig0':>7s} {'chi2':>12s}   MC/ALICE at 0-1, 1-2, 2-3 GeV")
for m in COMPARE:
    for sp in ("d", "dbar"):
        if (m, sp) not in fits7: continue
        x, c2, n, r, pt7 = fits7[(m, sp)]
        print(f"{m:18s} {sp:5s} {x:7.2f} {c2:7.1f}/{n:<3d}   " + "  ".join(f"{r[(pt7 >= a) & (pt7 < a+1)].mean():.2f}" for a in (0, 1, 2)))
'''),
md(r'''
**What to take from it.** Larger transverse kicks in the cluster decays (`KT_0` 2.5–3.2 GeV against 1.77 in s3_053) bring the protons from 0.3 to 3.5 GeV within 10 % of ALICE (Pythia: 1.28 at 0.3–0.6 GeV, 0.84 at 1–2 GeV). The last bin, 3.5–4 GeV, is still high (s6_008 1.38, s6_010 1.32, against 1.56 for s3_053 and 1.14 for Pythia). Because the deuterons are made from these nucleons, their $\chi^2$ at 7 TeV drops by a third ($d$) and a half ($\bar d$) against s3_053, close to Pythia. What remains is the deuterons above 2 GeV (0.56–0.70 of ALICE, Pythia 0.71–0.77): above 2 GeV Sherpa's protons have fewer neutrons within $k<0.5$ GeV than Pythia's, a correlation of the nucleons that the cluster parameters scanned here do not change.
'''),
])

# =============================================================================
NOTEBOOKS["ALICE_2019_Multiplicity.ipynb"] = notebook([
md(r'''
# ALICE 2019: deuterons versus the number of particles in the event (pp, 7 TeV)

## The experiment in plain words

Collisions at the same energy differ a lot: some produce few particles, others many. ALICE sorted its 7 TeV pp collisions by **how many charged particles** they had and measured deuterons and protons separately in each **multiplicity class** (Phys. Lett. B 794 (2019) 50, [arXiv:1902.09290](https://arxiv.org/abs/1902.09290), HEPData [ins1721729](https://www.hepdata.net/record/ins1721729)).

* **How events are sorted.** By the signal in the forward V0 detectors (the "V0M estimator"), which count charged particles at $2.8\le\eta\le5.1$ and $-3.7\le\eta\le-1.7$, away from the central region where the deuterons are measured. Events are ranked and split into classes by percentile: 0–4.7 % (the highest multiplicities), 4.7–9.5 %, 9.5–19 %, 19–38 %, 38–100 %. Only **INEL>0** events are used: at least one charged particle at $|\eta|<1$ ($\eta$, the pseudorapidity, is the angle-based version of rapidity).
* **What characterises a class.** Its mean number of charged particles per unit of pseudorapidity at $|\eta|<0.5$, $\langle dN_{ch}/d\eta\rangle$ (ALICE's Table 1: 17.5 in the highest class, 3.3 in the lowest).
* **What is measured per class.** Deuteron spectra (HEPData Tables 1–5), the deuteron-to-proton ratio $2d/(p+\bar p)$ (Table 12) and $B_2$ (Tables 6–10).

**Why it matters.** In the cross-section model, a pair's formation probability does not depend on the rest of the event, so the number of deuterons grows like the number of pairs, roughly like the square of the multiplicity. ALICE sees $d/p$ grow more slowly and $B_2$ fall with multiplicity, which it interprets as the deuteron source getting larger: an effect of the **distance** between the nucleons, which this model does not use.

**Generators in this notebook.** Sherpa with the original run cards (`Sherpa.Deuterons.ALICE_7TeV.yaml`, Amisic with fixed $p_{T,0}=0.5$ GeV) and the tuned Sherpa (`...Tuned.yaml`), five runs of 80k events each. Pythia 8.318: four runs of 1M non-diffractive events (`DeuteronFormationFromXSecs/pythia/ALICE_7TeV_multiplicity/`, made by `../PythiaReference/pythia_deuterons.cc` with the same V0 counting); its deuterons are the formed ones, its $\eta,\eta'$ stay undecayed (a few per cent fewer charged particles).
'''),
SETUP,
md(r'''
## 1. Building the classes

**What the next cell does.** For each Sherpa setup it adds up the five runs' `multiplicity_*.dat` files. These have one row per value of `nV0` (number of charged particles in the V0 acceptance), with the number of INEL>0 events, their summed $dN_{ch}$ at $|\eta|<0.5$, the summed $dN/dy$ of $d$, $\bar d$, $p$, $\bar p$ and the summed spectra. Events are ranked by `nV0` and split at ALICE's percentiles (`wts` = fraction of each row's events in each class); per class it computes $\langle dN_{ch}/d\eta\rangle$, the yields and spectra.

**What to take from it.** The table compares $\langle dN_{ch}/d\eta\rangle$ per class with ALICE (Sherpa simulates only non-diffractive events, which have somewhat more particles than ALICE's INEL>0 sample). The original Sherpa is ~30 % too high in every class. The tuned Sherpa has the right average over all events (5.95 vs 5.96) but too narrow a spread: too few particles in the highest classes (14.1 vs 17.5) and too many in the lowest (3.9 vs 3.3). The multiplicity *distribution* is therefore something the tuning does not yet describe.
'''),
code(r'''
CLASSES = [("I+II", 0., 4.7, 17.47, 0.52), ("III", 4.7, 9.5, 13.50, 0.40), ("IV+V", 9.5, 19., 10.76, 0.30),
           ("VI+VII", 19., 38., 7.54, 0.23), ("VIII+IX+X", 38., 100., 3.30, 0.13)]   # arXiv:1902.09290, Table 1
SPECIES = ["d", "dbar", "p", "pbar"]
SETUPS = {"Sherpa original cards": MULT_DIR, "Sherpa tuned": os.path.join(TUNE_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa master": os.path.join(MASTER_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa master + our tune": os.path.join(MASTER_TUNED_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa retuned (new branch)": os.path.join(RETUNE_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa s3_053": os.path.join(SCAN_S3_053_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa s4_005": os.path.join(SCAN_S4_005_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa s6_008": os.path.join(SCAN_S6_008_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa s6_010": os.path.join(SCAN_S6_010_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa Lund": os.path.join(LUND_S6_008_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa Lund tuned": os.path.join(LUND_TUNE1_DIR, "ALICE_7TeV_multiplicity"),
          "Sherpa Lund tune2": os.path.join(LUND_TUNE2_DIR, "ALICE_7TeV_multiplicity"),
          "Pythia 8.318": os.path.join(PYTHIA_DIR, "ALICE_7TeV_multiplicity")}
SETUPS = {k: v for k, v in SETUPS.items() if k in SHOWN}   # only the shown setups (deuteron_common.SHOWN)
SETUP_STYLE = {"Sherpa original cards": STYLE["Sherpa before"], "Sherpa tuned": STYLE["Sherpa tuned"],
               "Sherpa master": STYLE["Sherpa master"], "Sherpa master + our tune": STYLE["Sherpa master + our tune"],
               "Sherpa retuned (new branch)": STYLE["Sherpa retuned (new branch)"], "Sherpa s3_053": STYLE["Sherpa s3_053"], "Sherpa s4_005": STYLE["Sherpa s4_005"],
               "Sherpa s6_008": STYLE["Sherpa s6_008"], "Sherpa s6_010": STYLE["Sherpa s6_010"], "Sherpa Lund": STYLE["Sherpa Lund"], "Sherpa Lund tuned": STYLE["Sherpa Lund tuned"], "Sherpa Lund tune2": STYLE["Sherpa Lund tune2"], "Pythia 8.318": STYLE["Pythia 8.318"]}
def load_rows(path): return {int(r[0]): np.array(r[1:]) for r in np.atleast_2d(np.loadtxt(path))}
def build_classes(base):
    seeds = sorted(d for d in glob.glob(os.path.join(base, "seed_*")) if os.path.isdir(d))   # not the .log files
    ev, yl, y2 = {}, {sp: {} for sp in SPECIES}, {sp: {} for sp in SPECIES}
    for sd in seeds:
        for k, v in load_rows(os.path.join(sd, "multiplicity_events.dat")).items(): ev[k] = ev.get(k, 0) + v
        for sp in SPECIES:
            for k, v in load_rows(os.path.join(sd, f"multiplicity_{sp}_invyield_midy.dat")).items(): yl[sp][k] = yl[sp].get(k, 0) + v
            for k, v in load_rows(os.path.join(sd, f"multiplicity_{sp}_invyield_midy_sq.dat")).items(): y2[sp][k] = y2[sp].get(k, 0) + v
    nv0 = np.array(sorted(ev)); nev = np.array([ev[k][0] for k in nv0])
    order = np.argsort(-nv0); cum_hi = np.cumsum(nev[order])/nev.sum()*100; cum_lo = cum_hi - nev[order]/nev.sum()*100
    out = []
    for name, plo, phi, anch, anche in CLASSES:
        w = np.zeros(len(nv0)); w[order] = np.clip(np.minimum(cum_hi, phi) - np.maximum(cum_lo, plo), 0, None)/np.maximum(cum_hi-cum_lo, 1e-12)
        N = (w*nev).sum()
        c = {"name": name, "N": N, "nch": (w*np.array([ev[k][1] for k in nv0])).sum()/N, "alice_nch": anch, "alice_nch_err": anche}
        for j, sp in enumerate(SPECIES):
            c[sp+"_dndy"] = (w*np.array([ev[k][2+j] for k in nv0])).sum()/N
            c[sp] = (w[:, None]*np.array([yl[sp][k] for k in nv0])).sum(axis=0)/N
            c[sp+"_err"] = np.sqrt((w[:, None]*np.array([y2[sp][k] for k in nv0])).sum(axis=0))/N
        out.append(c)
    return out, inv_sigma0(read_summary(seeds[0])[1]), nev.sum()
CLS = {}
for name, base in SETUPS.items():
    CLS[name] = build_classes(base)
    print(f"{name}: {CLS[name][2]:.0f} INEL>0 events")
print(f"\n{'class':10s} " + "".join(f"{n:>24s}" for n in SETUPS) + f"{'ALICE':>16s}")
for i, c in enumerate(CLASSES):
    print(f"{c[0]:10s} " + "".join(f"{CLS[n][0][i]['nch']:24.2f}" for n in SETUPS) + f"{c[3]:10.2f} +- {c[4]:.2f}")
'''),
md(r'''
## 2. Deuteron spectra per class

**What the next cell does.** Per class, the deuteron spectrum $\frac{1}{N_{ev}}\frac{d^2N}{dy\,dp_T}$ (HEPData Tables 1–5) against each Sherpa setup. **One** $1/\sigma_0$ per setup is fitted to all five classes together, so the change from class to class is a prediction. $\chi^2$ per class in the panel.
'''),
code(r'''
pt_c = np.arange(0.05, 5.0, 0.1)
def class_dndpt(c, sp, l, h):
    m = (pt_c > l) & (pt_c < h)
    return (2*np.pi*pt_c[m]*c[sp][m]*0.1).sum()/(h-l)
DATA = []
for ic in range(len(CLASSES)):
    b = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table{ic+1}.csv"))[0] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in b]) for i in (1, 2, 3))
    err = np.hypot(0.5*(np.array([r[4] for r in b])-np.array([r[5] for r in b])), 0.5*(np.array([r[6] for r in b])-np.array([r[7] for r in b])))
    DATA.append((lo, hi, val, err))
SCALE = {}
fig, axes = plt.subplots(1, 5, figsize=(20, 4.4), sharey=True)
for name, (cls, inv0, _) in CLS.items():
    mcs = [np.array([class_dndpt(c, "d", l, h) for l, h in zip(lo, hi)]) for c, (lo, hi, val, err) in zip(cls, DATA)]
    s, _, _ = fit_scale(np.concatenate([d[2] for d in DATA]), np.concatenate([d[3] for d in DATA]), np.concatenate(mcs))
    SCALE[name] = s
    print(f"{name}: global fit 1/sigma_0 = {inv0*s:.2f}/b")
    for ax, c, (lo, hi, val, err), mc in zip(axes, cls, DATA, mcs):
        chi2 = (((val - s*mc)/err)**2).sum()
        ax.plot(0.5*(lo+hi), s*mc, marker="s", ms=3, label=f"{name}, $\\chi^2$ = {chi2:.0f}/{len(val)}", **SETUP_STYLE[name])
for ax, c, (lo, hi, val, err) in zip(axes, CLASSES, DATA):
    ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE", **dict(STYLE["data"], ms=3))
    ax.set_yscale("log"); ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_title(f"class {c[0]}"); ax.legend(fontsize=7)
axes[0].set_ylabel(r"$\frac{1}{N_{ev}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-1}$]")
plt.tight_layout(); plt.show()
'''),
md(r'''
## 3. $d/p$ and $B_2$ versus multiplicity

**What the next cell shows.**
* Left: $2d/(p+\bar p)$, the ratio of the $dN/dy$ of deuterons and protons, versus $\langle dN_{ch}/d\eta\rangle$ (ALICE Table 12, paper Fig. 8), each Sherpa setup at its global $1/\sigma_0$.
* Middle: the same, each normalised to its lowest-multiplicity value: the trend alone, independent of $\sigma_0$.
* Right: $B_2$ at $0.7<p_T/A<0.8$ GeV, $\frac{(d+\bar d)/2}{((p+\bar p)/2)^2}$ (ALICE Tables 6–10, paper Fig. 6).

**What to take from it.** ALICE's $d/p$ rises by ×1.7 from the lowest to the highest class and its $B_2$ falls; a formation model based only on momenta gives a steeper rise and a constant $B_2$ (see the introduction).
'''),
code(r'''
t12 = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, "ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table12.csv"))[0] if r[3] is not None]
fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(18, 4.5))
a1.errorbar([r[0] for r in t12], [r[3] for r in t12], yerr=[np.hypot(r[4], r[6]) for r in t12], label="ALICE", **STYLE["data"])
a2.plot([r[0] for r in t12], np.array([r[3] for r in t12])/t12[0][3], "o", color="black", label="ALICE")
b2data = []
for ic in range(len(CLASSES)):
    r_ = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table{6+ic}.csv"))[0] if r[3] is not None]
    sel = [r for r in r_ if r[1] >= 0.7-1e-9 and r[2] <= 0.8+1e-9]
    if sel: b2data.append((CLASSES[ic][3], np.mean([r[3] for r in sel]), np.mean([np.hypot(r[4], r[6]) for r in sel])))
a3.errorbar([b[0] for b in b2data], [b[1] for b in b2data], yerr=[b[2] for b in b2data], label="ALICE", **STYLE["data"])
for name, (cls, inv0, _) in CLS.items():
    s = SCALE[name]; nch = np.array([c["nch"] for c in cls])
    dp = np.array([2*s*c["d_dndy"]/(c["p_dndy"]+c["pbar_dndy"]) for c in cls])
    a1.plot(nch, dp, marker="s", label=name, **SETUP_STYLE[name])
    a2.plot(nch, dp/dp[-1], marker="s", label=name, **SETUP_STYLE[name])
    md_, mp_ = (pt_c > 1.4) & (pt_c < 1.6), (pt_c > 0.7) & (pt_c < 0.8)
    b2 = [s*0.5*(c["d"][md_]+c["dbar"][md_]).mean()/(0.5*(c["p"][mp_]+c["pbar"][mp_]).mean())**2 for c in cls]
    a3.plot(nch, b2, marker="s", label=name, **SETUP_STYLE[name])
    print(f"{name}: d/p highest/lowest class = {dp[0]/dp[-1]:.2f};  B2 [1e-3 GeV^2] per class (highest first): " + " ".join(f"{b*1e3:.1f}" for b in b2))
for a, yl, t in [(a1, r"$2d/(p+\bar p)$", "d/p vs multiplicity (cf. Fig. 8)"), (a2, "relative to the lowest class", r"trend only (no $\sigma_0$)"),
                 (a3, r"$B_2$ at $0.7<p_T/A<0.8$ GeV [GeV$^2$]", "B2 vs multiplicity (cf. Fig. 6)")]:
    a.set_xlabel(r"$\langle dN_{ch}/d\eta\rangle_{|\eta|<0.5}$"); a.set_ylabel(yl); a.set_title(t); a.legend(fontsize=8)
a1.set_ylim(0, None); a3.set_ylim(0, None)
plt.tight_layout(); plt.show()
'''),
])

# =============================================================================
NOTEBOOKS["Coalescence_Models.ipynb"] = notebook([
md(r'''
# The coalescence models: HardDisk and Exponential vs the cross-section model

The cross-section model of Dal & Raklev is one way to decide which nucleon pairs become deuterons. The traditional alternative is **coalescence**: a proton and a neutron merge if they are close in momentum (and possibly in space). Two versions are implemented as options of the same Sherpa phase (`MODEL: HardDisk` or `Exponential`):
* **HardDisk**: $P=1$ if $k<p_{\rm coal}$ (and, optionally, the pair's impact parameter $b<b_{\rm coal}$), else 0: the classic model, Eq. (1) of the paper.
* **Exponential**: $P=0.75\,e^{-b^2/C_b}\,e^{-k^2/C_p}$ with $C_{p,b}=p^2_{\rm coal},b^2_{\rm coal}/(0.75\pi)$: a smooth version with the form of the Wigner-function weight of a Gaussian deuteron.

Here $k=|\vec p_p-\vec p_n|$ in the pair's rest frame, $b$ the distance of closest approach of the two nucleons' straight paths in that frame (from Sherpa's hadron positions). The deuteron is emitted with a photon ($pn\to d\gamma$ kinematics). The runs used $p_{\rm coal}=0.2$ GeV; `_k` = momentum condition only, `_kb` = also $b<3$ fm. Like $1/\sigma_0$, $p_{\rm coal}$ is fitted to data; the yield scales roughly as $p_0^3$ (paper Sec. IV A), so a fitted normalisation $s$ means $p_0\simeq0.2\,{\rm GeV}\times s^{1/3}$, to be compared with the coalescence column of the paper's Tables VII, VIII.

All runs here use the **original** Sherpa run cards (no colour reconnections, fixed Amisic $p_{T,0}$), the same events for all three models.
'''),
SETUP,
md(r'''
**What the next cell does.** Defines where each model's output is (`model_dir`), and prints the antideuteron rate in ALEPH's window at LEP for each model, with the $p_0$ that would reproduce ALEPH ($p_0=0.2\,{\rm GeV}\times(5.9\times10^{-6}/{\rm rate})^{1/3}$).
'''),
code(r'''
SETUP_DIR = {"LEP": DATA_DIR, "BaBar_10.58GeV": BABAR_DIR, "ALICE_0.9TeV": os.path.join(OLD_ALICE, "ALICE_0.9TeV"),
             "ALICE_2.76TeV": os.path.join(OLD_ALICE, "ALICE_2.76TeV"), "ALICE_7TeV": os.path.join(OLD_ALICE, "ALICE_7TeV")}
CMODELS = [("CrossSection", "cross section"), ("HardDisk_k", "HardDisk, k only"), ("Exponential_k", "Exponential, k only")]
CSTYLE = {"CrossSection": STYLE["Sherpa tuned"], "HardDisk_k": STYLE["Pythia 8.318"], "Exponential_k": dict(color="#009E73", ls=":", lw=1.8)}
def model_dir(model, setup): return SETUP_DIR[setup] if model == "CrossSection" else os.path.join(COAL_DIR, model, setup)
for model in ["CrossSection", "HardDisk_k", "Exponential_k", "HardDisk_kb", "Exponential_kb"]:
    s, _ = read_summary(model_dir(model, "LEP")); db = s.get("expected_dbar_ALEPH_window", (0, 0))[0]
    extra = f"   -> p0 = {0.2*(5.9e-6/db)**(1/3):.3f} GeV" if model != "CrossSection" and db > 0 else "   (1/sigma_0 = 1.8/b)"
    print(f"{model:16s}: dbar in ALEPH window per Z decay {db:.2e}{extra}")
print("paper, row 'dbar, LEP': p0 = 152 MeV (Herwig++), 183 MeV (Pythia 8)")
'''),
md(r'''
**What the next cell shows.** For each model, the momentum difference $k$ of the pairs that form deuterons (normalised), at the Z pole and in pp at 7 TeV.

**What to take from it.** The essential difference: coalescence uses pairs with $k\lesssim0.3$ GeV only, the cross-section model mainly pairs near the $\Delta$ resonance at $k\approx1$ GeV.
'''),
code(r'''
fig, axes = plt.subplots(1, 2, figsize=(13, 4.5))
for ax, setup in zip(axes, ["LEP", "ALICE_7TeV"]):
    for model, lab in CMODELS:
        d = model_dir(model, setup); tot = None
        for ch in CHANNELS:
            if not os.path.exists(os.path.join(d, f"expected_k_{ch}.dat")): continue
            x, y, e, w = load_hist(f"expected_k_{ch}", d); tot = y if tot is None else tot + y
        ax.plot(x, tot/(tot.sum()*w), label=lab, **CSTYLE[model])
    ax.set_xlabel(r"$k$ [GeV]"); ax.set_ylabel("deuteron-forming pairs, normalised"); ax.set_xlim(0, 3); ax.set_title(setup); ax.legend(fontsize=8)
plt.tight_layout(); plt.show()
'''),
md(r'''
**What the next cell shows.** For each ALICE energy, the deuteron spectra of the three models, each at its own best-fit normalisation, against ALICE, with $\chi^2$; the fitted parameter is $1/\sigma_0$ for the cross-section model and $p_0$ for the coalescence models.

**What to take from it.** The cross-section model describes the shape best at all energies; the coalescence spectra are too steep and noisy (few pairs at very small $k$).
'''),
code(r'''
fig, axes = plt.subplots(3, 2, figsize=(13, 13))
for row, (rs, sqrts, table) in enumerate(ENERGIES):
    fnd, _ = f_nd(sqrts)
    for col, (sp, block, lab) in enumerate([("d", 0, "d"), ("dbar", 1, r"\bar d")]):
        ax = axes[row, col]; lo, hi, val, err = alice_spectrum(table, block)
        ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE", zorder=5, **STYLE["data"])
        for model, name in CMODELS:
            d = model_dir(model, f"ALICE_{rs}TeV")
            mc, _ = in_bins(f"expected_{sp}_invyield_midy", lo, hi, d, scale=fnd)
            s, _, chi2 = fit_scale(val, err, mc)
            par = f"$1/\\sigma_0$ = {1.8*s:.2f} b$^{{-1}}$" if model == "CrossSection" else f"$p_0$ = {0.2*s**(1/3):.3f} GeV"
            ax.plot(0.5*(lo+hi), s*mc, marker="s", ms=3, label=f"{name}: {par}, $\\chi^2$ = {chi2:.0f}/{len(val)}", **CSTYLE[model])
        ax.set_yscale("log"); ax.set_xlabel(r"$p_T$ [GeV]"); ax.set_title(f"${lab}$, pp {rs} TeV"); ax.legend(fontsize=7)
        ax.set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]")
plt.tight_layout(); plt.show()
'''),
md(r'''
**What the next cell shows.** The effect of the spatial condition $b<3$ fm (printed: yield with and without it), and the distribution of the impact parameter $b$ of $pn$ pairs at LEP and at 7 TeV.

**What to take from it.** At LEP most pairs that are close in momentum are also close in space (70–90 % survive the cut); at 7 TeV only 20–30 %, because Sherpa spreads its hadrons over tens of fm in pp collisions, whereas measurements of particle correlations (femtoscopy) find a source of about 1 fm. Spatial models need realistic hadron positions first.
'''),
code(r'''
for setup in ["LEP", "ALICE_7TeV"]:
    for m in ["HardDisk", "Exponential"]:
        a = read_summary(model_dir(m+"_k", setup))[0]["expected_d_total"][0]; b = read_summary(model_dir(m+"_kb", setup))[0]["expected_d_total"][0]
        print(f"{setup:11s} {m:28s}: d per event, k only {a:.2e}, with b < 3 fm {b:.2e}  (kept {b/a:.0%})")
fig, ax = plt.subplots(figsize=(7, 4.5))
for setup, st in [("LEP", STYLE["Pythia 8.318"]), ("ALICE_7TeV", STYLE["Sherpa tuned"])]:
    x, y, e, w = load_hist("pairs_b_pn", model_dir("HardDisk_k", setup)); ax.plot(x, y/(y.sum()*w), label=setup, **st)
ax.axvline(3, color="k", ls=":", lw=1); ax.set_xlabel(r"impact parameter $b$ of (anti)$pn$ pairs [fm]"); ax.set_ylabel("normalised"); ax.legend()
plt.tight_layout(); plt.show()
'''),
])

# =============================================================================
NOTEBOOKS["Sherpa_vs_Pythia_Diagnostics.ipynb"] = notebook([
md(r'''
# Why Sherpa and Pythia differ: the diagnostics

The deuteron yield per event is the sum over nucleon pairs of $P(k)=\sigma(k)/\sigma_0$. The cross sections are the same in both generators (Sherpa's agree with the paper's formulas to $10^{-6}$), so any difference in the deuterons comes from **which nucleon pairs each generator produces**. This notebook collects the checks that located the problem:

1. **Pairs per $k$**: how many nucleon pairs each generator has at each momentum difference, as ratios Sherpa/Pythia, before and after the tuning.
2. **Where the protons come from**: directly from the hadronisation or from decays of heavier baryons ($\Delta$, ...), and how hard their spectra are.
3. **Hadron rescattering on/off**: rescattering is the only step in which the hadrons' positions and times could change their momenta.

The conclusion of these checks: the protons Sherpa produced directly in the hadronisation were too soft and too many at high energy, because Sherpa had no colour reconnections between its many partonic interactions. Switching them on (after fixing the bug described in `../colour_reconnection_fix.md`) and adjusting two hadronisation settings gave the tuned Sherpa.
'''),
SETUP,
md(r'''
## 1. Nucleon pairs per $k$

**What the next cell shows.** For the three ALICE energies, the ratio Sherpa/Pythia of the number of pairs per event per $k$, for $pn$, $pp$ and $nn$ pairs together, in five $k$ ranges, for Sherpa before and Sherpa tuned. The table also gives the total number of $pn$ pairs per event.

**What to take from it.** The $\Delta$ window ($0.7<k<1.2$ GeV) makes most deuterons. Sherpa before had 0.8× Pythia's pairs there at 0.9 TeV and 1.5× at 7 TeV, an excess growing with energy; the tuned Sherpa's ratio is nearly flat in $k$ (the same distribution as Pythia, just fewer pairs: 0.5× at 0.9 TeV, 0.7× at 7 TeV), which is what improved the spectra; the smaller number is absorbed by a larger $1/\sigma_0$.
'''),
code(r'''
RANGES = [(0, .3), (.3, .7), (.7, 1.2), (1.2, 2.), (2., 5.)]
def pairs_in(model, setup):
    tot = None
    for t in ("pn", "pp", "nn"):
        x, y, e, w = load_hist(f"pairs_k_{t}", model_path(model, setup)); tot = y if tot is None else tot + y
    return np.array([(tot[(x > a) & (x < b)]*w).sum() for a, b in RANGES])
print("ratio to Pythia 8.318 of NN + NbarNbar pairs per event, k ranges [GeV]: " + "  ".join(f"{a}-{b}" for a, b in RANGES))
fig, axes = plt.subplots(1, 3, figsize=(16, 4.2), sharey=True)
for ax, (rs, _, _) in zip(axes, ENERGIES):
    setup = f"ALICE_{rs}TeV"; py = pairs_in("Pythia 8.318", setup)
    for m in [m for m in SHOWN if m.startswith("Sherpa") and available(m, setup)]:
        r = pairs_in(m, setup)/py
        print(f"{rs:>4s} TeV {m:28s}: " + "  ".join(f"{v:5.2f}" for v in r) + f"   (pn pairs per event: {read_summary(model_path(m, setup))[0]['pairs_pn'][0]:.1f}, Pythia {read_summary(model_path('Pythia 8.318', setup))[0]['pairs_pn'][0]:.1f})")
        ax.plot([0.5*(a+b) for a, b in RANGES], r, marker=MARKER[m], label=m, **STYLE[m])
    ax.axhline(1, color="k", lw=0.8); ax.axvspan(0.7, 1.2, color="grey", alpha=0.15)
    ax.set_xlabel(r"$k$ [GeV]"); ax.set_title(f"pp {rs} TeV"); ax.legend(fontsize=8)
axes[0].set_ylabel("pairs / Pythia 8.318")
plt.tight_layout(); plt.show()
'''),
md(r'''
## 2. Where the protons come from

**What the next cell shows.** For mid-rapidity ($|y|<0.5$) protons and antiprotons, the share coming **directly** from the hadronisation and from **$\Delta(1232)$** decays (other sources are below 1 %), and the mean $p_T$ of the direct ones, for each generator at 0.9 and 7 TeV. Recorded by `ProtonOrigin` in `Deuteron_Formation.C` and by the same classification in the Pythia program: the immediate parent of each proton.

**What to take from it.** The $\Delta$ share is about the same in all generators (~25–32 %); what differed was the direct protons' mean $p_T$: Pythia's rises from 0.65 to 0.81 GeV between 0.9 and 7 TeV (as in data), Sherpa before stayed at 0.64. The rise in Pythia comes from colour reconnections between the many partonic interactions of a high-energy collision. The tuned Sherpa reaches Pythia's value at 7 TeV (0.82 GeV), but is already at 0.78 GeV at 0.9 TeV: its reconnections act too strongly at low energy. This is the remaining energy dependence of $1/\sigma_0$ (ALICE notebook, Sec. 3), and the natural next step is a reconnection strength that depends on the collision energy.
'''),
code(r'''
SRC = {"Sherpa before": lambda rs: os.path.join(ORIGIN_DIR, f"sherpa_ALICE_{rs}TeV"),
       "Sherpa tuned":  lambda rs: os.path.join(TUNE_DIR, f"ALICE_{rs}TeV"),
       "Sherpa master": lambda rs: os.path.join(MASTER_DIR, f"ALICE_{rs}TeV"),
       "Sherpa master + our tune": lambda rs: os.path.join(MASTER_TUNED_DIR, f"ALICE_{rs}TeV"),
       "Sherpa retuned (new branch)": lambda rs: os.path.join(RETUNE_DIR, f"ALICE_{rs}TeV"),
       "Sherpa s3_053": lambda rs: os.path.join(SCAN_S3_053_DIR, f"ALICE_{rs}TeV"),
       "Sherpa s4_005": lambda rs: os.path.join(SCAN_S4_005_DIR, f"ALICE_{rs}TeV"),
       "Sherpa s6_008": lambda rs: os.path.join(SCAN_S6_008_DIR, f"ALICE_{rs}TeV"),
       "Sherpa s6_010": lambda rs: os.path.join(SCAN_S6_010_DIR, f"ALICE_{rs}TeV"),
       "Sherpa Lund":   lambda rs: os.path.join(LUND_S6_008_DIR, f"ALICE_{rs}TeV"),
       "Sherpa Lund tuned": lambda rs: os.path.join(LUND_TUNE1_DIR, f"ALICE_{rs}TeV"),
       "Sherpa Lund tune2": lambda rs: os.path.join(LUND_TUNE2_DIR, f"ALICE_{rs}TeV"),
       "Pythia 8.318":  lambda rs: os.path.join(ORIGIN_DIR, f"pythia_ALICE_{rs}TeV")}
SRC = {k: v for k, v in SRC.items() if k in SHOWN}
print(f"{'':14s} {'':7s} {'p+pbar/event':>13s} {'direct':>8s} {'Delta':>7s} {'<pT> direct':>12s}")
for m, f in SRC.items():
    for rs in ("0.9", "7"):
        d = f(rs)
        if not os.path.exists(os.path.join(d, "proton_origin_direct_invyield_midy.dat")): continue
        s = read_summary(d)[0]
        tot = sum(s.get(f"proton_midy_from_{o}", (0, 0))[0] for o in ("direct", "Delta", "Nstar", "strangebaryon", "other"))
        x, y, e, w = load_hist("proton_origin_direct_invyield_midy", d); dn = 2*np.pi*x*y*w
        print(f"{m:28s} {rs+' TeV':7s} {tot:13.3f} {s['proton_midy_from_direct'][0]/tot:8.0%} {s['proton_midy_from_Delta'][0]/tot:7.0%} {(x*dn).sum()/dn.sum():9.2f} GeV")
'''),
md(r'''
## 3. Hadron rescattering on vs off

**What the next cell shows.** The ratio (rescattering off)/(on) of the numbers of protons, antiprotons, nucleon pairs and expected deuterons per event, at LEP and the three ALICE energies (original run cards; runs with `HADRON_RESCATTERING: {MODE: false}` in `norescattering/`).

**What to take from it.** All ratios are 1 within about 1 %: rescattering, and with it the hadrons' positions and formation times, does not affect the nucleons or deuterons. (With Sherpa's current hadron positions it happens about 0.6 times per 7 TeV event.)
'''),
code(r'''
keys = ["nucleon_p_total", "nucleon_pbar_total", "pairs_pn", "pairs_pp", "pairs_nn", "expected_d_total", "expected_dbar_total"]
setups = [("LEP", DATA_DIR), ("ALICE_0.9TeV", os.path.join(OLD_ALICE, "ALICE_0.9TeV")),
          ("ALICE_2.76TeV", os.path.join(OLD_ALICE, "ALICE_2.76TeV")), ("ALICE_7TeV", os.path.join(OLD_ALICE, "ALICE_7TeV"))]
print(f"{'off/on':22s}" + "".join(f"{s:>18s}" for s, _ in setups))
for k in keys:
    line = f"{k:22s}"
    for s, on_dir in setups:
        a, b = read_summary(on_dir)[0][k], read_summary(os.path.join(NORESC_DIR, s))[0][k]
        r = b[0]/a[0]; line += f"{r:10.3f} +- {r*np.hypot(a[1]/a[0], b[1]/b[0]):.3f}"
    print(line)
'''),
])

if __name__ == "__main__":
    # python3 make_notebooks.py [name.ipynb ...]: only the given notebooks (default: all)
    wanted = sys.argv[1:] or list(NOTEBOOKS)
    for name, nb in NOTEBOOKS.items():
        if name not in wanted: continue
        json.dump(nb, open(os.path.join(XSECS, name), "w"), indent=1, ensure_ascii=False)
        print("wrote", name)
