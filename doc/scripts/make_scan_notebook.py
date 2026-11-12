"""Writes TuneScan.ipynb: every scanned Sherpa setting compared with the ALICE 7 TeV data.
Run after each round:   python3 make_scan_notebook.py && jupyter nbconvert --to notebook --execute --inplace ../examples/RemnantTests/SimpleLEP/TuneScan/TuneScan.ipynb
The rounds are listed in ROUNDS in the first code cell of the notebook."""
import os

import nbformat as nbf

from paths import TUNE_SCAN

nb = nbf.v4.new_notebook()
md = lambda s: nbf.v4.new_markdown_cell(s.strip())
code = lambda s: nbf.v4.new_code_cell(s.strip())

nb.cells = [
md(r'''
# Sherpa parameter scan against ALICE, pp at 7 TeV

Every point of the parameter scan (`tune_scan.py`, configs `scan_config*.yaml`) compared with ALICE:
the identified-hadron spectra $\pi^\pm$, $K^\pm$, $p+\bar p$ at $|y|<0.5$ (HEPData ins1357424, Table 1) and the
(anti)deuteron spectra (ins1625294, Table 3). Each point is one Sherpa run of 100k non-diffractive events on the
`hadronicRescattering` branch; its settings are in `runs*/<point>/card.yaml`.

**Normalisation.** Sherpa generates non-diffractive events only; its spectra are multiplied by the non-diffractive
fraction $f_{ND}$ of inelastic events (ALICE, ins1181770), as in `XSecs/ALICE_2018_pp.ipynb`. The deuteron curves
use the $1/\sigma_0$ fitted to each spectrum (only the data errors in the fit, as in Dal & Raklev).

**Score.** As in `tune_scan.py score`: $\chi^2$ per data point for each observable (data errors plus a 5 % model
uncertainty for $\pi$, $K$, $p$; the notebook definition for deuterons), total = $\pi$ + $K$ + 2$p$ + 2$d$.
'''),
code(r'''
import os, sys, csv, glob, warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)   # ALICE bins beyond the 5 GeV end of the histograms
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
sys.path.insert(0, '../../../../scripts')   # doc/scripts: deuteron_common.py, paths.py
from deuteron_common import *            # helpers and data paths of the XSecs notebooks
apply_style()

# the rounds of the scan: label -> (results file, runs directory, short description)
ROUNDS = {
    'round 1': ('results.csv',        'runs',        'Amisic + CR + BARYON_FRACTION + PT_MAX'),
    'round 2': ('results_round2.csv', 'runs_round2', 'AHADIC: PT_MAX, BARYON_FRACTION, ALPHA/BETA/GAMMA_D, KT_0, Delta weight'),
    'round 3': ('results_round3.csv', 'runs_round3', 'all parameters together'),
}
# reference setups drawn in every figure
REFS = {'Sherpa retuned (new branch)': os.path.join(RETUNE_DIR, 'ALICE_7TeV'),
        'Sherpa tuned': os.path.join(TUNE_DIR, 'ALICE_7TeV'),
        'Pythia 8.318': os.path.join(PYTHIA_DIR, 'ALICE_7TeV')}
N_SHOW = 8                                 # best points of each round shown one by one
FND, _ = f_nd(7000.)
'''),
md(r'''
## 1. All scanned points, best first

Columns: the parameters of the point, $\chi^2$ per point of each observable, the total, $\langle dN_{ch}/d\eta\rangle$
($|\eta|<0.5$, INEL>0 non-diffractive events) and the proton MC/ALICE averaged in four $p_T$ ranges. The last
column, the *tail index* $p(3\text{–}4\,\mathrm{GeV})/p(1\text{–}2\,\mathrm{GeV})$, is 1 for ALICE (and $\approx$1 for Pythia).
'''),
code(r'''
POINTS = []      # (round, point, row, directory)
for rnd, (res, runs, desc) in ROUNDS.items():
    if not os.path.exists(res): continue
    for row in csv.DictReader(open(res)):
        POINTS.append((rnd, row['point'], row, os.path.join(runs, row['point'], 'out')))
POINTS.sort(key=lambda p: float(p[2]['TOTAL']))
skip = {'point', 'TOTAL'}
print(f"{'round':8s} {'point':8s} {'total':>6s} {'chi2 pi/K/p/d':>22s} {'Nch':>5s} {'p 0.3-0.6 / 1-2 / 2.2-2.8 / 3-4':>32s} {'tail':>5s}   parameters")
for rnd, name, r, d in POINTS:
    pars = ', '.join(f"{k.split(':')[-1]}={r[k]}" for k in r if ':' in k)
    tail = float(r['p_3.0-4.0'])/float(r['p_1.0-2.0'])
    print(f"{rnd:8s} {name:8s} {float(r['TOTAL']):6.1f} {r['chi2_pions']:>5s}/{r['chi2_kaons']:>5s}/{r['chi2_protons']:>5s}/{r['chi2_deuterons']:>5s}"
          f" {r['dNch_deta']:>5s} {r['p_0.3-0.6']:>7s} {r['p_1.0-2.0']:>7s} {r['p_2.2-2.8']:>7s} {r['p_3.0-4.0']:>7s} {tail:5.2f}   {pars}")
'''),
md(r'''
## 2. Overview: proton MC/ALICE of every point

Each thin line is one scanned point, coloured by its total score (dark = best); the thick lines are the references.
This shows the range the scanned parameters can reach: none of the points bends the curve down above 3 GeV
while keeping it at 1 around 1–2 GeV.
'''),
code(r'''
blocks = read_hepdata_blocks(os.path.join(HEPDATA, 'ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv'))
def alice_pikp(block):
    "ALICE pi, K, p (particles + antiparticles): pT bin edges, dN/dy dpT per inelastic event, total error"
    rows = [r for r in blocks[block] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in rows]) for i in (1, 2, 3))
    err = np.hypot(0.5*(np.array([r[4] for r in rows]) - np.array([r[5] for r in rows])),
                   0.5*(np.array([r[6] for r in rows]) - np.array([r[7] for r in rows])))
    return lo, hi, val, err
SPECIES = {'pions': (0, ['meson_pi_plus', 'meson_pi_minus'], r'$\pi^++\pi^-$'),
           'kaons': (1, ['meson_K_plus', 'meson_K_minus'], r'$K^++K^-$'),
           'protons': (2, ['nucleon_p', 'nucleon_pbar'], r'$p+\bar p$')}

def mc_dndpt(d, names):
    "Sherpa/Pythia: sum of invariant-yield histograms -> dN/dy dpT per inelastic event (centres, values, errors, width)"
    tot = None
    for n in names:
        x, y, e, w = load_hist(n + '_invyield_midy', d)
        v, ve = 2*np.pi*x*y*FND, 2*np.pi*x*e*FND
        tot = (x, v, ve, w) if tot is None else (x, tot[1] + v, np.hypot(tot[2], ve), w)
    return tot

def ratio_to_alice(d, species):
    "MC/ALICE in the ALICE bins, or None if the run has no histogram for this species"
    block, names, _ = SPECIES[species]
    if not all(os.path.exists(os.path.join(d, n + '_invyield_midy.dat')) for n in names): return None
    lo, hi, val, err = alice_pikp(block)
    mc, _ = in_bins(mc_dndpt(d, names), lo, hi)
    return 0.5*(lo + hi), mc/val, err/val

def deuteron_fit(d, sp, block):
    "MC deuteron spectrum per unit 1/sigma_0, fitted to ALICE 7 TeV: centres, curve, 1/sigma_0, chi2, n"
    prefix = 'formed' if 'pythia' in d else 'expected'
    lo, hi, val, err = alice_spectrum(3, block)
    x, y, e, w = load_hist(f'{prefix}_{sp}_invyield_midy', d)
    s = inv_sigma0(read_summary(d)[1])
    mc, _ = in_bins((x, y/s, e/s, w), lo, hi, scale=FND)
    inv, _, chi2 = fit_scale(val, err, mc)
    return x, inv*FND*y/s, inv, chi2, len(val)

REF_STYLE = {'Sherpa retuned (new branch)': STYLE['Sherpa retuned (new branch)'],
             'Sherpa tuned': STYLE['Sherpa tuned'], 'Pythia 8.318': STYLE['Pythia 8.318']}

fig, ax = plt.subplots(figsize=(9, 4.8))
tots = np.array([float(p[2]['TOTAL']) for p in POINTS])
cmap = plt.get_cmap('viridis')
for (rnd, name, r, d), t in zip(POINTS[::-1], tots[::-1]):
    rr = ratio_to_alice(d, 'protons')
    if rr is None: continue
    ax.plot(rr[0], rr[1], color=cmap((t - tots.min())/(tots.max() - tots.min() + 1e-9)), lw=0.8, alpha=0.6)
for lab, d in REFS.items():
    rr = ratio_to_alice(d, 'protons')
    if rr is not None: ax.plot(rr[0], rr[1], label=lab, **REF_STYLE[lab])
ax.axhline(1, color='k', lw=0.8); ax.set_xlim(0, 5); ax.set_ylim(0, 3)
ax.set_xlabel(r'$p_T$ [GeV]'); ax.set_ylabel(r'$p+\bar p$: MC / ALICE')
ax.set_title(f'all {len(POINTS)} scanned points (colour: total score, dark = best)'); ax.legend(fontsize=9)
plt.tight_layout(); plt.show()
'''),
md(r'''
## 3. The best points, one by one

For each round, the `N_SHOW` best points. Each figure: **left**, the proton spectrum against ALICE with MC/ALICE
below; **middle**, MC/ALICE for pions, kaons and protons together; **right**, the deuteron and antideuteron
spectra with their fitted $1/\sigma_0$ and $\chi^2$ (same fit as `ALICE_2018_pp`). The references are drawn in
every panel: the retuned new-branch setting, the old-branch tune and Pythia.
'''),
code(r'''
def plot_point(rnd, name, r, d):
    fig = plt.figure(figsize=(17, 5.2))
    gs = GridSpec(2, 3, height_ratios=[3, 1.2], hspace=0.05, wspace=0.28)
    title = f"{rnd}, {name}: total {float(r['TOTAL']):.1f}  (" + ', '.join(f"{k.split(':')[-1]} = {r[k]}" for k in r if ':' in k) + ")"
    fig.suptitle(title, fontsize=11)
    curves = [(f'{name}', d, dict(color='black', ls='-', lw=1.8))] + [(lab, dd, REF_STYLE[lab]) for lab, dd in REFS.items()]

    # left: protons, spectrum and ratio
    a0, a1 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
    lo, hi, val, err = alice_pikp(2); c = 0.5*(lo + hi)
    a0.errorbar(c, val/(2*np.pi*c), yerr=err/(2*np.pi*c), xerr=0.5*(hi - lo), label=r'ALICE $p+\bar p$', **STYLE['data'])
    for lab, dd, st in curves:
        x, v, _, _ = mc_dndpt(dd, SPECIES['protons'][1]); ok = v > 0
        a0.plot(x[ok], v[ok]/(2*np.pi*x[ok]), label=lab, **st)
        rr = ratio_to_alice(dd, 'protons'); a1.plot(rr[0], rr[1], **st)
    a1.fill_between(c, 1 - err/val, 1 + err/val, color='0.85', zorder=0)
    a0.set_yscale('log'); a0.set_xlim(0, 5); a1.set_xlim(0, 5); a1.set_ylim(0, 3); a1.axhline(1, color='k', lw=0.8)
    a0.set_xticklabels([]); a0.set_ylabel(r'$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]')
    a1.set_xlabel(r'$p_T$ [GeV]'); a1.set_ylabel('MC / ALICE'); a0.legend(fontsize=8)
    a0.set_title(fr"protons: $\chi^2$/pt = {r['chi2_protons']}")

    # middle: pi, K, p ratios of this point (solid) and of the retuned reference (dashed)
    b = fig.add_subplot(gs[:, 1])
    for (sp, (_, _, lab)), col in zip(SPECIES.items(), ['#0072B2', '#009E73', '#D55E00']):
        rr = ratio_to_alice(d, sp); b.plot(rr[0], rr[1], color=col, lw=1.8, label=f"{lab}: $\\chi^2$/pt = {r['chi2_' + sp]}")
        rr = ratio_to_alice(REFS['Sherpa retuned (new branch)'], sp)
        if rr is not None: b.plot(rr[0], rr[1], color=col, lw=1.0, ls='--')
    b.axhline(1, color='k', lw=0.8); b.set_xlim(0, 5); b.set_ylim(0, 3); b.set_xlabel(r'$p_T$ [GeV]')
    b.set_ylabel('MC / ALICE'); b.legend(fontsize=8, title='solid: this point, dashed: retuned', title_fontsize=8)
    b.set_title(fr"$\langle dN_{{ch}}/d\eta\rangle$ = {r['dNch_deta']}")

    # right: deuterons and antideuterons
    e = fig.add_subplot(gs[:, 2])
    for (sp, block, lab, mk) in (('d', 0, r'$d$', 'o'), ('dbar', 1, r'$\bar d$', 's')):
        lo, hi, val, err = alice_spectrum(3, block); scale = 1 if sp == 'd' else 0.1
        e.errorbar(0.5*(lo + hi), scale*val, yerr=scale*err, xerr=0.5*(hi - lo), color='black', fmt=mk, ms=4,
                   label=f'ALICE {lab}' + ('' if sp == 'd' else r' $\times 0.1$'))
        for i, (labc, dd, st) in enumerate(curves):
            x, yfit, inv, chi2, n = deuteron_fit(dd, sp, block); ok = yfit > 0
            e.plot(x[ok], scale*yfit[ok], label=(f'{labc}: ' if sp == 'd' else '') + f'{lab} $1/\\sigma_0$ = {inv:.2f}/b, $\\chi^2$ = {chi2:.0f}/{n}'
                   if i == 0 or sp == 'd' else None, **st)
    e.set_yscale('log'); e.set_xlim(0, 3.5); e.set_ylim(1e-8, 1e-3); e.set_xlabel(r'$p_T$ [GeV]')
    e.set_ylabel(r'$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]'); e.legend(fontsize=7, loc='lower left')
    e.set_title(fr"deuterons: $\chi^2$/pt = {r['chi2_deuterons']}")
    plt.show()

for rnd in ROUNDS:
    best = [p for p in POINTS if p[0] == rnd][:N_SHOW]
    for p in best: plot_point(*p)
'''),
]
nbf.write(nb, os.path.join(TUNE_SCAN, 'TuneScan.ipynb'))
print('wrote TuneScan.ipynb')
