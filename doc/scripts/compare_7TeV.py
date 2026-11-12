"""Compare Sherpa settings and Pythia at 7 TeV with ALICE, as numbers: pi, K, p spectra
(MC / ALICE in pT ranges), the deuteron spectrum (fitted 1/sigma_0, chi2, MC / ALICE in
pT ranges) and the deuteron partners of the protons.

  python3 compare_7TeV.py

Each setup is one output directory (seeds merged with merge_seeds.py).  Uses the
helpers of deuteron_common.py, so the numbers are those of ALICE_2018_pp.ipynb."""
import os, sys
import numpy as np
from paths import XSECS
os.chdir(XSECS)    # the notebooks' folder, as when they run
from deuteron_common import (DATA_DIR, HEPDATA, PYTHIA_DIR, ENERGIES, read_hepdata_blocks, load_hist,
                             read_summary, inv_sigma0, alice_spectrum, f_nd, in_bins, fit_scale)

SETUPS = {
    's3_053':               os.path.join(DATA_DIR, 'scanTune_s3_053', 'ALICE_7TeV'),
    's4_005':               os.path.join(DATA_DIR, 'scanTune_s4_005', 'ALICE_7TeV'),
    's5_003 (400k)':        os.path.join(DATA_DIR, 's5_003_seeds', 'merged'),
    's6_008 (400k)':        os.path.join(DATA_DIR, 's6_008_seeds', 'merged'),
    's6_010 (400k)':        os.path.join(DATA_DIR, 's6_010_seeds', 'merged'),
    'Pythia 8.318':         os.path.join(DATA_DIR, 'pythia_anti', 'ALICE_7TeV'),   # 1M, seed 14, with meson and partner histograms
}
SETUPS = {k: d for k, d in SETUPS.items() if os.path.exists(os.path.join(d, 'summary.dat'))}
# better statistics for single quantities: the deuterons of the Pythia run of ALICE_2018_pp.ipynb
# (1/sigma_0 = 18/b in the run, so many more deuterons are formed), and ScanTune's 9 merged seeds
# for the partners (0.9M events)
DEUTERON_DIR = {'Pythia 8.318': os.path.join(PYTHIA_DIR, 'ALICE_7TeV')}
PARTNER_DIR  = {'s3_053': os.path.join(DATA_DIR, 'partners_anti', 'merged_scanTune')}
W = 16
f7, _ = f_nd(7000.)

def header(title):
    print(f'\n{title}\n{"":10s}' + ''.join(f'{k:>{W}s}' for k in SETUPS))

# pi, K, p: ALICE EPJC 75 (2015) 226, 1/N_inel d2N/dy dpT of the sum of both charges
blocks = read_hepdata_blocks(os.path.join(HEPDATA, 'ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv'))
SPECIES = [('pions',   0, ('meson_pi_plus', 'meson_pi_minus'), [(0.1, 0.5), (0.5, 1), (1, 2), (2, 3)]),
           ('kaons',   1, ('meson_K_plus', 'meson_K_minus'),   [(0.2, 0.5), (0.5, 1), (1, 2), (2, 3), (3, 4)]),
           ('protons', 2, ('nucleon_p', 'nucleon_pbar'),       [(0.3, 0.6), (0.6, 1), (1, 1.5), (1.5, 2), (2, 2.5), (2.5, 3), (3, 4)])]
for label, block, hists, ranges in SPECIES:
    rows = [r for r in blocks[block] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in rows]) for i in (1, 2, 3))
    pt = 0.5*(lo + hi); inv = val/(2*np.pi*pt)
    curves = {}
    for k, d in SETUPS.items():
        tot = None
        if not os.path.exists(os.path.join(d, hists[0] + '_invyield_midy.dat')): curves[k] = None; continue   # e.g. no pi/K in Pythia
        for h in hists:
            x, y, e, w = load_hist(h + '_invyield_midy', d); tot = y if tot is None else tot + y
        curves[k] = (x, f7*tot)
    header(f'{label}, MC / ALICE')
    for a, b in ranges:
        m = (pt > a) & (pt < b)
        if not m.any() or b > hi[-1] + 1e-9: continue
        line = f'{a:.1f}-{b:<6.1f}'
        for k, c in curves.items():
            if c is None: line += f'{"--":>{W}s}'; continue
            x, t = c
            mm = (x > a - 1e-9) & (x < b - 1e-9)
            line += f'{t[mm].mean()/inv[m].mean():{W}.2f}'
        print(line)

# deuterons and antideuterons at their best-fit 1/sigma_0 (ALICE PRC 97 (2018) 024615)
table = [t for rs, _, t in ENERGIES if rs == '7'][0]
for sp, block in (('d', 0), ('dbar', 1)):
    lo, hi, val, err = alice_spectrum(table, block)
    pt = 0.5*(lo + hi); res = {}
    for k, d in SETUPS.items():
        d = DEUTERON_DIR.get(k, d)
        prefix = 'formed' if 'Pythia' in k else 'expected'
        x, y, e, w = load_hist(f'{prefix}_{sp}_invyield_midy', d)
        s = inv_sigma0(read_summary(d)[1])
        mc, _ = in_bins((x, y/s, e/s, w), lo, hi, scale=f7)
        fit, fe, chi2 = fit_scale(val, err, mc)
        res[k] = (fit, chi2, fit*mc/val)
    header(f'{sp} spectrum: fitted 1/sigma_0 [1/b], chi2/{len(val)}, MC / ALICE')
    print(f'{"1/sig0":10s}' + ''.join(f'{r[0]:{W}.2f}' for r in res.values()))
    print(f'{"chi2":10s}' + ''.join(f'{r[1]:{W}.1f}' for r in res.values()))
    for a, b in ((0, 1), (1, 2), (2, 3)):
        m = (pt >= a) & (pt < b)
        print(f'{a}-{b} GeV'.ljust(10) + ''.join(f'{r[2][m].mean():{W}.2f}' for r in res.values()))

# deuteron partners: fraction of protons with a same-sign neutron within k < 0.5 GeV
def raw(d, name):
    rows = [l.split() for l in open(os.path.join(d, name + '.dat'))
            if l.strip() and not l.startswith('#') and len(l.split()) in (2, 3)]
    a = np.array([[float(r[0]), float(r[1])] for r in rows])[:-1]
    nev = float(open(os.path.join(d, 'summary.dat')).readline().split('events = ')[1].split()[0])
    return a[:, 0], a[:, 1]*(a[1, 0] - a[0, 0])*nev
header('partner rate (%): protons with a neutron within k < 0.5 GeV')
for a, b in ((0, 0.6), (0.6, 1), (1, 1.6), (1.6, 2), (2, 3)):
    line = f'{a:.1f}-{b:<6.1f}'
    for k, d in SETUPS.items():
        d = PARTNER_DIR.get(k, d)
        if not os.path.exists(os.path.join(d, 'partner_protons_k05_pT.dat')): line += f'{"--":>{W}s}'; continue
        x, P = raw(d, 'partner_protons_pT'); _, K = raw(d, 'partner_protons_k05_pT')
        m = (x >= a - 1e-9) & (x < b - 1e-9)
        line += f'{100*K[m].sum()/P[m].sum():{W}.2f}'
    print(line)
