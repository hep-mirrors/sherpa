"""Scorecard: one generator setting against ALICE, next to Pythia 8.318, as in Dal & Raklev
(one generator setting for everything, 1/sigma_0 fitted separately to each spectrum).

  python3 scorecard.py <setup dir> [<setup dir> ...]

A setup dir holds ALICE_0.9TeV/, ALICE_2.76TeV/, ALICE_7TeV/ (Deuteron_Formation output of
one setting at the three energies), e.g. SimpleLEP/DeuteronFormationFromXSecs/scan_s3_053.
Missing energies are skipped.  Lower chi2 is better; a check mark means better than Pythia."""
import os, sys, warnings
import numpy as np
warnings.filterwarnings('ignore', category=RuntimeWarning)
import matplotlib; matplotlib.use('Agg')
from deuteron_common import (HEPDATA, PYTHIA_DIR, alice_spectrum, f_nd, fit_scale, in_bins, inv_sigma0,
                             load_hist, read_hepdata_blocks, read_summary)
ENERGIES = [('0.9', 900., 1), ('2.76', 2760., 2), ('7', 7000., 3)]
MUNC = 0.05

def prefix(d): return 'formed' if 'pythia' in os.path.abspath(d) else 'expected'

def deuterons(d, sqrts, table):
    "(1/sigma_0, chi2, n) for d and dbar; chi2 with the data errors only, as in the paper"
    fnd, _ = f_nd(sqrts); s = inv_sigma0(read_summary(d)[1]); out = {}
    for sp, block in (('d', 0), ('dbar', 1)):
        lo, hi, val, err = alice_spectrum(table, block)
        x, y, e, w = load_hist(f'{prefix(d)}_{sp}_invyield_midy', d)
        mc, _ = in_bins((x, y/s, e/s, w), lo, hi, scale=fnd)
        inv, _, chi2 = fit_scale(val, err, mc)
        out[sp] = (inv, chi2, len(val))
    return out

blocks = read_hepdata_blocks(os.path.join(HEPDATA, 'ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv'))
def alice_pikp(block):
    rows = [r for r in blocks[block] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in rows]) for i in (1, 2, 3))
    err = np.hypot(0.5*(np.array([r[4] for r in rows]) - np.array([r[5] for r in rows])),
                   0.5*(np.array([r[6] for r in rows]) - np.array([r[7] for r in rows])))
    sel = (lo >= 0.3 - 1e-9) & (hi <= 5.0 + 1e-9)
    return lo[sel], hi[sel], val[sel], err[sel]

def spectrum(d, names, block):
    "chi2 per point of a pi/K/p spectrum at 7 TeV (data errors + 5 %), MC/ALICE, bin lows"
    if not all(os.path.exists(os.path.join(d, n + '_invyield_midy.dat')) for n in names): return None
    fnd, _ = f_nd(7000.); lo, hi, val, err = alice_pikp(block); tot = None
    for n in names:
        x, y, e, w = load_hist(n + '_invyield_midy', d)
        v, ve = 2*np.pi*x*y*fnd, 2*np.pi*x*e*fnd
        tot = (x, v, ve, w) if tot is None else (x, tot[1] + v, np.hypot(tot[2], ve), w)
    mc, mce = in_bins(tot, lo, hi)
    return ((mc - val)**2/(err**2 + (MUNC*val)**2 + mce**2)).mean(), mc/val, lo

b2rows = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, 'ALICE_pp_deuterons_HEPData_ins1625294_Table8.csv'))[0] if r[3] is not None]
b2_pt = np.array([r[0] for r in b2rows]); b2_val = np.array([r[3] for r in b2rows])
b2_err = np.hypot(np.array([r[4] for r in b2rows]), np.array([r[6] for r in b2rows]))
def b2(d, inv_d):
    fnd, _ = f_nd(7000.); s = inv_sigma0(read_summary(d)[1])
    xd, yd, _, _ = load_hist(f'{prefix(d)}_d_invyield_midy', d); xp, yp, _, _ = load_hist('nucleon_p_invyield_midy', d)
    mc = (np.interp(2*b2_pt, xd, yd)/s*inv_d*fnd)/(np.interp(b2_pt, xp, yp)*fnd)**2; r = mc/b2_val
    return ((mc - b2_val)**2/(b2_err**2 + (MUNC*b2_val)**2)).mean(), r[(b2_pt > 1.2) & (b2_pt < 1.45)].mean()/r[b2_pt < 0.6].mean()

def card(setup):
    "all numbers of one setup -> dict"
    c = {}
    for rs, sqrts, table in ENERGIES:
        d = os.path.join(setup, f'ALICE_{rs}TeV')
        if os.path.exists(os.path.join(d, 'summary.dat')): c[f'deut_{rs}'] = deuterons(d, sqrts, table)
    d7 = os.path.join(setup, 'ALICE_7TeV')
    if os.path.exists(os.path.join(d7, 'summary.dat')):
        for sp, names, block in (('protons', ['nucleon_p', 'nucleon_pbar'], 2), ('pions', ['meson_pi_plus', 'meson_pi_minus'], 0),
                                 ('kaons', ['meson_K_plus', 'meson_K_minus'], 1)):
            c[sp] = spectrum(d7, names, block)
        c['b2'] = b2(d7, c['deut_7']['d'][0])
        mf = os.path.join(d7, 'multiplicity_events.dat')
        if os.path.exists(mf):
            m = np.loadtxt(mf, comments='#'); c['nch'] = m[:, 2].sum()/m[:, 1].sum()
    return c

def show(name, c, ref=None):
    ok = lambda mine, theirs: '' if ref is None else (' ✓' if mine < theirs else ' ✗')
    print(f'\n=== {name}')
    for rs, _, _ in ENERGIES:
        if f'deut_{rs}' not in c: print(f'  Fig. 7, {rs:>4s} TeV: (not run)'); continue
        parts = []
        for sp in ('d', 'dbar'):
            inv, chi2, n = c[f'deut_{rs}'][sp]
            parts.append(f'{sp:4s} chi2 {chi2:6.1f}/{n:<2d}{ok(chi2, ref[f"deut_{rs}"][sp][1]) if ref else ""}  1/sigma0 {inv:.2f}/b')
        print(f'  Fig. 7, {rs:>4s} TeV: ' + '   '.join(parts))
    if 'deut_0.9' in c and 'deut_7' in c:
        drift = c['deut_0.9']['d'][0]/c['deut_7']['d'][0]
        rdrift = ref['deut_0.9']['d'][0]/ref['deut_7']['d'][0] if ref else None
        print(f'  1/sigma0 drift 0.9/7 TeV (d): x{drift:.2f}' + ('' if ref is None else (' ✓' if abs(drift - 1) < abs(rdrift - 1) else ' ✗')) + '   (1 = same at all energies)')
    if c.get('protons'):
        chi2, r, lo = c['protons']; tail = r[(lo >= 3) & (lo < 4)].mean()/r[(lo >= 1) & (lo < 2)].mean()
        print(f'  Fig. 8, protons 7 TeV: chi2/pt {chi2:.2f}{ok(chi2, ref["protons"][0]) if ref else ""}   MC/ALICE 0.3-0.6 {r[lo < 0.6].mean():.2f}, 1-2 {r[(lo >= 1) & (lo < 2)].mean():.2f},'
              f' 2.2-2.8 {r[(lo >= 2.2) & (lo < 2.8)].mean():.2f}, 3-4 {r[(lo >= 3) & (lo < 4)].mean():.2f}   tail index {tail:.2f}')
    if c.get('b2'):
        print(f'  B2, 7 TeV: chi2/pt {c["b2"][0]:.2f}{ok(c["b2"][0], ref["b2"][0]) if ref else ""}   rise vs ALICE {c["b2"][1]:.2f} (1 = as in data)')
    extra = []
    for sp in ('pions', 'kaons'):
        if c.get(sp): extra.append(f'{sp} chi2/pt {c[sp][0]:.1f}')
    if 'nch' in c: extra.append(f'<dNch/deta> {c["nch"]:.2f}')
    if extra: print('  sanity: ' + ', '.join(extra))

pythia = card(PYTHIA_DIR)
show('Pythia 8.318 (reference)', pythia)
for setup in sys.argv[1:]:
    show(os.path.basename(os.path.normpath(setup)), card(setup), pythia)
