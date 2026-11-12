"""Where the deuteron-forming nucleon pairs come from: Sherpa settings against Pythia 8.318 at 7 TeV.

  python3 partner_origin.py

Reads the partner_* histograms written by Deuteron_Formation.C (FillPartners) and by
PythiaReference/pythia_deuterons.cc, both for protons and antiprotons at |y| < 0.5:
  partner rate      fraction of protons with a same-sign neutron within k < 0.5 GeV
  cone pT           summed charged pT within DeltaR < 0.5 of the proton (jet activity),
                    for all protons and for the protons with a partner
  pair origin       proton and nearest partner made directly / one in a hadron decay / both
Several runs of one setting (seeds) are added up.  All quantities are ratios within one
generator, so the normalisation of the runs does not matter."""
import glob, os
import numpy as np

from paths import DEUTERON_DATA
W = DEUTERON_DATA
SETUPS = {
    'Sherpa ScanTune (s3_053)': sorted(glob.glob(os.path.join(W, 'partners', 'scanTune_7TeV'))
                                       + glob.glob(os.path.join(W, 'partners', 'scanTune_seed*'))),
    'Sherpa retuned':           sorted(glob.glob(os.path.join(W, 'partners', 'retuned_seed*'))),
    'Pythia 8.318':             [os.path.join(W, 'pythia_partners', 'ALICE_7TeV')],
}
BINS = [(0, 0.6), (0.6, 1), (1, 1.6), (1.6, 2), (2, 3), (3, 5)]   # on the 0.2 GeV bin edges

def counts(dirs, name):
    "raw counts per proton-pT bin, added over the runs (histograms are per event per GeV)"
    tot = np.zeros(len(BINS))
    for d in dirs:
        f = os.path.join(d, name + '.dat')
        if not os.path.exists(f): continue
        nev = float(open(os.path.join(d, 'summary.dat')).readline().split('events = ')[1].split()[0])
        rows = [l.split() for l in open(f) if l.strip() and not l.startswith('#') and len(l.split()) in (2, 3)]
        a = np.array([[float(r[0]), float(r[1])] for r in rows])[:-1]       # drop overflow
        w = a[1, 0] - a[0, 0]
        tot += np.array([(a[:, 1][(a[:, 0] >= lo - 1e-9) & (a[:, 0] < hi - 1e-9)]).sum()*w*nev for lo, hi in BINS])
    return tot

for label, dirs in SETUPS.items():
    if not dirs: continue
    nev = sum(float(open(os.path.join(d, 'summary.dat')).readline().split('events = ')[1].split()[0]) for d in dirs)
    P, K = counts(dirs, 'partner_protons_pT'), counts(dirs, 'partner_protons_k05_pT')
    C, PC = counts(dirs, 'partner_protons_cone_pT'), counts(dirs, 'partner_pairs_cone_pT')
    O = {o: counts(dirs, f'partner_pairs_{o}_pT') for o in ('direct_direct', 'direct_decay', 'decay_decay')}
    print(f'\n=== {label}: {len(dirs)} run(s), {nev:.3g} events')
    print(f"{'proton pT':>10} {'protons':>8} {'partnered':>9} {'partner rate':>18} {'cone pT, all':>13} {'cone pT, partnered':>19} {'origin dir-dir/dir-dec/dec-dec':>31}")
    for i, (lo, hi) in enumerate(BINS):
        if P[i] <= 0: continue
        rate, err = K[i]/P[i], np.sqrt(max(K[i], 1.))/P[i]
        o = '/'.join(f'{O[k][i]/K[i]:.2f}' if K[i] > 0 else ' -- ' for k in ('direct_direct', 'direct_decay', 'decay_decay'))
        pc = f'{PC[i]/K[i]:.2f} GeV' if K[i] > 0 else '--'
        print(f'{lo:4.1f}-{hi:<4.1f} {P[i]:9.0f} {K[i]:9.0f}   {100*rate:6.3f} +- {100*err:5.3f} %  {C[i]/P[i]:7.2f} GeV  {pc:>14}      {o}')
