"""Are the deuteron-forming pairs made next to an antibaryon (B ... Bbar ... B)?  Sherpa against Pythia at 7 TeV.

  python3 antibaryon_near.py

Reads the partner_*anti* histograms written by Deuteron_Formation.C (FillPartners) and by
PythiaReference/pythia_deuterons.cc, for (anti)protons at |y| < 0.5:
  anti within 1 GeV   fraction of protons with a baryon of opposite baryon number within
                      k < 1 GeV, for all protons and for the protons with a partner
                      (same-sign neutron within k < 0.5 GeV); their ratio is the enhancement
  nearest-antibaryon k  for proton pT 1.5-5 GeV, the fraction of protons whose nearest
                      antibaryon lies below k = 0.5, 1, 2 GeV, all protons and partnered ones
In Pythia the high-pT pairs come from two nearby diquark breaks in one string
(string_pairs_* diagnostic); if Sherpa's large cluster kicks smear this, its partnered
protons should show a smaller enhancement than Pythia's.  Several runs of one setting
(seeds) are added up; all quantities are ratios within one generator."""
import glob, os
import numpy as np

from paths import DEUTERON_DATA
W = DEUTERON_DATA
SETUPS = {
    'Sherpa ScanTune (s3_053)':        sorted(glob.glob(os.path.join(W, 'partners_anti', 'scanTune_seed*'))),
    'Sherpa retuned':                  sorted(glob.glob(os.path.join(W, 'partners_anti', 'retuned_seed*'))),
    'Sherpa ScanTune, default kicks':  sorted(glob.glob(os.path.join(W, 'partners_anti', 'smallKicks_seed*'))),
    'Pythia 8.318':                    [os.path.join(W, 'pythia_anti', 'ALICE_7TeV')],
}
BINS = [(0, 0.6), (0.6, 1), (1, 1.6), (1.6, 2), (2, 3), (3, 5)]   # on the 0.2 GeV bin edges
KCUTS = [0.5, 1., 2.]

def events(d):
    return float(open(os.path.join(d, 'summary.dat')).readline().split('events = ')[1].split()[0])

def raw(d, name):
    "bin lower edges and raw counts (histograms are per event per GeV); the last bin is the overflow"
    rows = [l.split() for l in open(os.path.join(d, name + '.dat'))
            if l.strip() and not l.startswith('#') and len(l.split()) in (2, 3)]
    a = np.array([[float(r[0]), float(r[1])] for r in rows])
    return a[:-1, 0], a[:-1, 1]*(a[1, 0] - a[0, 0])*events(d), a[-1, 1]*(a[1, 0] - a[0, 0])*events(d)

def binned(dirs, name, bins):
    "raw counts in the given ranges of the histogram variable, added over the runs"
    tot = np.zeros(len(bins))
    for d in dirs:
        if not os.path.exists(os.path.join(d, name + '.dat')): continue
        x, n, over = raw(d, name)
        tot += np.array([n[(x >= lo - 1e-9) & (x < hi - 1e-9)].sum() for lo, hi in bins])
    return tot

for label, dirs in SETUPS.items():
    dirs = [d for d in dirs if os.path.exists(os.path.join(d, 'partner_protons_anti1_pT.dat'))]
    if not dirs: continue
    P,  K  = binned(dirs, 'partner_protons_pT', BINS),       binned(dirs, 'partner_protons_k05_pT', BINS)
    PA, KA = binned(dirs, 'partner_protons_anti1_pT', BINS), binned(dirs, 'partner_pairs_anti1_pT', BINS)
    print(f'\n=== {label}: {len(dirs)} run(s), {sum(events(d) for d in dirs):.3g} events')
    print(f"{'proton pT':>10} {'protons':>8} {'anti<1GeV':>10} {'partnered':>9} {'anti<1GeV':>16} {'enhancement':>12}")
    for i, (lo, hi) in enumerate(BINS):
        if P[i] <= 0 or K[i] <= 0: continue
        fa, fk = PA[i]/P[i], KA[i]/K[i]
        ek = np.sqrt(max(KA[i], 1.))/K[i]
        print(f'{lo:4.1f}-{hi:<4.1f} {P[i]:9.0f} {fa:9.3f}  {K[i]:9.0f}  {fk:6.3f} +- {ek:5.3f}  {fk/fa:8.2f}')
    # nearest-antibaryon k for proton pT 1.5-5 GeV
    # all protons at pT 1.5-5 GeV: every one enters the k histogram once, k > 5 GeV in the overflow
    NP, NK = [sum(sum(raw(d, n)[1]) + raw(d, n)[2] for d in dirs)
              for n in ('partner_protons_kanti_pT1.5-5', 'partner_pairs_kanti_pT1.5-5')]
    cuts = [(0, c) for c in KCUTS]
    CP, CK = binned(dirs, 'partner_protons_kanti_pT1.5-5', cuts), binned(dirs, 'partner_pairs_kanti_pT1.5-5', cuts)
    print(f'  proton pT 1.5-5 GeV, nearest antibaryon within k <', ', '.join(f'{c:g}' for c in KCUTS), 'GeV:')
    print(f'    all protons ({NP:.0f}):      ' + '  '.join(f'{c/NP:.3f}' for c in CP))
    print(f'    partnered protons ({NK:.0f}):  ' + '  '.join(f'{c/NK:.3f}' for c in CK))
