"""Search for a setting that beats Pythia 8.318 on the ALICE 7 TeV data.

  python3 fit_vs_pythia.py [config]        # default scan_config_round3.yaml

Same smooth fits of log(chi2) as fit_scan.py, but a different target: minimise the WORST
ratio Sherpa chi2 / Pythia chi2 over the observables Pythia is compared on (protons,
deuterons, B2), so that no observable is traded against another.  A value below 1 means
the predicted setting is better than Pythia on all of them.  Pions and kaons (no Pythia
reference in our runs) must stay below MESON_CAP chi2 per point.
Pythia's values are computed with exactly the scan's rules (tune_scan.py score)."""
import csv, os, sys
import numpy as np
import yaml
from scipy.optimize import minimize

PYTHIA = {'chi2_protons': 2.69, 'chi2_deuterons': (99.8 + 54.0)/41, 'chi2_b2': 1.85}
MESON_CAP = 6.0

from paths import TUNE_SCAN
cfg_file = sys.argv[1] if len(sys.argv) > 1 else 'scan_config_round3.yaml'
cfg = yaml.safe_load(open(os.path.join(TUNE_SCAN, cfg_file)))
pars = list(cfg['PARAMETERS'])
lo = np.array([cfg['PARAMETERS'][p]['RANGE'][0] for p in pars])
hi = np.array([cfg['PARAMETERS'][p]['RANGE'][1] for p in pars])
rows = list(csv.DictReader(open(os.path.join(TUNE_SCAN, cfg['RESULTS']))))
U = (np.array([[float(r[p]) for p in pars] for r in rows]) - lo)/(hi - lo)*2 - 1
P = len(pars)

def features(u):
    return np.stack([np.ones(len(u))] + [u[:, i] for i in range(P)] + [u[:, i]**2 for i in range(P)], axis=1)
F = features(U); LAM = 0.3; models = {}
for c in list(PYTHIA) + ['chi2_pions', 'chi2_kaons']:
    y = np.log(np.array([float(r[c]) for r in rows]))
    models[c] = np.linalg.solve(F.T@F + LAM*np.eye(F.shape[1]), F.T@y)
pred = lambda u, c: float(np.exp(features(u[None, :])@models[c])[0])

def worst(u):
    w = max(pred(u, c)/v for c, v in PYTHIA.items())
    pen = sum(max(0., pred(u, c) - MESON_CAP)**2 for c in ('chi2_pions', 'chi2_kaons'))
    return w + 0.1*pen

# the scanned points themselves, by the same measure
print('scanned points closest to beating Pythia (worst ratio Sherpa/Pythia over p, d, B2):')
meas = sorted(rows, key=lambda r: max(float(r[c])/v for c, v in PYTHIA.items()))
for r in meas[:5]:
    print(f"  {r['point']}: worst ratio {max(float(r[c])/v for c, v in PYTHIA.items()):.2f}  "
          f"(p {float(r['chi2_protons'])/PYTHIA['chi2_protons']:.2f}, d {float(r['chi2_deuterons'])/PYTHIA['chi2_deuterons']:.2f}, "
          f"B2 {float(r['chi2_b2'])/PYTHIA['chi2_b2']:.2f}; pi {r['chi2_pions']}, K {r['chi2_kaons']})")

best = None; rng = np.random.default_rng(2)
for _ in range(400):
    res = minimize(worst, rng.uniform(-1, 1, P), bounds=[(-1, 1)]*P, method='Powell')
    if best is None or res.fun < best.fun: best = res
x = lo + (best.x + 1)/2*(hi - lo)
print(f'\npredicted: worst ratio Sherpa/Pythia = {max(pred(best.x, c)/v for c, v in PYTHIA.items()):.2f}  (< 1 beats Pythia on all)')
for c, v in PYTHIA.items(): print(f'  {c:15s} {pred(best.x, c):6.2f}   Pythia {v:.2f}')
for c in ('chi2_pions', 'chi2_kaons'): print(f'  {c:15s} {pred(best.x, c):6.2f}   (cap {MESON_CAP})')
print('setting:')
for p, v, l, h in zip(pars, x, lo, hi):
    print(f'  {p:36s} {v:7.3f}   (range {l} - {h}){"   <-- at the edge" if min(v - l, h - v) < 0.03*(h - l) else ""}')
vals = {p: round(float(v), 3) for p, v in zip(pars, x)}
yaml.safe_dump(vals, open(os.path.join(TUNE_SCAN, 'fit_vs_pythia_best.yaml'), 'w'), sort_keys=False)
fc = dict(cfg); fc.update(MODE='grid', RUNS_DIR='runs_fit_vs_pythia', RESULTS='results_fit_vs_pythia.csv',
                          PARAMETERS={p: {'VALUES': [v]} for p, v in vals.items()})
yaml.safe_dump(fc, open(os.path.join(TUNE_SCAN, 'scan_config_fit_vs_pythia.yaml'), 'w'), sort_keys=False)
print('-> fit_vs_pythia_best.yaml; run it with: python3 tune_scan.py make/run/score scan_config_fit_vs_pythia.yaml')
