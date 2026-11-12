"""Fit a smooth model to the scan scores and find the best setting (in the spirit of Professor).

  python3 fit_scan.py [config]          # default scan_config_round3.yaml

Each chi2 part of results (pions, kaons, protons, deuterons, b2) is fitted separately as
log(chi2) = linear + quadratic terms in the parameters (scaled to [-1, 1]); with N points and
P parameters the full quadratic needs (P+1)(P+2)/2 coefficients, so the cross terms are used
only if there are at least 1.5 times that many points, else the diagonal model.  A ridge term
keeps the fit stable.  The weighted sum of the fitted parts (weights of the config) is then
minimised inside the scanned ranges.  Writes the best setting to fit_<results>.yaml and a
config that runs it (scan_config_fit.yaml, RUNS_DIR runs_fit)."""
import csv, os, sys
import numpy as np
import yaml
from scipy.optimize import minimize

from paths import TUNE_SCAN
cfg_file = sys.argv[1] if len(sys.argv) > 1 else 'scan_config_round3.yaml'
cfg = yaml.safe_load(open(os.path.join(TUNE_SCAN, cfg_file)))
pars = list(cfg['PARAMETERS'])
lo = np.array([cfg['PARAMETERS'][p]['RANGE'][0] for p in pars])
hi = np.array([cfg['PARAMETERS'][p]['RANGE'][1] for p in pars])
weights = {f'chi2_{k}': v['WEIGHT'] for k, v in cfg['SCORE']['OBSERVABLES'].items()}
rows = list(csv.DictReader(open(os.path.join(TUNE_SCAN, cfg['RESULTS']))))
U = (np.array([[float(r[p]) for p in pars] for r in rows]) - lo)/(hi - lo)*2 - 1
P = len(pars); full = len(rows) >= 1.5*(P + 1)*(P + 2)/2

def features(u):
    cols = [np.ones(len(u))] + [u[:, i] for i in range(P)] + [u[:, i]**2 for i in range(P)]
    if full: cols += [u[:, i]*u[:, j] for i in range(P) for j in range(i + 1, P)]
    return np.stack(cols, axis=1)

F = features(U); LAM = 0.3; models = {}
print(f'{len(rows)} points, {P} parameters, {"full" if full else "diagonal"} quadratic ({F.shape[1]} coefficients)')
for c in weights:
    y = np.log(np.array([float(r[c]) for r in rows]))
    models[c] = np.linalg.solve(F.T@F + LAM*np.eye(F.shape[1]), F.T@y)
    idx = np.arange(len(y)); np.random.default_rng(0).shuffle(idx); q = []
    for k in range(5):                       # 5-fold check on points not used in the fit
        te = idx[k::5]; tr = np.setdiff1d(idx, te)
        w = np.linalg.solve(F[tr].T@F[tr] + LAM*np.eye(F.shape[1]), F[tr].T@y[tr])
        q += list(np.exp(F[te]@w - y[te]))
    q = np.array(q)
    print(f'  {c:15s} prediction/truth on unseen points: median {np.median(q):.2f}, 68% in {np.percentile(q, 16):.2f}-{np.percentile(q, 84):.2f}')

def total(u):
    f = features(u[None, :])
    return sum(w*np.exp(f@models[c])[0] for c, w in weights.items())

best = None; rng = np.random.default_rng(1)
for _ in range(300):
    res = minimize(total, rng.uniform(-1, 1, P), bounds=[(-1, 1)]*P, method='L-BFGS-B')
    if best is None or res.fun < best.fun: best = res
x = lo + (best.x + 1)/2*(hi - lo); f = features(best.x[None, :])
print(f'\npredicted best total {best.fun:.1f}  (best scanned: {min(float(r["TOTAL"]) for r in rows):.1f})')
for c in weights: print(f'  predicted {c}: {np.exp(f@models[c])[0]:.2f}')
print('setting:')
for p, v, l, h in zip(pars, x, lo, hi):
    print(f'  {p:36s} {v:7.3f}   (range {l} - {h}){"   <-- at the edge" if min(v - l, h - v) < 0.03*(h - l) else ""}')
best_values = {p: round(float(v), 3) for p, v in zip(pars, x)}
out = os.path.join(TUNE_SCAN, 'fit_' + os.path.splitext(cfg['RESULTS'])[0] + '.yaml')
yaml.safe_dump(best_values, open(out, 'w'), sort_keys=False)
fitcfg = dict(cfg); fitcfg.update(MODE='grid', RUNS_DIR='runs_fit_' + os.path.splitext(cfg['RESULTS'])[0],
                                  RESULTS='results_fit_' + cfg['RESULTS'], PARAMETERS={p: {'VALUES': [v]} for p, v in best_values.items()})
yaml.safe_dump(fitcfg, open(os.path.join(TUNE_SCAN, 'scan_config_fit_' + os.path.splitext(cfg['RESULTS'])[0] + '.yaml'), 'w'), sort_keys=False)
print(f'-> {out}; run it with: python3 tune_scan.py make/run/score scan_config_fit_{os.path.splitext(cfg["RESULTS"])[0]}.yaml')
