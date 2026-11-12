"""Parameter scan of Sherpa minimum-bias settings against ALICE data at 7 TeV.

  python3 tune_scan.py make   [config]   write <RUNS_DIR>/<point>/card.yaml for every point
  python3 tune_scan.py run    [config]   run every point that has not finished yet, PARALLEL at a time
  python3 tune_scan.py score  [config]   score every finished point -> <RESULTS> (best first)
  python3 tune_scan.py status [config]   show which points are finished
[config] defaults to scan_config.yaml (round 1); each round has its own
RUNS_DIR and RESULTS, set in its config (configs, base cards, runs and
results are in SimpleLEP/TuneScan).

Each point is a copy of BASE_CARD with the scanned parameters set; its
Deuteron_Formation output goes to runs/<point>/out.  The score compares the
pi, K, p spectra with ALICE (ins1357424) and the deuteron spectra with ALICE
(ins1625294, 1/sigma_0 fitted as in XSecs/ALICE_2018_pp.ipynb)."""
import csv, os, subprocess, sys, time
import numpy as np
import yaml

from paths import TUNE_SCAN, XSECS
CONFIG_FILE = sys.argv[2] if len(sys.argv) > 2 else 'scan_config.yaml'
CONFIG = yaml.safe_load(open(os.path.join(TUNE_SCAN, CONFIG_FILE)))
RUNS = os.path.join(TUNE_SCAN, CONFIG.get('RUNS_DIR', 'runs'))


# --------------------------------------------------------------------------
# the points
# --------------------------------------------------------------------------
def points():
    """[(name, {parameter: value})] in a fixed order; adding points keeps the old ones."""
    pars = CONFIG['PARAMETERS']
    if CONFIG['MODE'] == 'grid':
        grids = np.meshgrid(*[pars[p]['VALUES'] for p in pars], indexing='ij')
        rows = np.stack([g.ravel() for g in grids], axis=1)
    else:
        rng = np.random.default_rng(CONFIG['SEED'])
        lo = np.array([pars[p]['RANGE'][0] for p in pars])
        hi = np.array([pars[p]['RANGE'][1] for p in pars])
        rows = lo + (hi - lo)*rng.random((CONFIG['N_POINTS'], len(pars)))
    prefix = 'g' if CONFIG['MODE'] == 'grid' else f"s{CONFIG['SEED']}_"
    # EXTRA_POINTS: named points with given values (e.g. a reference setting), run and scored as the others
    extra = [(name, {p: float(v) for p, v in vals.items()}) for name, vals in (CONFIG.get('EXTRA_POINTS') or {}).items()]
    return extra + [(f'{prefix}{i+1:03d}', {p: round(float(v), 3) for p, v in zip(pars, row)})
                    for i, row in enumerate(rows)]

def set_value(card, path, value):
    "set card['A']['B']... = value for path 'A:B:...'"
    keys = path.split(':')
    node = card
    for k in keys[:-1]:
        node = node.setdefault(k, {})
    node[keys[-1]] = value

def out_dir(name):
    return os.path.join(RUNS, name, 'out')

def finished(name):
    "the point's run has written all its events"
    summary = os.path.join(out_dir(name), 'summary.dat')
    if not os.path.exists(summary): return False
    header = open(summary).readline()
    want = float(str(CONFIG['EVENTS']).replace('k', 'e3').replace('M', 'e6'))
    return 'events = ' in header and float(header.split('events = ')[1].split()[0]) >= want


# --------------------------------------------------------------------------
# make and run
# --------------------------------------------------------------------------
def make():
    base = open(os.path.join(TUNE_SCAN, CONFIG['BASE_CARD'])).read()
    for name, values in points():
        card = yaml.safe_load(base)
        for path, value in values.items():
            set_value(card, path, value)
        card['EVENTS'] = CONFIG['EVENTS']
        set_value(card, 'DEUTERON_FORMATION:OUTPUT', out_dir(name))
        os.makedirs(os.path.join(RUNS, name), exist_ok=True)
        with open(os.path.join(RUNS, name, 'card.yaml'), 'w') as f:
            f.write(f'# TuneScan point {name}, made from {CONFIG["BASE_CARD"]} by tune_scan.py\n')
            for path, value in values.items():
                f.write(f'#   {path} = {value}\n')
            yaml.safe_dump(card, f, sort_keys=False)
    print(f'{len(points())} cards in {RUNS}')

def run():
    todo = [name for name, _ in points() if not finished(name)]
    print(f'{len(todo)} points to run, {CONFIG["PARALLEL"]} at a time')
    running = {}
    while todo or running:
        for name, proc in list(running.items()):
            if proc.poll() is not None:
                del running[name]
                print(time.strftime('%H:%M'), name, 'finished' if finished(name) else f'FAILED (exit {proc.returncode})', flush=True)
        while todo and len(running) < CONFIG['PARALLEL']:
            name = todo.pop(0)
            d = os.path.join(RUNS, name)
            log = open(os.path.join(d, 'run.log'), 'w')
            running[name] = subprocess.Popen(['caffeinate', '-i', CONFIG['SHERPA'], '-f', 'card.yaml'],
                                             cwd=d, stdout=log, stderr=subprocess.STDOUT)
            print(time.strftime('%H:%M'), name, 'started', flush=True)
        time.sleep(20)

def status():
    for name, values in points():
        print(name, 'done' if finished(name) else '-', values)


# --------------------------------------------------------------------------
# score
# --------------------------------------------------------------------------
def score():
    os.chdir(XSECS)
    import matplotlib; matplotlib.use('Agg')
    from deuteron_common import (HEPDATA, alice_spectrum, f_nd, fit_scale, in_bins,
                                 inv_sigma0, load_hist, read_hepdata_blocks, read_summary)
    fnd, _ = f_nd(float(CONFIG['SQRTS']))
    sc = CONFIG['SCORE']; obs = sc['OBSERVABLES']; munc = sc['MODEL_UNCERTAINTY']

    # ALICE pi, K, p (sum of particles and antiparticles), dN/dy dpT per inelastic event
    blocks = read_hepdata_blocks(os.path.join(HEPDATA, 'ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv'))
    def alice_pikp(block):
        rows = [r for r in blocks[block] if r[3] is not None]
        lo, hi, val = (np.array([r[i] for r in rows]) for i in (1, 2, 3))
        err = np.hypot(0.5*(np.array([r[4] for r in rows]) - np.array([r[5] for r in rows])),
                       0.5*(np.array([r[6] for r in rows]) - np.array([r[7] for r in rows])))
        return lo, hi, val, err
    SPECIES = {'pions': (0, ['meson_pi_plus', 'meson_pi_minus']),
               'kaons': (1, ['meson_K_plus', 'meson_K_minus']),
               'protons': (2, ['nucleon_p', 'nucleon_pbar'])}

    def mc_dndpt(d, names):
        "sum of invariant-yield histograms -> dN/dy dpT per inelastic event"
        tot = None
        for n in names:
            x, y, e, w = load_hist(n + '_invyield_midy', d)
            v, ve = 2*np.pi*x*y*fnd, 2*np.pi*x*e*fnd
            tot = (x, v, ve, w) if tot is None else (x, tot[1] + v, np.hypot(tot[2], ve), w)
        return tot

    def chi2_spectrum(d, species, pt_range=None):
        block, names = SPECIES[species]
        lo, hi, val, err = alice_pikp(block)
        a, b = pt_range or obs[species]['PT_RANGE']
        sel = (lo >= a - 1e-9) & (hi <= b + 1e-9)
        lo, hi, val, err = lo[sel], hi[sel], val[sel], err[sel]
        mc, mce = in_bins(mc_dndpt(d, names), lo, hi)
        r = (mc - val)**2/(err**2 + (munc*val)**2 + mce**2)
        return r.sum()/len(r), mc/val

    def chi2_deuterons(d):
        "as in ALICE_2018_pp: 1/sigma_0 fitted to d and dbar separately, chi2 summed"
        s = inv_sigma0(read_summary(d)[1]); c2, n, fits = 0., 0, []
        for sp, block in (('d', 0), ('dbar', 1)):
            lo, hi, val, err = alice_spectrum(3, block)
            x, y, e, w = load_hist(f'expected_{sp}_invyield_midy', d)
            mc, _ = in_bins((x, y/s, e/s, w), lo, hi, scale=fnd)
            inv, _, chi2 = fit_scale(val, err, mc)
            c2 += chi2; n += len(val); fits.append(inv)
        return c2/n, fits

    # ALICE B2 of deuterons at 7 TeV (ins1625294, Table 8, block 0), per nucleon pT
    b2rows = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, 'ALICE_pp_deuterons_HEPData_ins1625294_Table8.csv'))[0]
              if r[3] is not None]
    b2_pt = np.array([r[0] for r in b2rows]); b2_val = np.array([r[3] for r in b2rows])
    b2_err = np.hypot(np.array([r[4] for r in b2rows]), np.array([r[6] for r in b2rows]))

    def chi2_b2(d, inv_d):
        """B2(pT/A) = d(2 pT/A) / p(pT/A)^2, invariant yields per inelastic event; the deuteron
        yield uses the 1/sigma_0 fitted to the deuteron spectrum (inv_d), so B2 tests the
        pairing of the nucleons for a given deuteron normalisation"""
        s = inv_sigma0(read_summary(d)[1])
        xd, yd, _, _ = load_hist('expected_d_invyield_midy', d)
        xp, yp, _, _ = load_hist('nucleon_p_invyield_midy', d)
        dval = np.interp(2*b2_pt, xd, yd)/s*inv_d*fnd
        pval = np.interp(b2_pt, xp, yp)*fnd
        mc = dval/pval**2
        r = (mc - b2_val)**2/(b2_err**2 + (munc*b2_val)**2)
        return r.sum()/len(r), mc/b2_val

    names = [p for p in CONFIG['PARAMETERS']]
    rows = []
    for name, values in points():
        if not finished(name): continue
        d = out_dir(name); row = {'point': name, **values}; total = 0.
        for species in ('pions', 'kaons', 'protons'):
            c2, ratio = chi2_spectrum(d, species)
            row[f'chi2_{species}'] = round(c2, 2); total += obs[species]['WEIGHT']*c2
        if 'proton_tail' in obs:
            # the protons above 3 GeV on their own, so that the tail is not diluted by the many
            # bins at lower pT (fine-tune round)
            c2, _ = chi2_spectrum(d, 'protons', obs['proton_tail']['PT_RANGE'])
            row['chi2_ptail'] = round(c2, 2); total += obs['proton_tail']['WEIGHT']*c2
        c2, fits = chi2_deuterons(d)
        row['chi2_deuterons'] = round(c2, 2); total += obs['deuterons']['WEIGHT']*c2
        row['invsigma0_d'], row['invsigma0_dbar'] = round(fits[0], 2), round(fits[1], 2)
        c2, b2ratio = chi2_b2(d, fits[0])
        row['chi2_b2'] = round(c2, 2)
        if 'b2' in obs: total += obs['b2']['WEIGHT']*c2
        # rise of B2 from pT/A ~0.45 to ~1.35 GeV, MC relative to ALICE (1 = same rise)
        row['b2_rise'] = round(float(b2ratio[(b2_pt > 1.2) & (b2_pt < 1.45)].mean()/b2ratio[b2_pt < 0.6].mean()), 2)
        m = np.loadtxt(os.path.join(d, 'multiplicity_events.dat'), comments='#')
        row['dNch_deta'] = round(m[:, 2].sum()/m[:, 1].sum(), 2)
        # proton shape at a glance: MC/ALICE in a few pT ranges
        _, ratio = chi2_spectrum(d, 'protons')
        lo, hi, _, _ = alice_pikp(2)
        a, b = obs['protons']['PT_RANGE']; sel = (lo >= a - 1e-9) & (hi <= b + 1e-9); lo = lo[sel]
        for l, h in ((0.3, 0.6), (1.0, 2.0), (2.2, 2.8), (3.0, 4.0), (4.0, 5.0)):
            row[f'p_{l}-{h}'] = round(float(ratio[(lo >= l) & (lo < h)].mean()), 2)
        row['TOTAL'] = round(total, 2)
        rows.append(row)
    rows.sort(key=lambda r: r['TOTAL'])
    out = os.path.join(TUNE_SCAN, CONFIG.get('RESULTS', 'results.csv'))
    with open(out, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    short = {p: p.split(':')[-1].replace('BARYON_FRACTION', 'BF').replace('SIGMA_ND_NORM', 'SIG_ND')
               .replace('MULTI_WEIGHT_R0L0_DELTA_3/2', 'W_DELTA') for p in names}
    cols = (['point'] + names + ['chi2_pions', 'chi2_kaons', 'chi2_protons'] + (['chi2_ptail'] if 'proton_tail' in obs else [])
            + ['chi2_deuterons', 'chi2_b2', 'b2_rise', 'TOTAL', 'dNch_deta', 'p_0.3-0.6', 'p_1.0-2.0', 'p_2.2-2.8', 'p_3.0-4.0', 'p_4.0-5.0'])
    print(f'{len(rows)} finished points, best first (chi2 per point; p_a-b = protons MC/ALICE) -> {out}')
    print(' '.join(f'{short.get(c, c).replace("chi2_", "x2_"):>9s}' for c in cols))
    for r in rows:
        print(' '.join(f'{str(r[c]):>9s}' for c in cols))


if __name__ == '__main__':
    {'make': make, 'run': run, 'score': score, 'status': status}[sys.argv[1]]()
