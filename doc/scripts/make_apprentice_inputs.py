"""Turns finished TuneScan points into apprentice input (https://github.com/HEPonHPC/apprentice,
here the fork in ~/code_software/apprentice).

    python3 make_apprentice_inputs.py <runs dir> <output dir>
        e.g. ../examples/RemnantTests/SimpleLEP/TuneScan/runs_round3 ../examples/RemnantTests/SimpleLEP/TuneScan/apprentice/round3

writes
    <output dir>/runs/<point>/params.dat    the point's parameters (from the header of card.yaml)
    <output dir>/runs/<point>/sherpa.yoda   the observables below, on ALICE's bins
    <output dir>/data.json                  ALICE's values and errors, {"/obs#bin": [value, error]}
    <output dir>/weights_pikp.txt           tune on pions, kaons, protons (weights as scan_config)
    <output dir>/weights_deuterons.txt      the deuterons, only to build their surrogates

Observables, all pp 7 TeV at |y|<0.5 and per inelastic event (non-diffractive yields x f_ND),
with the same binning and normalisation as scorecard.py:
    /ALICE_7TeV/pions, /kaons, /protons   1/N_inel d2N/dy dpT [1/GeV], (pi+ + pi-), (K+ + K-), (p + pbar),
                                          ALICE HEPData ins1357424 Table 1, 0.3 < pT < 5 GeV
    /ALICE_7TeV/d, /dbar                  1/(2 pi pT N_inel) d2N/dy dpT [1/GeV^2] PER UNIT OF 1/sigma_0
                                          (1/b), ALICE HEPData ins1625294 Table 3; multiply by a
                                          fitted 1/sigma_0 to compare with data
The data errors in data.json are stat (+) syst (+) 5 % of the value, as MODEL_UNCERTAINTY in
scan_config.yaml.  The weight of a bin is sqrt(W/n) (W = observable weight, n = its bins), so that
apprentice's chi2 = sum over observables of W x chi2/n, the TOTAL of tune_scan.py."""
import json, os, re, sys
import numpy as np
import yoda
from paths import XSECS
ARGS = [os.path.abspath(a) for a in sys.argv[1:]]           # before the chdir below
os.chdir(XSECS)                                             # the notebooks' folder, as when they run
from deuteron_common import alice_spectrum, f_nd, in_bins, inv_sigma0, load_hist, read_summary
import contextlib, io
with contextlib.redirect_stdout(io.StringIO()):               # scorecard.py prints its Pythia card on import
    from scorecard import alice_pikp

MODEL_UNCERTAINTY = 0.05
SQRTS = 7000.
# observable -> (histograms summed, ALICE block of ins1357424 Table 1, weight in the tune)
SPECTRA = {"pions":   (["meson_pi_plus", "meson_pi_minus"], 0, 1.0),
           "kaons":   (["meson_K_plus", "meson_K_minus"],   1, 1.0),
           "protons": (["nucleon_p", "nucleon_pbar"],       2, 2.0)}
DEUTERONS = {"d": 0, "dbar": 1}            # block of ins1625294 Table 3 (7 TeV)
# ALICE V0M multiplicity classes at 7 TeV (arXiv:1902.09290, Table 1): percentile range (busiest
# first) and <dNch/deta> at |eta|<0.5 with its error; weight of the observable in the tune
NCH_CLASSES = [(0., 4.7, 17.47, 0.52), (4.7, 9.5, 13.50, 0.40), (9.5, 19., 10.76, 0.30),
               (19., 38., 7.54, 0.23), (38., 100., 3.30, 0.13)]
NCH_WEIGHT = 2.0
# card path in tune_scan.py -> short name in params.dat
SHORT = {"AMISIC:SIGMA_ND_NORM": "SIGMA_ND_NORM", "AMISIC:PT_0(ref)": "PT_0",
         "COLOUR_RECONNECTIONS:ETA_Q": "ETA_Q", "AHADIC:PT_MAX": "PT_MAX", "AHADIC:KT_0": "KT_0",
         "AHADIC:BARYON_FRACTION": "BARYON_FRACTION", "AHADIC:MULTI_WEIGHT_R0L0_DELTA_3/2": "DELTA_WEIGHT",
         "AHADIC:ALPHA_L": "ALPHA_L", "AHADIC:BETA_L": "BETA_L", "AHADIC:GAMMA_L": "GAMMA_L",
         "AHADIC:TRANSITION_THRESHOLD": "TRANSITION_THRESHOLD", "AHADIC:DECAY_THRESHOLD": "DECAY_THRESHOLD"}

def card_parameters(card):
    "the '#   <card path> = <value>' lines tune_scan.py writes at the top of card.yaml"
    params = {}
    for line in open(card):
        if not line.startswith("#"): break
        match = re.match(r"#\s+(\S+) = (\S+)", line)
        if match: params[SHORT.get(match.group(1), match.group(1))] = float(match.group(2))
    return params

def estimate(path, lo, hi, values, errors):
    "a YODA Estimate1D on the bins [lo_i, hi_i] (ALICE's bins are contiguous)"
    edges = list(lo) + [hi[-1]]
    assert np.allclose(lo[1:], hi[:-1]), path
    est = yoda.Estimate1D(edges, path)
    for i, (v, e) in enumerate(zip(values, errors)):
        est.bin(i+1).setVal(float(v)); est.bin(i+1).setErr((-float(e), float(e)))
    return est

def nch_classes(outdir):
    """<dNch/deta> (|eta|<0.5) in ALICE's V0M classes: INEL>0 events ranked by their V0M count
    and cut at ALICE's percentiles (fractional rows at the edges), as ALICE_2019_Multiplicity"""
    m = np.atleast_2d(np.loadtxt(os.path.join(outdir, "multiplicity_events.dat")))
    order = np.argsort(-m[:, 0]); nev, nch = m[order, 1], m[order, 2]
    hi = np.cumsum(nev)/nev.sum()*100; lo = hi - nev/nev.sum()*100
    values = []
    for plo, phi, _, _ in NCH_CLASSES:
        w = np.clip(np.minimum(hi, phi) - np.maximum(lo, plo), 0, None)/np.maximum(hi - lo, 1e-12)
        values.append((w*nch).sum()/(w*nev).sum())
    return np.array(values)

def observables(outdir):
    "the YODA objects of one finished point"
    fnd, _ = f_nd(SQRTS)
    objects = []
    for name, (hists, block, _) in SPECTRA.items():
        lo, hi, _, _ = alice_pikp(block)
        total = None
        for h in hists:
            x, y, e, w = load_hist(h + "_invyield_midy", outdir)
            v, ve = 2*np.pi*x*y*fnd, 2*np.pi*x*e*fnd           # invariant yield -> d2N/dy dpT
            total = (x, v, ve, w) if total is None else (x, total[1] + v, np.hypot(total[2], ve), w)
        mc, mce = in_bins(total, lo, hi)
        objects.append(estimate(f"/ALICE_7TeV/{name}", lo, hi, mc, mce))
    s = inv_sigma0(read_summary(outdir)[1])
    for sp, block in DEUTERONS.items():
        lo, hi, _, _ = alice_spectrum(3, block)
        x, y, e, w = load_hist(f"expected_{sp}_invyield_midy", outdir)
        mc, mce = in_bins((x, y/s, e/s, w), lo, hi, scale=fnd)
        objects.append(estimate(f"/ALICE_7TeV/{sp}", lo, hi, mc, mce))
    if os.path.exists(os.path.join(outdir, "multiplicity_events.dat")):
        v = nch_classes(outdir)                     # bins 0..4 = classes I+II ... VIII+IX+X
        objects.append(estimate("/ALICE_7TeV/nch_classes", np.arange(5.), np.arange(1., 6.), v, 0.01*v))
    return objects

def main(runs_dir, out):
    os.makedirs(os.path.join(out, "runs"), exist_ok=True)
    used = 0
    for point in sorted(os.listdir(runs_dir)):
        card, outdir = os.path.join(runs_dir, point, "card.yaml"), os.path.join(runs_dir, point, "out")
        if not os.path.exists(os.path.join(outdir, "summary.dat")): continue
        params = card_parameters(card)
        if not params: continue
        target = os.path.join(out, "runs", point)
        os.makedirs(target, exist_ok=True)
        with open(os.path.join(target, "params.dat"), "w") as f:
            for name, value in params.items(): f.write(f"{name} {value}\n")
        yoda.write(observables(outdir), os.path.join(target, "sherpa.yoda"))
        used += 1
    data = {}
    weights = {"pikp": [], "deuterons": []}
    for name, (_, block, weight) in SPECTRA.items():
        lo, hi, val, err = alice_pikp(block)
        for i, (v, e) in enumerate(zip(val, err)):
            data[f"/ALICE_7TeV/{name}#{i}"] = [float(v), float(np.hypot(e, MODEL_UNCERTAINTY*v))]
        weights["pikp"].append(f"/ALICE_7TeV/{name} {np.sqrt(weight/len(val)):.6f}")
    for i, (_, _, v, e) in enumerate(NCH_CLASSES):
        data[f"/ALICE_7TeV/nch_classes#{i}"] = [v, float(np.hypot(e, MODEL_UNCERTAINTY*v))]
    weights["pikp"].append(f"/ALICE_7TeV/nch_classes {np.sqrt(NCH_WEIGHT/len(NCH_CLASSES)):.6f}")
    for sp, block in DEUTERONS.items():
        lo, hi, val, err = alice_spectrum(3, block)
        for i, (v, e) in enumerate(zip(val, err)):
            data[f"/ALICE_7TeV/{sp}#{i}"] = [float(v), float(e)]
        weights["deuterons"].append(f"/ALICE_7TeV/{sp} 1")
    with open(os.path.join(out, "data.json"), "w") as f: json.dump(data, f, indent=1)
    for key, lines in weights.items():
        with open(os.path.join(out, f"weights_{key}.txt"), "w") as f: f.write("\n".join(lines) + "\n")
    print(f"{used} points -> {out}/runs, data.json, weights_pikp.txt, weights_deuterons.txt")

if __name__ == "__main__":
    if len(ARGS) != 2: sys.exit(__doc__)
    main(*ARGS)
