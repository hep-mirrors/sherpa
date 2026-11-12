"""Shared helpers for the deuteron notebooks (examples/RemnantTests/SimpleLEP/XSecs) and scripts.

Every notebook does `from deuteron_common import *`.  This module only reads
files and does arithmetic; it runs no generator.  Contents:

* paths    : where the Sherpa and Pythia outputs and the HEPData tables are
* readers  : load_hist, rebin, read_summary, inv_sigma0, nevents,
             read_hepdata_blocks, alice_spectrum, f_nd
* models   : MODELS / model_path / model_hist: one interface to the three
             generator setups (Sherpa tuned, Sherpa before, Pythia 8.318)
* fitting  : fit_scale, the least-squares fit of one normalisation
* style    : apply_style, STYLE, savefig
"""
import csv, glob, os
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

from paths import DEUTERON_DATA, DEUTERON_REPORT, TUNE_SCAN

# --------------------------------------------------------------------------
# paths (doc/scripts/paths.py)
# --------------------------------------------------------------------------
DATA_DIR   = DEUTERON_DATA                      # Sherpa LEP run (Sherpa.Deuterons.yaml) and parent of all runs
HEPDATA    = os.path.join(DATA_DIR, "data")     # HEPData tables (csv) and the MCPlots file
VAL_DIR    = os.path.join(DATA_DIR, "validation_invsigma0_18")
BABAR_DIR  = os.path.join(DATA_DIR, "BaBar_10.58GeV")
OLD_ALICE  = DATA_DIR                            # ALICE_*TeV with the original cards (fixed pT_0 = 0.5)
ALICE_RUNS = os.path.join(DATA_DIR, "amisic_default")   # "Sherpa before": Amisic defaults, no reconnections
TUNE_DIR   = os.path.join(DATA_DIR, "tuneD")     # "Sherpa tuned": colour reconnections + AHADIC settings (BARYON_FRACTION 0.20)
MASTER_DIR = os.path.join(DATA_DIR, "master_CRon")   # branch hadronicRescattering: origin/master + deuteron phase, official TUNE MPI-CRon
MASTER_TUNED_DIR = os.path.join(DATA_DIR, "master_ourTune")   # same + our tuneD CR/AHADIC settings
RETUNE_DIR = os.path.join(DATA_DIR, "retuneQ7")   # master + our settings, ETA_Q 0.63, PT_Min(ref) 1.4
SCAN_DOPT_DIR = os.path.join(DATA_DIR, "scan_dopt")     # TuneScan: best single setting for the deuterons (fit_vs_pythia)
SCAN_S3_053_DIR = os.path.join(DATA_DIR, "scanTune_s3_053") # frozen ScanTune cards (deuteron_tune_report/runcards/ScanTune_s3_053), all experiments
SCAN_S3_053_KICKS_DIR = os.path.join(DATA_DIR, "scanTune_defaultKicks") # ScanTune with master's kicks PT_MAX 0.68, KT_0 1.21 (7 TeV, 3 seeds merged)
SCAN_S4_005_DIR = os.path.join(DATA_DIR, "scanTune_s4_005") # frozen cards deuteron_tune_report/runcards/ScanTune_s4_005 (fine-scan point s4_005), all experiments
SCAN_S6_008_DIR = os.path.join(DATA_DIR, "scanTune_s6_008") # frozen cards deuteron_tune_report/runcards/ScanTune_s6_008 (round 6, larger kicks); ALICE_7TeV = 4 seeds merged (400k)
SCAN_S6_010_DIR = os.path.join(DATA_DIR, "scanTune_s6_010") # frozen cards deuteron_tune_report/runcards/ScanTune_s6_010 (round 6, larger kicks); ALICE_7TeV = 4 seeds merged (400k)
LUND_TUNE2_DIR = os.path.join(DATA_DIR, "lund_tune2")       # Lund_tune1 with MPI + colour reconnections re-tuned (lund_tuning_report/runcards/Lund_tune2: SIGMA_ND_NORM 0.674, PT_0 0.955, ETA_Q 0.477, RESHUFFLE 0.416)
LUND_TUNE3_DIR = os.path.join(DATA_DIR, "lund_tune3")       # Lund_tune1 with Lund strings re-tuned to ALEPH (lund_tuning_report/runcards/Lund_tune3: aLund 0.907, bLund 0.784, aExtraDiquark 0.752)
LUND_AS_DIRS = {f: os.path.join(DATA_DIR, f"lund_as{f}") for f in ("0.5", "0.35")}   # Lund_tune1 with FS_AS_FAC = 0.5, 0.35 (LEP, BaBar only so far)
LUND_TUNE1_DIR = os.path.join(DATA_DIR, "lund_tune1")       # Lund strings with the activity re-tuned (lund_tune/runcards/Lund_tune1: SIGMA_ND_NORM 0.614, PT_0 1.21); LEP, BaBar = lund_s6_008
LUND_S6_008_DIR = os.path.join(DATA_DIR, "lund_s6_008")     # s6_008 with Pythia 8's Lund strings instead of AHADIC (lund_tuning_report/runcards/Lund_s6_008); ALICE_7TeV = 4 seeds merged (400k)
FINE_DIR   = os.path.join(TUNE_SCAN, "runs_fine")   # fine scan around s3_053 (TuneScan/scan_config_fine.yaml), 7 TeV only, 100k events per point
TUNE_V1    = os.path.join(DATA_DIR, "tuneC")     # first tuned version (BARYON_FRACTION 0.23), kept for reference
PYTHIA_DIR = os.path.join(DATA_DIR, "pythia")    # Pythia 8.318, SimpleLEP/PythiaReference/pythia_deuterons.cc
COAL_DIR   = os.path.join(DATA_DIR, "coalescence")
MULT_DIR   = os.path.join(DATA_DIR, "ALICE_7TeV_multiplicity")
NORESC_DIR = os.path.join(DATA_DIR, "norescattering")
ORIGIN_DIR = os.path.join(DATA_DIR, "origin")

CHANNELS = ["pn_dgamma", "pn_dpi0", "pp_dpiplus", "nn_dpiminus",
            "pn_dpi0pi0", "pn_dpipluspiminus", "pp_dpipluspi0", "nn_dpiminuspi0"]
LABELS = {"pn_dgamma": r"$pn\to d\gamma$", "pn_dpi0": r"$pn\to d\pi^0$",
          "pp_dpiplus": r"$pp\to d\pi^+$", "nn_dpiminus": r"$nn\to d\pi^-$",
          "pn_dpi0pi0": r"$pn\to d\pi^0\pi^0$", "pn_dpipluspiminus": r"$pn\to d\pi^+\pi^-$",
          "pp_dpipluspi0": r"$pp\to d\pi^+\pi^0$", "nn_dpiminuspi0": r"$nn\to d\pi^-\pi^0$"}
# ALICE energies: label, sqrt(s) in GeV, HEPData ins1625294 table of the (anti)deuteron spectra
ENERGIES = [("0.9", 900., 1), ("2.76", 2760., 2), ("7", 7000., 3)]

# --------------------------------------------------------------------------
# readers for the generator output
# --------------------------------------------------------------------------
def load_dat(path):
    """Text file with 2-3 numeric columns -> x, y, err (err None if absent).
    Skips comment lines and Sherpa's histogram header (more than 3 numbers)."""
    xs, ys, es = [], [], []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith(("#", "//", "%")): continue
        tok = line.split()
        if not 2 <= len(tok) <= 3: continue
        try: v = [float(t) for t in tok]
        except ValueError: continue
        xs.append(v[0]); ys.append(v[1]); es.append(v[2] if len(v) == 3 else np.nan)
    x, y, e = map(np.array, (xs, ys, es))
    o = np.argsort(x)
    return x[o], y[o], (None if np.all(np.isnan(e)) else e[o])

def load_hist(name, directory):
    """Histogram <directory>/<name>.dat (Sherpa or Pythia) -> bin centres, values,
    errors, bin width.  Values are per event and per unit of x; the overflow row
    is dropped."""
    x, y, e = load_dat(os.path.join(directory, name + ".dat"))
    x, y = x[:-1], y[:-1]
    e = np.zeros_like(y) if e is None else e[:-1]
    w = x[1] - x[0]
    return x + 0.5*w, y, e, w

def rebin(name, edges, directory, scale=1.):
    """Average a per-unit histogram over new bins `edges` (which must be whole
    multiples of the original bins); returns centres, values, errors (x scale)."""
    x, y, e, w = load_hist(name, directory)
    edges = np.asarray(edges); vals, errs = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x > lo) & (x < hi)
        vals.append((y[m]*w).sum()/(hi-lo)); errs.append(np.sqrt(((e[m]*w)**2).sum())/(hi-lo))
    return 0.5*(edges[1:]+edges[:-1]), scale*np.array(vals), scale*np.array(errs)

def read_summary(directory):
    """summary.dat -> ({key: (mean per event, error of the mean)}, header lines)."""
    values, header = {}, []
    for line in open(os.path.join(directory, "summary.dat")):
        if line.startswith("#"): header.append(line[1:].strip()); continue
        tok = line.split()
        if len(tok) == 3: values[tok[0]] = (float(tok[1]), float(tok[2]))
    return values, header

def inv_sigma0(header):
    "1/sigma_0 [1/b] the run was made with, from the summary header"
    return float(header[0].split("1/sigma_0 = ")[1].split()[0])

def nevents(header):
    "number of events of the run, from the summary header"
    return float(header[0].split("events = ")[1].split()[0])

# --------------------------------------------------------------------------
# readers for the data (HEPData csv files)
# --------------------------------------------------------------------------
def read_hepdata_blocks(path):
    "HEPData csv -> list of blocks (one per table in the file), each a list of numeric rows ('-' -> None)"
    blocks, rows = [], []
    for line in open(path):
        line = line.strip()
        if not line or line.startswith("#"):
            if rows: blocks.append(rows); rows = []
            continue
        tok = next(csv.reader([line]))
        try: rows.append([None if v.strip() == "-" else float(v) for v in tok])
        except ValueError: continue          # column header line
    if rows: blocks.append(rows)
    return blocks

def alice_spectrum(table, block):
    """ALICE (anti)deuteron spectrum, HEPData ins1625294 Table `table` (1: 0.9, 2: 2.76,
    3: 7 TeV), block 0: d, 1: dbar -> pT bin edges lo, hi, value, total error.
    Value: 1/(2 pi pT N_inel) d2N/dy dpT at |y|<0.5 [GeV^-2]; error: stat (+) syst."""
    rows = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp_deuterons_HEPData_ins1625294_Table{table}.csv"))[block]
            if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in rows]) for i in (1, 2, 3))
    err = np.hypot(0.5*(np.array([r[4] for r in rows]) - np.array([r[5] for r in rows])),
                   0.5*(np.array([r[6] for r in rows]) - np.array([r[7] for r in rows])))
    return lo, hi, val, err

def f_nd(sqrts):
    """Non-diffractive fraction of inelastic pp events, 1 - SD/INEL - DD/INEL, measured by
    ALICE (HEPData ins1181770, Tables 1, 2), and its error."""
    out, err2 = 1., 0.
    for t in (1, 2):
        for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp_diffraction_HEPData_ins1181770_Table{t}.csv"))[0]:
            if abs(r[0] - sqrts) < 1.:
                out -= r[3]; err2 += (0.5*(r[4] - r[5]))**2
    return out, np.sqrt(err2)

# --------------------------------------------------------------------------
# one interface to the generator setups
# --------------------------------------------------------------------------
# name -> (output directory of a setup, prefix of the deuteron histograms)
#   Sherpa writes "expected_*" (sum of probabilities, smooth) and "formed_*";
#   Pythia only forms deuterons by random draws, so its histograms are "formed_*".
ALL_MODELS = {
    "Sherpa tuned":  (lambda setup: os.path.join(TUNE_DIR, setup),   "expected"),
    "Sherpa before": (lambda setup: (DATA_DIR if setup == "LEP" else BABAR_DIR if setup == "BaBar_10.58GeV"
                                      else os.path.join(ALICE_RUNS, setup)), "expected"),
    "Pythia 8.318":  (lambda setup: os.path.join(PYTHIA_DIR, setup), "formed"),
    "Sherpa master": (lambda setup: os.path.join(MASTER_DIR, setup), "expected"),
    "Sherpa master + our tune": (lambda setup: os.path.join(MASTER_TUNED_DIR, setup), "expected"),
    "Sherpa retuned (new branch)": (lambda setup: os.path.join(RETUNE_DIR, setup), "expected"),
    "Sherpa d-optimised": (lambda setup: os.path.join(SCAN_DOPT_DIR, setup), "expected"),
    "Sherpa s3_053": (lambda setup: os.path.join(SCAN_S3_053_DIR, setup), "expected"),
    "Sherpa s3_053, default kicks": (lambda setup: os.path.join(SCAN_S3_053_KICKS_DIR, setup), "expected"),
    # fine-scan point, run at 7 TeV only (another setup gives a directory that does not exist)
    "Sherpa s4_005": (lambda setup: os.path.join(SCAN_S4_005_DIR, setup), "expected"),
    "Sherpa s6_008": (lambda setup: os.path.join(SCAN_S6_008_DIR, setup), "expected"),
    "Sherpa s6_010": (lambda setup: os.path.join(SCAN_S6_010_DIR, setup), "expected"),
    "Sherpa Lund":   (lambda setup: os.path.join(LUND_S6_008_DIR, setup), "expected"),
    "Sherpa Lund tuned": (lambda setup: os.path.join(LUND_TUNE1_DIR, setup), "expected"),
    "Sherpa Lund tune2": (lambda setup: os.path.join(LUND_TUNE2_DIR, setup), "expected"),
    "Sherpa Lund tune3": (lambda setup: os.path.join(LUND_TUNE3_DIR, setup), "expected"),
    "Sherpa Lund FS_AS_FAC 0.5": (lambda setup: os.path.join(LUND_AS_DIRS["0.5"], setup), "expected"),
    "Sherpa Lund FS_AS_FAC 0.35": (lambda setup: os.path.join(LUND_AS_DIRS["0.35"], setup), "expected"),
    "Sherpa s4_024": (lambda setup: os.path.join(FINE_DIR, "s4_024", "out") if setup == "ALICE_7TeV" else "", "expected"),
}
# the setups shown in the plots and tables; the others stay available in ALL_MODELS
# (e.g. "Sherpa before", "Sherpa tuned", "Sherpa d-optimised") -- add a name here to show it again
SHOWN  = ["Sherpa s3_053", "Sherpa s4_005", "Sherpa s6_008", "Sherpa s6_010", "Pythia 8.318"]
MODELS = {m: ALL_MODELS[m] for m in SHOWN}
def model_path(model, setup):
    "output directory of `model` for `setup` ('LEP', 'BaBar_10.58GeV', 'ALICE_7TeV', ...)"
    return ALL_MODELS[model][0](setup)

def available(model, setup):
    return os.path.exists(os.path.join(model_path(model, setup), "summary.dat"))

def model_scale(model, setup):
    "1/sigma_0 [1/b] the model's run used: divide its yields by this to get yields per unit of 1/sigma_0"
    return inv_sigma0(read_summary(model_path(model, setup))[1])

def deuteron_hist(model, setup, sp, var):
    """(Anti)deuteron histogram of a model per unit of 1/sigma_0: sp in 'd', 'dbar';
    var e.g. 'invyield_midy' or 'p'.  Returns centres, values, errors, width."""
    d = model_path(model, setup); prefix = ALL_MODELS[model][1]
    x, y, e, w = load_hist(f"{prefix}_{sp}_{var}", d)
    s = model_scale(model, setup)
    return x, y/s, e/s, w

def deuteron_count(model, setup, key):
    """Per-event deuteron count from summary.dat per unit of 1/sigma_0, e.g. key
    'dbar_ALEPH_window', 'dbar_total' -> (mean, error)."""
    d = model_path(model, setup); prefix = ALL_MODELS[model][1]
    v, e = read_summary(d)[0].get(f"{prefix}_{key}", (0., 0.))
    s = model_scale(model, setup)
    return v/s, e/s

# --------------------------------------------------------------------------
# fitting
# --------------------------------------------------------------------------
def fit_scale(data, err, model):
    """Least-squares normalisation x of `model` to `data` (only the data errors, as
    in Dal & Raklev, Sec. IV): minimises sum((data - x*model)/err)^2.
    Returns best x, its error, chi2 at the minimum."""
    data, err, model = map(np.asarray, (data, err, model))
    a = (model**2/err**2).sum()
    x = (data*model/err**2).sum()/a
    return x, 1./np.sqrt(a), (((data - x*model)/err)**2).sum()

def in_bins(name_or_xy, lo, hi, directory=None, scale=1.):
    """Average a histogram over each data bin [lo_i, hi_i] -> values, errors (x scale).
    Accepts a histogram name + directory, or a tuple (x, y, e, w) from load_hist."""
    x, y, e, w = load_hist(name_or_xy, directory) if isinstance(name_or_xy, str) else name_or_xy
    vals, errs = [], []
    for l, h in zip(lo, hi):
        m = (x > l) & (x < h)
        if not m.any(): m = np.abs(x - 0.5*(l+h)) <= 0.5*w     # data bin narrower than a model bin
        vals.append((y[m]*w).sum()/(w*m.sum())); errs.append(np.sqrt(((e[m]*w)**2).sum())/(w*m.sum()))
    return scale*np.array(vals), scale*np.array(errs)

def at_points(name, directory, points):
    """Histogram `name` at the points `points`, interpolated in log(y) between bin centres.
    For data bins narrower than a histogram bin (ALICE's B2 at 7 TeV: 0.05 GeV against
    0.1 GeV), where in_bins would give the same value to two neighbouring data bins."""
    x, y, e, w = load_hist(name, directory)
    ok = y > 0
    return np.exp(np.interp(points, x[ok], np.log(y[ok])))

# --------------------------------------------------------------------------
# figure style
# --------------------------------------------------------------------------
OKABE_ITO = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#F0E442", "#000000"]
STYLE = {
    "data":          dict(color="black", fmt="o", ms=4, mfc="black"),
    "Sherpa tuned":  dict(color="#D55E00", ls="-",  lw=2.0),
    "Sherpa before": dict(color="0.45",    ls="--", lw=1.4),
    "Pythia 8.318":  dict(color="#0072B2", ls="-.", lw=1.6),
    "Sherpa master": dict(color="#009E73", ls=":",  lw=2.0),
    "Sherpa master + our tune": dict(color="#CC79A7", ls="-", lw=1.8),
    "Sherpa retuned (new branch)": dict(color="#009E73", ls="-", lw=2.2),
    "Sherpa d-optimised": dict(color="#D55E00", ls="-", lw=2.4),
    "Sherpa s3_053": dict(color="#009E73", ls="--", lw=1.8),
    "Sherpa s3_053, default kicks": dict(color="#CC79A7", ls=":", lw=2.4),
    "Sherpa s4_005": dict(color="#E69F00", ls="-", lw=2.0),
    "Sherpa s4_024": dict(color="#56B4E9", ls="-", lw=2.0),
    "Sherpa s6_008": dict(color="#CC79A7", ls=":", lw=2.4),
    "Sherpa s6_010": dict(color="#D55E00", ls="-", lw=2.4),
    "Sherpa Lund":   dict(color="#E69F00", ls="-", lw=2.2),
    "Sherpa Lund tuned": dict(color="#E69F00", ls="-", lw=2.2),
    "Sherpa Lund tune2": dict(color="#009E73", ls="-", lw=2.2),
    "Sherpa Lund FS_AS_FAC 0.5": dict(color="#0072B2", ls="-", lw=1.8),
    "Sherpa Lund FS_AS_FAC 0.35": dict(color="#CC79A7", ls="-", lw=1.8),
}
MARKER = {"Sherpa tuned": "s", "Sherpa before": "D", "Pythia 8.318": "^", "Sherpa master": "o", "Sherpa master + our tune": "v", "Sherpa retuned (new branch)": "P", "Sherpa d-optimised": "*", "Sherpa s3_053": "X", "Sherpa s3_053, default kicks": "h", "Sherpa s4_005": "d", "Sherpa s4_024": "p", "Sherpa s6_008": "v", "Sherpa s6_010": "s", "Sherpa Lund": "D", "Sherpa Lund tuned": "D", "Sherpa Lund tune2": "o", "Sherpa Lund FS_AS_FAC 0.5": "s", "Sherpa Lund FS_AS_FAC 0.35": "^"}   # for point plots

def apply_style():
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["Times New Roman", "Times", "STIXGeneral"],
        "mathtext.fontset": "stix",
        "font.size": 12, "axes.titlesize": 12, "axes.labelsize": 13, "legend.fontsize": 9,
        "xtick.labelsize": 11, "ytick.labelsize": 11,
        "axes.prop_cycle": cycler(color=OKABE_ITO),
        "axes.grid": True, "grid.alpha": 0.25, "grid.linestyle": ":",
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.minor.visible": True, "ytick.minor.visible": True,
        "lines.linewidth": 1.6, "errorbar.capsize": 0,
        "legend.frameon": False, "figure.dpi": 110, "savefig.dpi": 300, "savefig.bbox": "tight",
    })

SAVE_FIGS = False
def savefig(fig, name):
    "write deuteron_report/figures/<name>.pdf if SAVE_FIGS is set (deuteron_common.SAVE_FIGS = True)"
    if SAVE_FIGS: fig.savefig(os.path.join(DEUTERON_REPORT, "figures", name + ".pdf"))

apply_style()
