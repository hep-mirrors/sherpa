"""All figures of tunedLund.tex, one PDF each, into doc/lund_tune/figures/.

    python3 make_figures.py            all figures of tunedLund.tex -> doc/lund_tune/figures/
    python3 make_figures.py final      the final setting (Lund_mbtune), Lund_tune1 and Pythia
                                       -> doc/final_report/figures/deuterons/ (combined report)

Generators: Sherpa with Lund strings, Lund_tune1 and Lund_tune4, and Pythia 8.318.
Deuteron models: the Dal & Raklev cross-section model ("CrossSection") and the coalescence models
HardDisk and Exponential, with the momentum condition only (_k) or also b < 3 fm (_kb).
Every figure compares the generators with the data of one experiment; the free parameter
(1/sigma_0, or the coalescence momentum p_0) is fitted per data set, as in Dal & Raklev.

Inputs (doc/examples/RemnantTests/SimpleLEP):
  DeuteronFormationFromXSecs/lund_tune1, lund_tune4, pythia       cross-section model
  DeuteronFormationFromXSecs/coalescence_lund/<tune1|tune4|pythia>/<model>/<setup>
                                                                  coalescence (run_coalescence.py)
  DeuteronFormationFromXSecs/data                                 HEPData tables
  LEPRivet/lund_tune1, LEPRivet/runs_lund_lep1/app_tune           charged particles at the Z pole
The fit, chi2 and data helpers are those of XSecs/deuteron_common.py and TuneScan/scorecard.py, so
the numbers agree with the earlier notebooks and reports."""
import contextlib, glob, io, json, os, sys
import numpy as np

from paths import LUND_TUNE, XSECS, LEP_RIVET
FIG = os.path.join(LUND_TUNE, "figures")
os.chdir(XSECS)                                          # the notebooks' folder, as when they run
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from deuteron_common import (DATA_DIR, HEPDATA, LABELS, alice_spectrum, f_nd, fit_scale, in_bins, inv_sigma0,
                             load_hist, read_hepdata_blocks, read_summary)
with contextlib.redirect_stdout(io.StringIO()):
    import scorecard                                      # prints its Pythia scorecard on import

DATA = os.path.abspath(DATA_DIR)
COAL = os.path.join(DATA, "coalescence_lund")
GENERATORS = ["tune1", "tune4", "pythia"]
GEN_LABEL = {"tune1": "Sherpa Lund_tune1", "tune4": "Sherpa Lund_tune4", "pythia": "Pythia 8.318"}
GEN_STYLE = {"tune1": dict(color="#E69F00", ls="-", lw=2.0), "tune4": dict(color="#009E73", ls="-", lw=2.0),
             "pythia": dict(color="#0072B2", ls="-.", lw=1.8)}
GEN_MARKER = {"tune1": "D", "tune4": "o", "pythia": "^"}
# the final setting: Lund_tune1 with the minimum-bias tune (make_mbtune_cards.py validate6 p018)
GEN_LABEL["mbtune"] = "Sherpa final setting"
GEN_STYLE["mbtune"] = dict(color="#CC0000", ls="-", lw=2.2)
GEN_MARKER["mbtune"] = "s"
FINAL = len(sys.argv) > 1 and sys.argv[1] == "final"
TABLE = os.path.join(LUND_TUNE, "summary_table.tex")
if FINAL:
    # one Sherpa setting only, the one used at LEP and in pp
    GENERATORS = ["mbtune", "pythia"]
    FIG = os.path.join(os.path.dirname(LUND_TUNE), "final_report", "figures", "deuterons")
    TABLE = os.path.join(os.path.dirname(LUND_TUNE), "final_report", "summary_table.tex")
MODELS = ["CrossSection", "HardDisk_k", "Exponential_k", "HardDisk_kb", "Exponential_kb"]
MODEL_LABEL = {"CrossSection": "cross-section model", "HardDisk_k": r"HardDisk, $k<p_0$",
               "Exponential_k": r"Exponential, $k$ only", "HardDisk_kb": r"HardDisk, $k<p_0$, $b<3$ fm",
               "Exponential_kb": r"Exponential, $k$ and $b$ ($b_{\rm coal}=3$ fm)"}
P_RUN = 0.2                                               # p_coal of the coalescence runs [GeV]
ENERGIES = [("0.9", 900., 1), ("2.76", 2760., 2), ("7", 7000., 3)]
SPECIES_LABEL = {"d": "d", "dbar": r"\bar d"}

# ------------------------------------------------------------------ runs
def run_dir(gen, model, setup):
    if model == "CrossSection":
        return os.path.join(DATA, {"tune1": "lund_tune1", "tune4": "lund_tune4", "mbtune": "lund_mbtune", "pythia": "pythia"}[gen], setup)
    merged = os.path.join(COAL, gen, model, setup + "_merged")     # first run + extra seeds (run_coalescence.py merge)
    if os.path.exists(os.path.join(merged, "summary.dat")): return merged
    return os.path.join(COAL, gen, model, setup)

def exists(d): return os.path.exists(os.path.join(d, "summary.dat"))

def prefix(gen, model): return "formed" if (gen == "pythia" and model == "CrossSection") else "expected"

def run_scale(model, d):
    "the run's own value of the model parameter: yields / run_scale = yields per unit of 1/sigma_0 (or per run)"
    return inv_sigma0(read_summary(d)[1]) if model == "CrossSection" else 1.

def parameter(model, s, se):
    "fitted scale -> (value, error, label): 1/sigma_0 [1/b], or p_0 = p_run s^(1/3) [GeV] (yield ~ p_0^3)"
    if model == "CrossSection": return s, se, rf"$1/\sigma_0$ = {s:.2f} b$^{{-1}}$"
    if s <= 0: return 0., 0., r"$p_0$ = 0"
    p = P_RUN*s**(1./3.)
    return p, p*se/(3.*s), rf"$p_0$ = {1e3*p:.0f} MeV"

def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf")); plt.close(fig); print("  figures/" + name + ".pdf")

# ------------------------------------------------------------------ ALICE 2018: spectra at 0.9, 2.76, 7 TeV
def alice_fit(gen, model, rs, sqrts, table, sp):
    d = run_dir(gen, model, f"ALICE_{rs}TeV")
    if not exists(d): return None
    fnd, _ = f_nd(sqrts); s0 = run_scale(model, d)
    lo, hi, val, err = alice_spectrum(table, 0 if sp == "d" else 1)
    x, y, e, w = load_hist(f"{prefix(gen, model)}_{sp}_invyield_midy", d)
    mc, mce = in_bins((x, y/s0, e/s0, w), lo, hi, scale=fnd)
    s, se, chi2_data = fit_scale(val, err, mc)
    # chi2 with the data errors and the generator's statistical error (small for the cross-section
    # model, not for coalescence, which uses only the few pairs at small k).  A bin without any MC
    # entry gets the error of one entry (the smallest non-zero error of the spectrum).
    floor = mce[mce > 0].min() if (mce > 0).any() else 0.
    mce = np.where(mce > 0, mce, floor)
    chi2 = (((val - s*mc)**2)/(err**2 + (s*mce)**2)).sum()
    return dict(lo=lo, hi=hi, val=val, err=err, mc=s*mc, mce=s*mce, s=s, se=se, chi2=chi2, chi2_data=chi2_data, n=len(val))

def fig_alice_spectra(model):
    for sp in ("d", "dbar"):
        fig, axes = plt.subplots(2, 3, figsize=(16, 6.4), sharex="col", gridspec_kw={"height_ratios": [3, 1.3]})
        any_run = False
        for col, (rs, sqrts, table) in enumerate(ENERGIES):
            ax, rx = axes[0, col], axes[1, col]
            data = None
            for gen in GENERATORS:
                f = alice_fit(gen, model, rs, sqrts, table, sp)
                if f is None: continue
                any_run = True; data = f
                c = 0.5*(f["lo"] + f["hi"])
                _, _, lab = parameter(model, f["s"], f["se"])
                ax.plot(c, f["mc"], marker=GEN_MARKER[gen], ms=3.5,
                        label=rf"{GEN_LABEL[gen]}: {lab}, $\chi^2$/n = {f['chi2']:.1f}/{f['n']}", **GEN_STYLE[gen])
                rx.plot(c, f["mc"]/f["val"], marker=GEN_MARKER[gen], ms=3.5, **GEN_STYLE[gen])
            if data is not None:
                c = 0.5*(data["lo"] + data["hi"])
                ax.errorbar(c, data["val"], yerr=data["err"], xerr=0.5*(data["hi"] - data["lo"]), fmt="o", color="black", ms=4, label="ALICE", zorder=5)
                rx.fill_between(c, 1 - data["err"]/data["val"], 1 + data["err"]/data["val"], color="0.85", step="mid", label="ALICE error")
            ax.set_yscale("log"); ax.legend(fontsize=7.5, loc="lower left")
            ax.set_title(rf"${SPECIES_LABEL[sp]}$, pp $\sqrt{{s}}$ = {rs} TeV, $|y|<0.5$")
            rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0, 2); rx.set_xlabel(r"$p_T$ [GeV]")
        axes[0, 0].set_ylabel(r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]"); axes[1, 0].set_ylabel("MC / ALICE")
        fig.suptitle(f"ALICE (anti)deuterons, {MODEL_LABEL[model]}", y=1.0)
        plt.tight_layout()
        if any_run: save(fig, f"alice_spectra_{sp}_{model}")
        else: plt.close(fig)

# ------------------------------------------------------------------ LEP (ALEPH, OPAL) and BaBar
A_ALEPH, SA_ALEPH = 5.9e-6, np.hypot(1.8e-6, 0.5e-6)     # dbar per hadronic Z in ALEPH's window
EPS_OPAL, N_OPAL = 0.234, 1.64e6                         # OPAL: efficiency, hadronic Z decays (no candidate)
def lep_fit(gen, model):
    d = run_dir(gen, model, "LEP")
    if not exists(d): return None
    v, _ = read_summary(d); s0 = run_scale(model, d); pre = prefix(gen, model)
    a = v.get(f"{pre}_dbar_ALEPH_window", (0., 0.))[0]/s0; b = v.get(f"{pre}_dbar_OPAL_window", (0., 0.))[0]/s0
    if a <= 0: return None
    x = max(0., (A_ALEPH - EPS_OPAL*N_OPAL*b*SA_ALEPH**2/(2*a))/a)
    return dict(s=x, se=SA_ALEPH/a, chi2=(A_ALEPH - a*x)**2/SA_ALEPH**2 + EPS_OPAL*N_OPAL*b*x, a=a, d=d, s0=s0, pre=pre)

babar_row = read_hepdata_blocks(os.path.join(HEPDATA, "BaBar_continuum_ratio_HEPData_ins1286317_Table5.csv"))[0][0]
R_BABAR = babar_row[1]; RE_BABAR = np.hypot(0.5*(babar_row[2] - babar_row[3]), 0.5*(babar_row[4] - babar_row[5]))
def babar_fit(gen, model):
    d = run_dir(gen, model, "BaBar_10.58GeV")
    if not exists(d): return None
    v, _ = read_summary(d); a = v.get(f"{prefix(gen, model)}_dbar_total", (0., 0.))[0]/run_scale(model, d)
    if a <= 0: return None
    return dict(s=R_BABAR/a, se=RE_BABAR/a)

def fig_lep_spectrum(model):
    "the antideuteron momentum spectrum at the Z pole at each generator's fitted parameter"
    edges = np.arange(0., 3.001, 0.1); c = 0.5*(edges[1:] + edges[:-1])
    fig, ax = plt.subplots(figsize=(8, 5)); any_run = False
    for gen in GENERATORS:
        f = lep_fit(gen, model)
        if f is None: continue
        any_run = True
        x, y, e, w = load_hist(f"{f['pre']}_dbar_p", f["d"])
        v = np.array([(y[(x > l) & (x < h)]*w).sum()/(h - l) for l, h in zip(edges[:-1], edges[1:])])
        _, _, lab = parameter(model, f["s"], f["se"])
        ax.plot(c, f["s"]*v/f["s0"], label=f"{GEN_LABEL[gen]}: {lab}", **GEN_STYLE[gen])
    ax.axvspan(0.62, 1.03, color="grey", alpha=0.15, label="ALEPH window")
    ax.errorbar(0.825, A_ALEPH/0.41/0.95, yerr=SA_ALEPH/0.41/0.95, xerr=0.205, fmt="o", color="black", ms=6, zorder=5,
                label=r"ALEPH, $(5.9\pm1.9)\times10^{-6}$ in the window")
    ax.set_xlabel(r"$p$ [GeV]"); ax.set_ylabel(r"$\bar d$ per hadronic Z decay per GeV"); ax.set_xlim(0, 3); ax.set_ylim(0, None)
    ax.set_title(rf"$e^+e^-\to Z\to$ hadrons, $\sqrt{{s}}$ = 91.2 GeV, {MODEL_LABEL[model]}"); ax.legend(fontsize=8)
    plt.tight_layout()
    if any_run: save(fig, f"lep_spectrum_{model}")
    else: plt.close(fig)

def fig_consistency(model):
    "the fitted parameter per data set and generator: one value for all would be a perfect model"
    sets = ["LEP (ALEPH+OPAL)", "BaBar"] + [f"ALICE {rs} TeV ${SPECIES_LABEL[sp]}$" for rs, _, _ in ENERGIES for sp in ("d", "dbar")]
    fig, ax = plt.subplots(figsize=(11, 4.6)); any_run = False; table = {}
    for ig, gen in enumerate(GENERATORS):
        vals, errs = [], []
        for f in [lep_fit(gen, model), babar_fit(gen, model)] + [alice_fit(gen, model, rs, sq, t, sp) for rs, sq, t in ENERGIES for sp in ("d", "dbar")]:
            p, pe, _ = parameter(model, f["s"], f["se"]) if f else (np.nan, np.nan, "")
            vals.append(p); errs.append(pe)
        if np.all(np.isnan(vals)): continue
        any_run = True; table[gen] = vals
        xs = np.arange(len(sets)) + (ig - 1)*0.18
        ax.errorbar(xs, vals, yerr=errs, fmt=GEN_MARKER[gen], color=GEN_STYLE[gen]["color"], ms=6, label=GEN_LABEL[gen])
    ax.set_xticks(np.arange(len(sets))); ax.set_xticklabels(sets, rotation=20, ha="right", fontsize=9)
    ax.set_ylabel(r"fitted $1/\sigma_0$ [b$^{-1}$]" if model == "CrossSection" else r"fitted $p_0$ [GeV]")
    ax.set_ylim(0, None); ax.legend(fontsize=8); ax.set_title(f"The model parameter fitted to each measurement, {MODEL_LABEL[model]}")
    plt.tight_layout()
    if any_run: save(fig, f"consistency_{model}")
    else: plt.close(fig)
    return sets, table

# ------------------------------------------------------------------ B2 at 7 TeV
def fig_b2(model):
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7.5, 6.2), sharex=True, gridspec_kw={"height_ratios": [3, 1.3]})
    fnd, _ = f_nd(7000.); any_run = False
    pt, val, err = scorecard.b2_pt, scorecard.b2_val, scorecard.b2_err
    for gen in GENERATORS:
        d = run_dir(gen, model, "ALICE_7TeV"); f = alice_fit(gen, model, "7", 7000., 3, "d")
        if f is None: continue
        any_run = True
        s0 = run_scale(model, d)
        xd, yd, _, _ = load_hist(f"{prefix(gen, model)}_d_invyield_midy", d); xp, yp, _, _ = load_hist("nucleon_p_invyield_midy", d)
        mc = (np.interp(2*pt, xd, yd)/s0*f["s"]*fnd)/(np.interp(pt, xp, yp)*fnd)**2
        chi2 = (((mc - val)/err)**2).sum()
        ax.plot(pt, mc*1e3, marker=GEN_MARKER[gen], ms=3.5, label=rf"{GEN_LABEL[gen]}, $\chi^2$/n = {chi2:.1f}/{len(val)}", **GEN_STYLE[gen])
        rx.plot(pt, mc/val, marker=GEN_MARKER[gen], ms=3.5, **GEN_STYLE[gen])
    ax.errorbar(pt, val*1e3, yerr=err*1e3, fmt="o", color="black", ms=4, label="ALICE", zorder=5)
    rx.fill_between(pt, 1 - err/val, 1 + err/val, color="0.85", step="mid")
    ax.set_ylabel(r"$B_2$ [$10^{-3}$ GeV$^2$]"); ax.set_ylim(0, None); ax.legend(fontsize=8)
    ax.set_title(rf"$B_2$ of deuterons, pp 7 TeV, {MODEL_LABEL[model]}")
    rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0, 2); rx.set_xlabel(r"$p_T/A$ [GeV]"); rx.set_ylabel("MC / ALICE")
    plt.tight_layout()
    if any_run: save(fig, f"b2_7TeV_{model}")
    else: plt.close(fig)

# ------------------------------------------------------------------ ALICE 2019: multiplicity classes at 7 TeV
CLASSES = [("I+II", 0., 4.7, 17.47, 0.52), ("III", 4.7, 9.5, 13.50, 0.40), ("IV+V", 9.5, 19., 10.76, 0.30),
           ("VI+VII", 19., 38., 7.54, 0.23), ("VIII+IX+X", 38., 100., 3.30, 0.13)]   # arXiv:1902.09290, Table 1
MULT_SPECIES = ["d", "dbar", "p", "pbar"]
PT_C = np.arange(0.05, 5.0, 0.1)

def load_rows(path): return {int(r[0]): np.array(r[1:]) for r in np.atleast_2d(np.loadtxt(path))}

def mult_dirs(gen, model):
    if model == "CrossSection":
        base = os.path.join(run_dir(gen, model, "ALICE_7TeV_multiplicity"))
        return sorted(d for d in glob.glob(os.path.join(base, "seed_*")) if os.path.isdir(d))
    first = os.path.join(COAL, gen, model, "ALICE_7TeV")
    dirs = [first] + sorted(glob.glob(os.path.join(COAL, gen, model, "ALICE_7TeV_seeds", "seed_*")))
    return [d for d in dirs if os.path.exists(os.path.join(d, "multiplicity_events.dat"))]

def build_classes(dirs):
    "INEL>0 events of the runs in `dirs`, cut into ALICE's V0M classes (as ALICE_2019_Multiplicity.ipynb)"
    ev, yl = {}, {sp: {} for sp in MULT_SPECIES}
    for sd in dirs:
        for k, v in load_rows(os.path.join(sd, "multiplicity_events.dat")).items(): ev[k] = ev.get(k, 0) + v
        for sp in MULT_SPECIES:
            for k, v in load_rows(os.path.join(sd, f"multiplicity_{sp}_invyield_midy.dat")).items(): yl[sp][k] = yl[sp].get(k, 0) + v
    nv0 = np.array(sorted(ev)); nev = np.array([ev[k][0] for k in nv0])
    order = np.argsort(-nv0); cum_hi = np.cumsum(nev[order])/nev.sum()*100; cum_lo = cum_hi - nev[order]/nev.sum()*100
    out = []
    for name, plo, phi, _, _ in CLASSES:
        w = np.zeros(len(nv0))
        w[order] = np.clip(np.minimum(cum_hi, phi) - np.maximum(cum_lo, plo), 0, None)/np.maximum(cum_hi - cum_lo, 1e-12)
        N = (w*nev).sum()
        c = {"nch": (w*np.array([ev[k][1] for k in nv0])).sum()/N}
        for j, sp in enumerate(MULT_SPECIES):
            c[sp + "_dndy"] = (w*np.array([ev[k][2 + j] for k in nv0])).sum()/N
            c[sp] = (w[:, None]*np.array([yl[sp][k] for k in nv0])).sum(axis=0)/N
        out.append(c)
    return out

def fig_multiplicity(model):
    data = []
    for ic in range(len(CLASSES)):
        b = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table{ic+1}.csv"))[0] if r[3] is not None]
        lo, hi, val = (np.array([r[i] for r in b]) for i in (1, 2, 3))
        err = np.hypot(0.5*(np.array([r[4] for r in b]) - np.array([r[5] for r in b])), 0.5*(np.array([r[6] for r in b]) - np.array([r[7] for r in b])))
        data.append((lo, hi, val, err))
    t12 = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, "ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table12.csv"))[0] if r[3] is not None]
    dp_x = np.array([r[0] for r in t12]); dp_v = np.array([r[3] for r in t12]); dp_e = np.array([np.hypot(r[4], r[6]) for r in t12])
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(17, 4.6)); any_run = False
    a1.errorbar(dp_x, dp_v, yerr=dp_e, fmt="o", color="black", ms=5, label="ALICE", zorder=5)
    lowest = int(np.argmin(dp_x))
    a2.plot(dp_x, dp_v/dp_v[lowest], "o", color="black", ms=5, label="ALICE")
    b2data = []
    for ic in range(len(CLASSES)):
        rows = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table{6+ic}.csv"))[0] if r[3] is not None]
        sel = [r for r in rows if r[1] >= 0.7 - 1e-9 and r[2] <= 0.8 + 1e-9]
        if sel: b2data.append((CLASSES[ic][3], np.mean([r[3] for r in sel]), np.mean([np.hypot(r[4], r[6]) for r in sel])))
    a3.errorbar([b[0] for b in b2data], [1e3*b[1] for b in b2data], yerr=[1e3*b[2] for b in b2data], fmt="o", color="black", ms=5, label="ALICE", zorder=5)
    for gen in GENERATORS:
        dirs = mult_dirs(gen, model)
        if not dirs: continue
        any_run = True
        cls = build_classes(dirs)
        dndpt = lambda c, l, h: (2*np.pi*PT_C[(PT_C > l) & (PT_C < h)]*c["d"][(PT_C > l) & (PT_C < h)]*0.1).sum()/(h - l)
        mcs = [np.array([dndpt(c, l, h) for l, h in zip(lo, hi)]) for c, (lo, hi, _, _) in zip(cls, data)]
        s, _, _ = fit_scale(np.concatenate([d[2] for d in data]), np.concatenate([d[3] for d in data]), np.concatenate(mcs))
        nch = np.array([c["nch"] for c in cls])
        dp = np.array([2*s*c["d_dndy"]/(c["p_dndy"] + c["pbar_dndy"]) for c in cls])
        chi2 = (((np.interp(dp_x, nch[::-1], dp[::-1]) - dp_v)/dp_e)**2).sum()
        a1.plot(nch, dp, marker=GEN_MARKER[gen], label=rf"{GEN_LABEL[gen]}, $\chi^2$/n = {chi2:.0f}/{len(dp_v)}", **GEN_STYLE[gen])
        a2.plot(nch, dp/dp[-1], marker=GEN_MARKER[gen], label=f"{GEN_LABEL[gen]}: rise x{dp[0]/dp[-1]:.1f}", **GEN_STYLE[gen])
        md, mp = (PT_C > 1.4) & (PT_C < 1.6), (PT_C > 0.7) & (PT_C < 0.8)
        b2 = [1e3*s*0.5*(c["d"][md] + c["dbar"][md]).mean()/(0.5*(c["p"][mp] + c["pbar"][mp]).mean())**2 for c in cls]
        a3.plot(nch, b2, marker=GEN_MARKER[gen], label=GEN_LABEL[gen], **GEN_STYLE[gen])
    a2.lines[0].set_label(f"ALICE: rise x{dp_v[int(np.argmax(dp_x))]/dp_v[lowest]:.1f}")
    for a, yl, t in [(a1, r"$2d/(p+\bar p)$", r"$d/p$ against multiplicity"), (a2, "relative to the lowest class", r"trend only (no $\sigma_0$ / $p_0$)"),
                     (a3, r"$B_2$ at $p_T/A$ = 0.7–0.8 GeV [$10^{-3}$ GeV$^2$]", r"$B_2$ against multiplicity")]:
        a.set_xlabel(r"$\langle dN_{\rm ch}/d\eta\rangle_{|\eta|<0.5}$"); a.set_ylabel(yl); a.set_title(t); a.legend(fontsize=8); a.set_ylim(0, None)
    fig.suptitle(f"ALICE multiplicity classes, pp 7 TeV, {MODEL_LABEL[model]}", y=1.02)
    plt.tight_layout()
    if any_run: save(fig, f"multiplicity_{model}")
    else: plt.close(fig)

# ------------------------------------------------------------------ inputs: sigma(k), pairs, charged particles
def fig_sigma_k():
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8))
    groups = [["pn_dgamma", "pn_dpi0", "pn_dpi0pi0", "pn_dpipluspiminus"], ["pp_dpiplus", "pp_dpipluspi0", "nn_dpiminus", "nn_dpiminuspi0"]]
    for ax, group, title in zip(axes, groups, [r"$pn\to dX$ (Dal & Raklev, Fig. 4)", r"$pp\to dX$, $nn\to dX$ (Fig. 5)"]):
        for ch in group:
            x, y, _, _ = load_hist("test_sigma_" + ch, DATA)
            ax.plot(x, y, label=LABELS[ch])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 4.); ax.set_ylim(1., 5e3)
        ax.set_xlabel(r"$k$ [GeV]"); ax.set_ylabel(r"$\sigma$ [$\mu$b]"); ax.set_title(title); ax.legend(fontsize=8)
    plt.tight_layout(); save(fig, "sigma_k")

def fig_pairs_k():
    "pn pairs per event per GeV of k, and which k form the deuterons in each model"
    setups = ["LEP", "ALICE_0.9TeV", "ALICE_2.76TeV", "ALICE_7TeV"]
    fig, axes = plt.subplots(1, 5, figsize=(24, 4.4))
    low_k = {gen: [] for gen in GENERATORS}               # pn pairs per event with k < 0.2 GeV
    for ax, setup in zip(axes, setups):
        for gen in GENERATORS:
            d = run_dir(gen, "CrossSection", setup)
            if not os.path.exists(os.path.join(d, "pairs_k_pn.dat")):
                low_k[gen].append(np.nan); continue
            x, y, _, _ = load_hist("pairs_k_pn", d)
            ax.plot(x, y, label=GEN_LABEL[gen], **GEN_STYLE[gen])
            width = x[1] - x[0]
            low_k[gen].append((y[x + 0.5*width <= P_RUN + 1e-9]*width).sum())
        ax.set_yscale("log"); ax.set_xlim(0, 3); ax.set_xlabel(r"$k$ [GeV]"); ax.set_title(setup.replace("_", " ")); ax.legend(fontsize=8)
    axes[0].set_ylabel(r"$pn$ pairs per event per GeV")
    # the last panel: the number of pairs that coalescence can use, per event, for every setup
    ax = axes[4]
    for gen in GENERATORS:
        ax.plot(range(len(setups)), low_k[gen], marker=GEN_MARKER[gen], label=GEN_LABEL[gen], **GEN_STYLE[gen])
    ax.set_xticks(range(len(setups))); ax.set_xticklabels(["Z pole", "0.9 TeV", "2.76 TeV", "7 TeV"])
    ax.set_yscale("log"); ax.set_ylabel(r"$pn$ pairs per event with $k<0.2$ GeV"); ax.set_title("pairs that can coalesce"); ax.legend(fontsize=8)
    for gen in GENERATORS: print(f"    pn pairs with k < 0.2 GeV per event, {gen}:", ", ".join(f"{v:.3g}" for v in low_k[gen]))
    plt.tight_layout(); save(fig, "pairs_k_pn")
    # k of the deuteron-forming pairs: cross-section model (all channels) against HardDisk and Exponential
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4)); any_run = False
    k_gen = "tune4"
    if FINAL:
        k_gen = "mbtune"
    for ax, setup in zip(axes, ["LEP", "ALICE_7TeV"]):
        for model, ls in (("CrossSection", "-"), ("HardDisk_k", "--"), ("Exponential_k", ":")):
            d = run_dir(k_gen, model, setup)
            files = sorted(glob.glob(os.path.join(d, "expected_k_*.dat")))
            if not files: continue
            any_run = True
            tot = None
            for fpath in files:
                x, y, _, _ = load_hist(os.path.basename(fpath)[:-4], d); tot = y if tot is None else tot + y
            ax.plot(x, tot/(tot.sum()*(x[1] - x[0])), ls=ls, color=GEN_STYLE[k_gen]["color"], lw=2, label=f"{GEN_LABEL[k_gen]}, {MODEL_LABEL[model]}")
        ax.set_xlim(0, 2.5); ax.set_xlabel(r"$k$ [GeV]"); ax.set_ylabel("deuteron-forming pairs, normalised"); ax.set_title(setup.replace("_", " ")); ax.legend(fontsize=8)
    plt.tight_layout()
    if any_run: save(fig, "k_of_forming_pairs")
    else: plt.close(fig)

def fig_charged_7TeV():
    fig, axes = plt.subplots(2, 3, figsize=(16, 6.2), sharex="col", gridspec_kw={"height_ratios": [3, 1.3]})
    fnd, _ = f_nd(7000.)
    for col, (title, names, block) in enumerate(((r"$\pi^++\pi^-$", ["meson_pi_plus", "meson_pi_minus"], 0),
                                                  (r"$K^++K^-$", ["meson_K_plus", "meson_K_minus"], 1),
                                                  (r"$p+\bar p$", ["nucleon_p", "nucleon_pbar"], 2))):
        lo, hi, val, err = scorecard.alice_pikp(block); c = 0.5*(lo + hi); ax, rx = axes[0, col], axes[1, col]
        for gen in GENERATORS:
            d = run_dir(gen, "CrossSection", "ALICE_7TeV")
            res = scorecard.spectrum(d, names, block)
            if res is None: continue
            chi2, ratio, _ = res
            ax.plot(c, ratio*val, marker=GEN_MARKER[gen], ms=3, label=rf"{GEN_LABEL[gen]}, $\chi^2$/n = {chi2:.1f}", **GEN_STYLE[gen])
            rx.plot(c, ratio, marker=GEN_MARKER[gen], ms=3, **GEN_STYLE[gen])
        ax.errorbar(c, val, yerr=err, xerr=0.5*(hi - lo), fmt="o", color="black", ms=3.5, label="ALICE", zorder=5)
        rx.fill_between(c, 1 - err/val, 1 + err/val, color="0.85", step="mid")
        ax.set_yscale("log"); ax.set_title(f"{title}, pp 7 TeV, $|y|<0.5$"); ax.legend(fontsize=8)
        rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0.5, 1.5); rx.set_xlabel(r"$p_T$ [GeV]")
    axes[0, 0].set_ylabel(r"$\frac{1}{N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-1}$]"); axes[1, 0].set_ylabel("MC / ALICE")
    plt.tight_layout(); save(fig, "charged_7TeV_pikp")

ALICE_NCH = {"0.9": (3.75, 0.06), "2.76": (4.76, 0.08), "7": (5.98, 0.09)}   # ins1394854 Table 15, INEL>0
def nch_of(d):
    m = np.atleast_2d(np.loadtxt(os.path.join(d, "multiplicity_events.dat"))); return m[:, 2].sum()/m[:, 1].sum()

def fig_nch_energy():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(12, 4.4))
    sq = np.array([0.9, 2.76, 7.])
    al = np.array([ALICE_NCH[r][0] for r in ("0.9", "2.76", "7")]); ae = np.array([ALICE_NCH[r][1] for r in ("0.9", "2.76", "7")])
    ax.errorbar(sq, al, yerr=ae, fmt="o", color="black", ms=5, label="ALICE (INEL>0)", zorder=5)
    bx.plot(sq, al/al[-1], "o", color="black", ms=5, label="ALICE")
    for gen in GENERATORS:
        vals = []
        for rs in ("0.9", "2.76", "7"):
            cands = [os.path.join(run_dir(gen, "CrossSection", f"ALICE_{rs}TeV")), os.path.join(DATA, f"lund_{gen}", "ALICE_7TeV_seeds", "card_seed"),
                     os.path.join(COAL, gen, "HardDisk_k", f"ALICE_{rs}TeV")]
            d = next((c for c in cands if os.path.exists(os.path.join(c, "multiplicity_events.dat")) and (rs in c or "card_seed" in c)), None)
            vals.append(nch_of(d) if d else np.nan)
        vals = np.array(vals)
        if np.all(np.isnan(vals)): continue
        ax.plot(sq, vals, marker=GEN_MARKER[gen], label=f"{GEN_LABEL[gen]} (non-diffractive)", **GEN_STYLE[gen])
        bx.plot(sq, vals/vals[-1], marker=GEN_MARKER[gen], label=GEN_LABEL[gen], **GEN_STYLE[gen])
    for a in (ax, bx): a.set_xscale("log"); a.set_xlabel(r"$\sqrt{s}$ [TeV]"); a.legend(fontsize=8)
    ax.set_ylabel(r"$\langle dN_{\rm ch}/d\eta\rangle_{|\eta|<0.5}$"); ax.set_title("charged-particle density")
    bx.set_ylabel("relative to 7 TeV"); bx.set_title("its energy dependence")
    plt.tight_layout(); save(fig, "nch_energy")

def fig_lep_charged():
    import compare_lep
    compare_lep.LABEL.update({"lund_tune1": "Sherpa Lund_tune1 (Monash strings)", "app_tune": "Sherpa Lund_tune4 (LEP-retuned strings)"})
    compare_lep.STYLE.update({"lund_tune1": GEN_STYLE["tune1"], "app_tune": GEN_STYLE["tune4"]})
    cwd = os.getcwd(); os.chdir(LEP_RIVET)
    try:
        if not os.path.exists("app_tune"): os.symlink("runs_lund_lep1/app_tune", "app_tune")
        with contextlib.redirect_stdout(io.StringIO()): compare_lep.main(FIG, ["lund_tune1", "app_tune"])
    finally: os.chdir(cwd)
    for name in ("lep_nch", "lep_xp"):
        os.remove(os.path.join(FIG, name + ".png")); print("  figures/" + name + ".pdf")

# ------------------------------------------------------------------ impact parameter and coalescence momentum
def fig_impact_parameter():
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.4)); any_run = False
    for ax, setup in zip(axes, ["LEP", "ALICE_7TeV"]):
        for gen in GENERATORS:
            d = run_dir(gen, "HardDisk_kb", setup)
            if not os.path.exists(os.path.join(d, "pairs_b_pn.dat")): continue
            any_run = True
            x, y, _, w = load_hist("pairs_b_pn", d)
            ax.plot(x, y/(y.sum()*w), label=GEN_LABEL[gen] + (" (Pythia vertices on)" if gen != "pythia" else ""), **GEN_STYLE[gen])
        d = run_dir("tune4", "HardDisk_k", setup)
        if os.path.exists(os.path.join(d, "pairs_b_pn.dat")):
            x, y, _, w = load_hist("pairs_b_pn", d)
            ax.plot(x, y/(y.sum()*w), color="0.5", ls="--", lw=1.4, label="Sherpa Lund_tune4, vertices off (default)")
        ax.axvline(3., color="k", ls=":", lw=1, label=r"$b_{\rm coal}$ = 3 fm")
        ax.set_xlim(0, 15); ax.set_yscale("log"); ax.set_xlabel(r"$b$ [fm]"); ax.set_ylabel(r"$pn$ pairs, normalised [fm$^{-1}$]")
        ax.set_title(f"{setup.replace('_', ' ')}: impact parameter of all $pn$ pairs"); ax.legend(fontsize=8)
    plt.tight_layout()
    if any_run: save(fig, "impact_parameter")
    else: plt.close(fig)
    # fraction of the momentum-only HardDisk yield kept by b < b_coal
    fig, ax = plt.subplots(figsize=(7.5, 4.4)); any_run = False
    for setup, ls in (("LEP", "--"), ("ALICE_7TeV", "-")):
        for gen in GENERATORS:
            ref = run_dir(gen, "HardDisk_k", setup)
            if not exists(ref): continue
            tot = lambda d: sum(read_summary(d)[0].get(f"expected_{sp}_total", (0, 0))[0] for sp in ("d", "dbar"))
            bs, fr = [], []
            for b, name in ((2., "HardDisk_kb_b2"), (3., "HardDisk_kb"), (5., "HardDisk_kb_b5")):
                d = run_dir(gen, name, setup)
                if exists(d): bs.append(b); fr.append(tot(d)/tot(ref))
            if not bs: continue
            any_run = True
            ax.plot(bs, fr, marker=GEN_MARKER[gen], color=GEN_STYLE[gen]["color"], ls=ls, label=f"{GEN_LABEL[gen]}, {setup.replace('_', ' ')}")
    ax.set_xlabel(r"$b_{\rm coal}$ [fm]"); ax.set_ylabel(r"yield($b<b_{\rm coal}$) / yield($k$ only)"); ax.set_ylim(0, 1.05)
    ax.set_title(r"HardDisk ($p_0$ = 0.2 GeV): share of the deuterons kept by the $b$ condition"); ax.legend(fontsize=7.5)
    plt.tight_layout()
    if any_run: save(fig, "impact_parameter_fraction")
    else: plt.close(fig)

def fig_pcoal_scan():
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(13, 4.4)); any_run = False
    for setup, ls in (("LEP", "--"), ("ALICE_7TeV", "-")):
        for gen in GENERATORS:
            ps, ys = [], []
            for p, name in ((0.15, "HardDisk_k_p0.15"), (0.2, "HardDisk_k"), (0.3, "HardDisk_k_p0.3")):
                d = run_dir(gen, name, setup)
                if exists(d):
                    ps.append(p); ys.append(sum(read_summary(d)[0].get(f"expected_{sp}_total", (0, 0))[0] for sp in ("d", "dbar")))
            if len(ps) < 2: continue
            any_run = True
            ys = np.array(ys)/ys[ps.index(0.2)] if 0.2 in ps else np.array(ys)
            ax.plot(ps, ys, marker=GEN_MARKER[gen], color=GEN_STYLE[gen]["color"], ls=ls, label=f"{GEN_LABEL[gen]}, {setup.replace('_', ' ')}")
    pp = np.linspace(0.14, 0.32, 20); ax.plot(pp, (pp/0.2)**3, color="k", lw=0.8, label=r"$\propto p_0^3$")
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlabel(r"$p_0$ [GeV]"); ax.set_ylabel(r"deuteron yield / yield at $p_0$ = 0.2 GeV")
    ax.set_title(r"HardDisk ($k$ only): yield against $p_0$"); ax.legend(fontsize=7.5)
    for p, name, ls in ((0.15, "HardDisk_k_p0.15", ":"), (0.2, "HardDisk_k", "-"), (0.3, "HardDisk_k_p0.3", "--")):
        d = run_dir("tune4", name, "ALICE_7TeV")
        if not exists(d): continue
        x, y, _, w = load_hist("expected_d_invyield_midy", d)
        yy = 2*np.pi*x*y; m = x < 4
        bx.plot(x[m], yy[m]/(yy[m].sum()*w), ls=ls, color="#009E73", lw=2, label=rf"Lund_tune4, $p_0$ = {p} GeV")
    bx.set_xlabel(r"$p_T$ [GeV]"); bx.set_ylabel(r"$dN/dp_T$ of $d$ at $|y|<0.5$, normalised"); bx.set_title(r"7 TeV: spectral shape against $p_0$"); bx.legend(fontsize=8)
    plt.tight_layout()
    if any_run: save(fig, "pcoal_scan")
    else: plt.close(fig)

# ------------------------------------------------------------------ summary table (LaTeX)
def summary_table():
    "chi2 of the ALICE spectra (d + dbar, three energies) and B2, and the parameter spread, per model and generator"
    lines = []
    for model in MODELS:
        for gen in GENERATORS:
            chi, chid, n, pars = 0., 0., 0, []
            for rs, sq, t in ENERGIES:
                for sp in ("d", "dbar"):
                    f = alice_fit(gen, model, rs, sq, t, sp)
                    if f: chi += f["chi2"]; chid += f["chi2_data"]; n += f["n"]; pars.append(parameter(model, f["s"], f["se"])[0])
            if n == 0: continue
            lf, bf = lep_fit(gen, model), babar_fit(gen, model)
            lep = parameter(model, lf["s"], lf["se"])[0] if lf else np.nan
            bab = parameter(model, bf["s"], bf["se"])[0] if bf else np.nan
            unit = 1. if model == "CrossSection" else 1e3
            fmt = lambda v: "--" if np.isnan(v) else (f"{v*unit:.2f}" if model == "CrossSection" else f"{v*unit:.0f}")
            lines.append(f"{MODEL_LABEL[model]} & {GEN_LABEL[gen].replace('_', r'\_')} & {chid:.0f} & {chi:.0f}/{n} & {fmt(min(pars))}--{fmt(max(pars))} & {fmt(lep)} & {fmt(bab)} \\\\")
    with open(TABLE, "w") as f:
        f.write("% written by make_figures.py\n" + "\n".join(lines) + "\n")
    print("  summary_table.tex")

if __name__ == "__main__" and FINAL:
    # the coalescence runs exist for Lund_tune1 and Pythia only, not for the final setting:
    # cross-section model only
    os.makedirs(FIG, exist_ok=True)
    MODELS = ["CrossSection"]
    print("inputs"); fig_sigma_k(); fig_pairs_k(); fig_charged_7TeV(); fig_nch_energy()
    for model in MODELS:
        print(model)
        fig_alice_spectra(model); fig_b2(model); fig_lep_spectrum(model); fig_consistency(model); fig_multiplicity(model)
    summary_table()
elif __name__ == "__main__":
    os.makedirs(FIG, exist_ok=True)
    print("inputs"); fig_sigma_k(); fig_pairs_k(); fig_charged_7TeV(); fig_nch_energy(); fig_lep_charged()
    for model in MODELS:
        print(model)
        fig_alice_spectra(model); fig_b2(model); fig_lep_spectrum(model); fig_consistency(model); fig_multiplicity(model)
    print("impact parameter, p_0"); fig_impact_parameter(); fig_pcoal_scan()
    summary_table()
