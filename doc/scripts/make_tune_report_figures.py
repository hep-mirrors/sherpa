"""Figures and numbers of deuteron_tune_report.tex: the two Sherpa settings s6_008 and s6_010
against the data and Pythia 8.318, for the plots of Dal & Raklev (arXiv:1504.07242) and the
further measurements.

  python3 make_figures.py

Every figure has one column per Sherpa setting; each column shows the measurement, Pythia
8.318 and that setting.  Writes figures/*.pdf and numbers.json (the numbers quoted in the
text), into doc/deuteron_tune_report/.  Reads only generator output, through deuteron_common.py, so the numbers
are those of the notebooks in SimpleLEP/XSecs."""
import os, sys, glob, json
os.environ["MPLBACKEND"] = "Agg"
from paths import XSECS, DEUTERON_TUNE_REPORT
os.chdir(XSECS)    # the notebooks' folder, as when they run
from deuteron_common import *
import matplotlib.ticker

FIG = os.path.join(DEUTERON_TUNE_REPORT, "figures")
os.makedirs(FIG, exist_ok=True)
SETTINGS = ["Sherpa s6_008", "Sherpa s6_010"]
PY = "Pythia 8.318"
LOOK = {"Sherpa s6_008": dict(color="#D55E00", ls="-", lw=1.8),
        "Sherpa s6_010": dict(color="#009E73", ls="-", lw=1.8),
        PY:              dict(color="#0072B2", ls="--", lw=1.6)}
STEP = dict(baseline=None)    # step curves without vertical lines down to zero at the ends
PAPER = dict(color="0.45", mfc="none", ls="none")
DATA = dict(color="black", fmt="o", ms=3.5, elinewidth=1.0, zorder=5)
INV_LABEL = r"$\frac{1}{2\pi p_T N_{\rm inel}}\frac{d^2N}{dy\,dp_T}$ [GeV$^{-2}$]"
NUM = {}

def short(m): return m.replace("Sherpa ", "")
def save(fig, name):
    fig.savefig(os.path.join(FIG, name + ".pdf")); plt.close(fig); print("wrote", name)

def edges_of(lo, hi):
    "bin edges of contiguous bins"
    return np.append(lo, hi[-1])

def data_band(rx, lo, hi, val, err):
    "grey band of the relative data error in the ratio panel"
    for l, h, v, e in zip(lo, hi, val, err):
        rx.fill_between([l, h], 1 - e/v, 1 + e/v, color="0.85", lw=0)
    rx.axhline(1, color="k", lw=0.8)

def binned_curve(ax, rx, lo, hi, val, mc, model, label):
    "a model averaged over the data bins: steps in the main and the ratio panel"
    ax.stairs(mc, edges_of(lo, hi), label=label, **STEP, **LOOK[model])
    rx.stairs(mc/val, edges_of(lo, hi), **STEP, **LOOK[model])

# ------------------------------------------------------------------------------------------
# paper Figs. 4, 5: the cross sections (the same in every setting)
# ------------------------------------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
groups = [["pn_dgamma", "pn_dpi0", "pn_dpi0pi0", "pn_dpipluspiminus"],
          ["pp_dpiplus", "pp_dpipluspi0", "nn_dpiminus", "nn_dpiminuspi0"]]
for ax, group, title in zip(axes, groups, [r"$np\to dX$ (cf. Fig. 4)", r"$pp\to dX$, $nn\to dX$ (cf. Fig. 5)"]):
    for ch in group:
        x, y, _, _ = load_hist("test_sigma_" + ch, DATA_DIR)
        ax.plot(x, y, label=LABELS[ch])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(0.1, 4.); ax.set_ylim(1., 5e3)
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"$k$ [GeV]"); ax.set_ylabel(r"$\sigma$ [$\mu$b]"); ax.set_title(title); ax.legend(fontsize=8)
plt.tight_layout(); save(fig, "fig_sigma")

# ------------------------------------------------------------------------------------------
# paper Fig. 7: ALICE (anti)deuteron spectra at 0.9, 2.76, 7 TeV
# ------------------------------------------------------------------------------------------
FITS = {}
for m in SETTINGS + [PY]:
    for rs, sqrts, table in ENERGIES:
        fnd, _ = f_nd(sqrts); setup = f"ALICE_{rs}TeV"
        for sp, block in (("d", 0), ("dbar", 1)):
            lo, hi, val, err = alice_spectrum(table, block)
            mc, _ = in_bins(deuteron_hist(m, setup, sp, "invyield_midy"), lo, hi, scale=fnd)
            x, xe, chi2 = fit_scale(val, err, mc)
            FITS[(m, rs, sp)] = dict(inv=x, err=xe, chi2=chi2, n=len(val), mc=x*mc, lo=lo, hi=hi, val=val, dval=err)
NUM["alice_fits"] = {f"{short(m)}|{rs}|{sp}": [round(v["inv"], 3), round(v["err"], 3), round(v["chi2"], 1), v["n"]]
                     for (m, rs, sp), v in FITS.items()}

for sp, lab in (("dbar", r"\bar d"), ("d", "d")):
    fig, axes = plt.subplots(6, 2, figsize=(11, 15), gridspec_kw={"height_ratios": [3, 1]*3, "hspace": 0.0})
    for row, (rs, sqrts, table) in enumerate(ENERGIES):
        for col, m in enumerate(SETTINGS):
            ax, rx = axes[2*row, col], axes[2*row + 1, col]
            f0 = FITS[(m, rs, sp)]; lo, hi, val, err = f0["lo"], f0["hi"], f0["val"], f0["dval"]
            ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE", **DATA)
            data_band(rx, lo, hi, val, err)
            for mm in (PY, m):
                f = FITS[(mm, rs, sp)]
                binned_curve(ax, rx, lo, hi, val, f["mc"], mm,
                             f"{mm}: $1/\\sigma_0$ = {f['inv']:.2f} b$^{{-1}}$, $\\chi^2$ = {f['chi2']:.1f}/{f['n']}")
            ax.set_yscale("log"); ax.set_xlim(lo[0] - 0.05, hi[-1] + 0.05); rx.set_xlim(ax.get_xlim())
            ax.set_xticklabels([]); ax.legend(fontsize=7.5, loc="lower left")
            ax.text(0.97, 0.92, f"${lab}$, pp {rs} TeV, $|y|<0.5$", transform=ax.transAxes, ha="right", va="top", fontsize=10)
            rx.set_ylim(0.3, 1.7); rx.set_ylabel("MC/data", fontsize=9)
            if col == 0: ax.set_ylabel(INV_LABEL, fontsize=10)
            if row == 0: ax.set_title(f"ALICE, Pythia and {m}")
            if row == 2: rx.set_xlabel(r"$p_T$ [GeV]")
    fig.subplots_adjust(hspace=0.0, wspace=0.25)
    for row in range(3):                                   # gap between the energies
        for col in range(2):
            for k, a in enumerate((axes[2*row, col], axes[2*row+1, col])):
                p = a.get_position(); a.set_position([p.x0, p.y0 - 0.025*row, p.width, p.height])
    save(fig, f"fig7_{sp}")

# MC/ALICE of the 7 TeV deuterons in 1-GeV ranges
NUM["alice7_ratio"] = {}
for m in SETTINGS + [PY]:
    for sp in ("d", "dbar"):
        f = FITS[(m, "7", sp)]; c = 0.5*(f["lo"] + f["hi"]); r = f["mc"]/f["val"]
        NUM["alice7_ratio"][f"{short(m)}|{sp}"] = [round(r[(c >= a) & (c < a+1)].mean(), 2) for a in (0, 1, 2)]

# ------------------------------------------------------------------------------------------
# paper Fig. 8: protons at 7 TeV, and the pions and kaons of the same measurement
# ------------------------------------------------------------------------------------------
f7, _ = f_nd(7000.)
blocks = read_hepdata_blocks(os.path.join(HEPDATA, "ALICE_pp7TeV_piKp_HEPData_ins1357424_Table1.csv"))
def alice_pikp(block, ptmax):
    r_ = [r for r in blocks[block] if r[3] is not None and r[2] <= ptmax + 1e-9]
    lo, hi, val = (np.array([r[i] for r in r_]) for i in (1, 2, 3)); err = np.hypot([r[4] for r in r_], [r[6] for r in r_])
    pt = 0.5*(lo+hi)
    return lo, hi, val/(2*np.pi*pt), err/(2*np.pi*pt)

def mc_sum(m, names):
    "f_ND x sum of the histograms `names` of model m at 7 TeV (x, y, e, w), or None"
    d = model_path(m, "ALICE_7TeV"); tot, e2 = None, None
    for n in names:
        if not os.path.exists(os.path.join(d, n + ".dat")): return None
        x, y, e, w = load_hist(n, d)
        tot = y if tot is None else tot + y; e2 = e**2 if e2 is None else e2 + e**2
    return x, f7*tot, f7*np.sqrt(e2), w

def ratio_ranges(c, pt, inv, ranges):
    "MC/ALICE averaged over pT ranges, as in ALICE_2018_pp.ipynb"
    x, t = c[0], c[1]; out = []
    for a, b in ranges:
        ma, mm = (pt > a) & (pt < b), (x > a - 1e-9) & (x < b - 1e-9)
        out.append(round(t[mm].mean()/inv[ma].mean(), 2))
    return out

SPECIES = {"protons": (2, ("nucleon_p_invyield_midy", "nucleon_pbar_invyield_midy"), r"$p+\bar p$",
                       [(0.3, 0.6), (0.6, 1), (1, 2), (2, 3), (3, 3.5), (3.5, 4)], 4.0),
           "pions":   (0, ("meson_pi_plus_invyield_midy", "meson_pi_minus_invyield_midy"), r"$\pi^++\pi^-$",
                       [(0.1, 0.5), (0.5, 1), (1, 2), (2, 3)], 3.0),
           "kaons":   (1, ("meson_K_plus_invyield_midy", "meson_K_minus_invyield_midy"), r"$K^++K^-$",
                       [(0.2, 0.5), (0.5, 1), (1, 2), (2, 3), (3, 4)], 4.0)}
NUM["pikp_ratio"] = {}
for name, (block, hists, lab, ranges, ptmax) in SPECIES.items():
    lo, hi, inv, inverr = alice_pikp(block, ptmax); pt = 0.5*(lo+hi)
    for m in SETTINGS + [PY]:
        c = mc_sum(m, hists)
        if c is not None: NUM["pikp_ratio"][f"{name}|{short(m)}"] = ratio_ranges(c, pt, inv, ranges)
    NUM["pikp_ranges|" + name] = ranges

def pikp_figure(names, fname):
    fig, axes = plt.subplots(2*len(names), 2, figsize=(11, 5.6*len(names)),
                             gridspec_kw={"height_ratios": [3, 1]*len(names)})
    for row, name in enumerate(names):
        block, hists, lab, ranges, ptmax = SPECIES[name]
        lo, hi, inv, inverr = alice_pikp(block, ptmax)
        for col, m in enumerate(SETTINGS):
            ax, rx = axes[2*row, col], axes[2*row + 1, col]
            ax.errorbar(0.5*(lo+hi), inv, yerr=inverr, xerr=0.5*(hi-lo), label=f"ALICE {lab}", **dict(DATA, ms=3))
            data_band(rx, lo, hi, inv, inverr)
            for mm in (PY, m):
                c = mc_sum(mm, hists)
                if c is None: continue
                mc, _ = in_bins(c, lo, hi)
                binned_curve(ax, rx, lo, hi, inv, mc, mm, mm)
            ax.set_yscale("log"); ax.set_xlim(0, ptmax); rx.set_xlim(0, ptmax); ax.set_xticklabels([])
            ax.legend(fontsize=8, loc="lower left"); rx.set_ylim(0.4, 1.8); rx.set_ylabel("MC/data", fontsize=9)
            rx.set_xlabel(r"$p_T$ [GeV]")
            ax.set_title(f"{lab}, pp 7 TeV, $|y|<0.5$: ALICE, Pythia, {m}" if name == "protons" else
                         f"{lab}, pp 7 TeV, $|y|<0.5$: ALICE, {m}", fontsize=10)
            if col == 0: ax.set_ylabel(INV_LABEL, fontsize=10)
    plt.tight_layout(); save(fig, fname)
pikp_figure(["protons"], "fig8_protons")
pikp_figure(["pions", "kaons"], "fig_pions_kaons")

# ------------------------------------------------------------------------------------------
# ALICE: B2(pT/A) at the three energies, and dbar/d
# ------------------------------------------------------------------------------------------
fig, axes = plt.subplots(2, 3, figsize=(15, 8.2), sharey=True)
NUM["alice_B2"] = {}
for col, (rs, sqrts, table) in enumerate(ENERGIES):
    fnd, _ = f_nd(sqrts); setup = f"ALICE_{rs}TeV"
    blks = read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp_deuterons_HEPData_ins1625294_Table{5+table}.csv"))
    r_ = [r for r in blks[0] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in r_]) for i in (1, 2, 3))
    err = np.hypot(0.5*(np.array([r[4] for r in r_])-np.array([r[5] for r in r_])), 0.5*(np.array([r[6] for r in r_])-np.array([r[7] for r in r_])))
    for row, m in enumerate(SETTINGS):
        ax = axes[row, col]
        ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label=r"ALICE ($d$)", **DATA)
        for mm in (PY, m):
            dv, _ = in_bins(deuteron_hist(mm, setup, "d", "invyield_midy"), 2*lo, 2*hi)
            # protons at the bin centre of pT/A, interpolated (log) in the 0.1 GeV histogram: the
            # B2 bins at 7 TeV are 0.05 GeV wide, half a histogram bin
            xp, yp, _, _ = load_hist("nucleon_p_invyield_midy", model_path(mm, setup))
            pv = np.exp(np.interp(0.5*(lo+hi), xp[yp > 0], np.log(yp[yp > 0])))
            b2 = FITS[(mm, rs, "d")]["inv"]*dv/(fnd*pv**2)
            ax.stairs(b2, edges_of(lo, hi), label=mm, **STEP, **LOOK[mm])
            NUM["alice_B2"][f"{short(mm)}|{rs}"] = [round(v*1e3, 2) for v in b2]
        NUM["alice_B2"][f"ALICE|{rs}"] = [round(v*1e3, 2) for v in val]
        ax.set_title(f"$B_2$, pp {rs} TeV: ALICE, Pythia, {m}", fontsize=10); ax.legend(fontsize=8, loc="upper left")
        if row == 1: ax.set_xlabel(r"$p_T/A$ [GeV]")
        if col == 0: ax.set_ylabel(r"$B_2$ [GeV$^2$]")
        ax.set_ylim(0, 0.05)
plt.tight_layout(); save(fig, "fig_alice_B2")

def dndy(val, lo, hi): return (2*np.pi*0.5*(lo+hi)*val*(hi-lo)).sum()
RAT = {}
for rs, sqrts, table in ENERGIES:
    fnd, _ = f_nd(sqrts); setup = f"ALICE_{rs}TeV"
    for sp, block in (("d", 0), ("dbar", 1)):
        lo, hi, val, err = alice_spectrum(table, block)
        RAT[("ALICE", rs, sp)] = (dndy(val, lo, hi), np.sqrt(((2*np.pi*0.5*(lo+hi)*err*(hi-lo))**2).sum()))
        for m in SETTINGS + [PY]:
            mc, _ = in_bins(deuteron_hist(m, setup, sp, "invyield_midy"), lo, hi, scale=fnd)
            RAT[(m, rs, sp)] = (dndy(mc, lo, hi), 0.)
E = np.array([0.9, 2.76, 7.])
ratio = lambda who, rs: RAT[(who, rs, "dbar")][0]/RAT[(who, rs, "d")][0]
rerr = lambda rs: ratio("ALICE", rs)*np.hypot(RAT[("ALICE", rs, "dbar")][1]/RAT[("ALICE", rs, "dbar")][0],
                                              RAT[("ALICE", rs, "d")][1]/RAT[("ALICE", rs, "d")][0])
NUM["dbar_over_d"] = {who if who == "ALICE" else short(who): [round(ratio(who, rs), 3) for rs, _, _ in ENERGIES]
                      for who in ["ALICE"] + SETTINGS + [PY]}
NUM["dbar_over_d_err"] = [round(rerr(rs), 3) for rs, _, _ in ENERGIES]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, m in zip(axes, SETTINGS):
    ax.errorbar(E, [ratio("ALICE", rs) for rs, _, _ in ENERGIES], yerr=[rerr(rs) for rs, _, _ in ENERGIES], label="ALICE", **DATA)
    for k, mm in enumerate((PY, m)):
        ax.plot(E*(1 + 0.04*(k+1)), [ratio(mm, rs) for rs, _, _ in ENERGIES], marker="s", ms=5, label=mm, **LOOK[mm])
    ax.axhline(1, color="0.5", lw=0.8); ax.set_xscale("log"); ax.set_ylim(0.6, 1.2)
    ax.set_xlabel(r"$\sqrt{s}$ [TeV]"); ax.set_title(f"$\\bar d/d$ (measured $p_T$ range): ALICE, Pythia, {m}", fontsize=10)
    ax.legend(fontsize=8)
axes[0].set_ylabel(r"$\bar d/d$")
plt.tight_layout(); save(fig, "fig_dbar_over_d")

# ------------------------------------------------------------------------------------------
# LEP (ALEPH, OPAL): the fit of the paper, the spectrum, B2
# ------------------------------------------------------------------------------------------
A, sA = 5.9e-6, np.hypot(1.8e-6, 0.5e-6)
eps, Nev = 0.234, 1.64e6
LEP = {}
for m in SETTINGS + [PY]:
    a, ae = deuteron_count(m, "LEP", "dbar_ALEPH_window")
    b, be = deuteron_count(m, "LEP", "dbar_OPAL_window")
    x = max(0., (A - eps*Nev*b*sA**2/(2*a))/a); xe = sA/a
    chi2 = (A - a*x)**2/sA**2 + eps*Nev*b*x
    s = read_summary(model_path(m, "LEP"))[0]
    npb = s.get("nucleon_pbar_ALEPH_halfwindow", (0., 0.))[0]
    dOmega = 2*np.pi*2*0.95
    inv_w = lambda n, mass, p1, p2: n*np.sqrt((0.5*(p1+p2))**2 + mass**2)/((0.5*(p1+p2))**2*dOmega*(p2-p1))
    B2 = inv_w(a*x, 1.875613, 0.62, 1.03)/inv_w(npb, 0.938272, 0.31, 0.515)**2
    LEP[m] = dict(a=a, b=b, inv=x, err=xe, chi2=chi2, aleph_only=A/a, B2=B2, pbar_half=npb, opal=eps*Nev*b*x)
NUM["lep"] = {short(m): {k: (float(f"{v:.4g}")) for k, v in d.items()} for m, d in LEP.items()}

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
edges = np.arange(0., 3.001, 0.1); cen = 0.5*(edges[1:] + edges[:-1])
for ax, m in zip(axes, SETTINGS):
    ax.axvspan(0.62, 1.03, color="grey", alpha=0.12)
    ax.errorbar([0.825], [A/0.41], yerr=[sA/0.41], xerr=[0.205], label=r"ALEPH, $|\cos\theta|<0.95$", **DATA)
    for mm in (PY, m):
        xs, y, e, w = deuteron_hist(mm, "LEP", "dbar", "p")
        v = np.array([(y[(xs > l) & (xs < h)]*w).sum()/(h-l) for l, h in zip(edges[:-1], edges[1:])])
        L = LEP[mm]
        ax.stairs(L["inv"]*v, edges, label=f"{mm}: $1/\\sigma_0$ = {L['inv']:.2f} b$^{{-1}}$", **STEP, **LOOK[mm])
        ax.plot([0.825], [L["a"]*L["inv"]/0.41], marker="s", ms=6, mfc="none", color=LOOK[mm]["color"], ls="none")
    ax.set_xlim(0, 3); ax.set_ylim(0, None); ax.set_xlabel(r"$p$ [GeV]")
    ax.set_title(f"$\\bar d$ in Z decays: ALEPH, Pythia, {m}", fontsize=10); ax.legend(fontsize=8)
axes[0].set_ylabel(r"$\bar d$ per hadronic Z decay per GeV")
plt.tight_layout(); save(fig, "fig_lep")

# ------------------------------------------------------------------------------------------
# paper Fig. 9: BaBar continuum
# ------------------------------------------------------------------------------------------
row = read_hepdata_blocks(os.path.join(HEPDATA, "BaBar_continuum_ratio_HEPData_ins1286317_Table5.csv"))[0][0]
R, R_err = row[1], np.hypot(0.5*(row[2]-row[3]), 0.5*(row[4]-row[5]))
sig_had = 9.63/R
BAB = {}
for m in SETTINGS + [PY]:
    a, ae = deuteron_count(m, "BaBar_10.58GeV", "dbar_total")
    BAB[m] = dict(a=a, inv=R/a, err=R_err/a, rate_at_lep=a*LEP[m]["inv"], xsec_at_lep=a*LEP[m]["inv"]*sig_had)
NUM["babar"] = {short(m): {k: float(f"{v:.4g}") for k, v in d.items()} for m, d in BAB.items()}
NUM["babar_rate"] = [R, R_err]
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
edges = np.arange(0.35, 2.2501, 0.1)
for ax, m in zip(axes, SETTINGS):
    for mm in (PY, m):
        xs, y, e, w = deuteron_hist(mm, "BaBar_10.58GeV", "dbar", "p")
        v = np.array([(y[(xs > l) & (xs < h)]*w).sum()/(h-l) for l, h in zip(edges[:-1], edges[1:])])
        ax.stairs(BAB[mm]["inv"]*v, edges, label=f"{mm}: $1/\\sigma_0$ = {BAB[mm]['inv']:.2f} b$^{{-1}}$", **STEP, **LOOK[mm])
    ax.text(0.97, 0.95, "BaBar: $(3.01\\pm0.36)\\times10^{-6}$ $\\bar d$ per event\n(all momenta, from BaBar's fit)",
            transform=ax.transAxes, ha="right", va="top", fontsize=8)
    ax.set_xlim(0.35, 2.25); ax.set_ylim(0, None); ax.set_xlabel(r"$p$ [GeV]")
    ax.set_title(f"$e^+e^-\\to q\\bar q$, 10.58 GeV: Pythia, {m}", fontsize=10); ax.legend(fontsize=8, loc="center right")
axes[0].set_ylabel(r"$\frac{1}{N_{\rm had}}\frac{dN_{\bar d}}{dp}$ [GeV$^{-1}$]")
plt.tight_layout(); save(fig, "fig9_babar")

# ------------------------------------------------------------------------------------------
# paper Fig. 6: antinucleon pairs at the Z pole; k distributions against Pythia
# ------------------------------------------------------------------------------------------
def load_map(d, name):
    a = np.loadtxt(os.path.join(d, "pairmap_" + name + ".dat"))
    n = int(round(np.sqrt(len(a)))); mm = a[:, 2].reshape(n, n)
    return mm.reshape(n//2, 2, n//2, 2).sum(axis=(1, 3))
fig, axes = plt.subplots(2, 3, figsize=(14, 8))
NUM["pairs_ratio"] = {}
for row, m in enumerate(SETTINGS):
    d = model_path(m, "LEP")
    maps = {n: load_map(d, n) for n in ("pbarnbar", "pbarpbar", "nbarnbar")}
    vmax = max(v.max() for v in maps.values())
    for ax, (n, lab) in zip(axes[row], [("pbarnbar", r"$\bar p\bar n$"), ("pbarpbar", r"$\bar p\bar p$"), ("nbarnbar", r"$\bar n\bar n$")]):
        im = ax.imshow(maps[n]/vmax, origin="lower", extent=[0, 5, 0, 5], vmin=0, vmax=1, cmap="YlGnBu_r", aspect="auto")
        ax.set_xlabel(r"$|\vec p_1+\vec p_2|$ [GeV]"); ax.set_ylabel(r"$k$ [GeV]"); ax.set_title(f"{m}: {lab} pairs", fontsize=10)
        ax.grid(False); ax.minorticks_off(); fig.colorbar(im, ax=ax)
    tot = {n: v.sum() for n, v in maps.items()}
    NUM["pairs_ratio"][short(m)] = [round(tot["pbarnbar"]/tot["pbarpbar"], 2), round(tot["pbarnbar"]/tot["nbarnbar"], 2)]
plt.tight_layout(); save(fig, "fig6_pairmaps")

fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
PAIRC = {"pn": "#0072B2", "pp": "#D55E00", "nn": "#009E73"}
NUM["pairs_per_event"] = {}
for ax, m in zip(axes, SETTINGS):
    for mm, ls in ((PY, "--"), (m, "-")):
        for pair in ("pn", "pp", "nn"):
            x, y, e, w = load_hist("pairs_k_" + pair, model_path(mm, "LEP"))
            ax.plot(x, y, color=PAIRC[pair], ls=ls, lw=1.5, label=f"{pair}, {mm}")
            NUM["pairs_per_event"][f"{short(mm)}|{pair}"] = round(float((y*w).sum()), 3)
            NUM["pairs_per_event"][f"{short(mm)}|{pair}|k<1.5"] = round(float((y[x < 1.5]*w).sum()), 4)
    ax.axvline(1.0, color="k", ls=":", lw=1)
    ax.set_xlim(0, 5); ax.set_xlabel(r"$k$ [GeV]"); ax.set_title(f"nucleon + antinucleon pairs, Z pole: Pythia, {m}", fontsize=10)
    ax.legend(fontsize=7.5, ncol=2)
axes[0].set_ylabel("pairs per event per GeV")
plt.tight_layout(); save(fig, "fig_pairs_k")

NUM["channels_lep"] = {}
for m in SETTINGS:
    s = read_summary(model_path(m, "LEP"))[0]
    tot = sum(s.get("expected_dbar_" + ch, (0, 0))[0] for ch in CHANNELS)
    NUM["channels_lep"][short(m)] = {ch: round(s.get("expected_dbar_" + ch, (0, 0))[0]/tot, 3) for ch in CHANNELS}

# ------------------------------------------------------------------------------------------
# the deuteron partners of the protons at 7 TeV (Sherpa vs Pythia)
# ------------------------------------------------------------------------------------------
PARTNER = {m: model_path(m, "ALICE_7TeV") for m in SETTINGS}
PARTNER[PY] = os.path.join(DATA_DIR, "pythia_anti", "ALICE_7TeV")
def raw(d, name):
    rows = [l.split() for l in open(os.path.join(d, name + ".dat"))
            if l.strip() and not l.startswith("#") and len(l.split()) in (2, 3)]
    a = np.array([[float(r[0]), float(r[1])] for r in rows])[:-1]
    nev = nevents(read_summary(d)[1])
    return a[:, 0], a[:, 1]*(a[1, 0] - a[0, 0])*nev
PBINS = [(0, 0.6), (0.6, 1), (1, 1.6), (1.6, 2), (2, 3), (3, 5)]
NUM["partners"] = {}
fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
for ax, m in zip(axes, SETTINGS):
    for mm in (PY, m):
        x, P = raw(PARTNER[mm], "partner_protons_pT"); _, K = raw(PARTNER[mm], "partner_protons_k05_pT")
        frac, ferr = [], []
        for a, b in PBINS:
            sel = (x >= a - 1e-9) & (x < b - 1e-9)
            frac.append(100*K[sel].sum()/P[sel].sum()); ferr.append(100*np.sqrt(max(K[sel].sum(), 1.))/P[sel].sum())
        NUM["partners"][short(mm)] = [round(f, 2) for f in frac]
        ax.stairs(frac, [b[0] for b in PBINS] + [PBINS[-1][1]], label=mm, **STEP, **LOOK[mm])
        ax.errorbar([0.5*(a+b) for a, b in PBINS], frac, yerr=ferr, fmt="none", color=LOOK[mm]["color"])
    ax.set_xlabel(r"proton $p_T$ [GeV]"); ax.set_xlim(0, 5); ax.set_ylim(0, None); ax.legend(fontsize=8)
    ax.set_title(f"protons with a neutron within $k<0.5$ GeV, 7 TeV: Pythia, {m}", fontsize=9)
axes[0].set_ylabel("fraction of protons [%]")
plt.tight_layout(); save(fig, "fig_partners")

# ------------------------------------------------------------------------------------------
# ALICE 2019: multiplicity classes at 7 TeV (Sherpa only: the Pythia runs do not record the V0 multiplicity)
# ------------------------------------------------------------------------------------------
CLASSES = [("I+II", 0., 4.7, 17.47, 0.52), ("III", 4.7, 9.5, 13.50, 0.40), ("IV+V", 9.5, 19., 10.76, 0.30),
           ("VI+VII", 19., 38., 7.54, 0.23), ("VIII+IX+X", 38., 100., 3.30, 0.13)]
MSPEC = ["d", "dbar", "p", "pbar"]
def load_rows(path): return {int(r[0]): np.array(r[1:]) for r in np.atleast_2d(np.loadtxt(path))}
def build_classes(base):
    seeds = sorted(glob.glob(os.path.join(base, "seed_*")))
    ev, yl = {}, {sp: {} for sp in MSPEC}
    for sd in seeds:
        for k, v in load_rows(os.path.join(sd, "multiplicity_events.dat")).items(): ev[k] = ev.get(k, 0) + v
        for sp in MSPEC:
            for k, v in load_rows(os.path.join(sd, f"multiplicity_{sp}_invyield_midy.dat")).items(): yl[sp][k] = yl[sp].get(k, 0) + v
    nv0 = np.array(sorted(ev)); nev = np.array([ev[k][0] for k in nv0])
    order = np.argsort(-nv0); cum_hi = np.cumsum(nev[order])/nev.sum()*100; cum_lo = cum_hi - nev[order]/nev.sum()*100
    out = []
    for name, plo, phi, anch, anche in CLASSES:
        w = np.zeros(len(nv0)); w[order] = np.clip(np.minimum(cum_hi, phi) - np.maximum(cum_lo, plo), 0, None)/np.maximum(cum_hi-cum_lo, 1e-12)
        N = (w*nev).sum()
        c = {"name": name, "nch": (w*np.array([ev[k][1] for k in nv0])).sum()/N}
        for j, sp in enumerate(MSPEC):
            c[sp+"_dndy"] = (w*np.array([ev[k][2+j] for k in nv0])).sum()/N
            c[sp] = (w[:, None]*np.array([yl[sp][k] for k in nv0])).sum(axis=0)/N
        out.append(c)
    return out, inv_sigma0(read_summary(seeds[0])[1])
pt_c = np.arange(0.05, 5.0, 0.1)
def class_dndpt(c, sp, l, h):
    sel = (pt_c > l) & (pt_c < h)
    return (2*np.pi*pt_c[sel]*c[sp][sel]*0.1).sum()/(h-l)
MDATA = []
for ic in range(len(CLASSES)):
    b = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table{ic+1}.csv"))[0] if r[3] is not None]
    lo, hi, val = (np.array([r[i] for r in b]) for i in (1, 2, 3))
    err = np.hypot(0.5*(np.array([r[4] for r in b])-np.array([r[5] for r in b])), 0.5*(np.array([r[6] for r in b])-np.array([r[7] for r in b])))
    MDATA.append((lo, hi, val, err))
MULT = {}
fig, axes = plt.subplots(2, 5, figsize=(19, 7.5), sharey=True)
for row, m in enumerate(SETTINGS):
    cls, inv0 = build_classes(os.path.join(model_path(m, "ALICE_7TeV_multiplicity")))
    mcs = [np.array([class_dndpt(c, "d", l, h) for l, h in zip(lo, hi)]) for c, (lo, hi, val, err) in zip(cls, MDATA)]
    s, _, _ = fit_scale(np.concatenate([d[2] for d in MDATA]), np.concatenate([d[3] for d in MDATA]), np.concatenate(mcs))
    chis = []
    for ax, c, (lo, hi, val, err), mc, cl in zip(axes[row], cls, MDATA, mcs, CLASSES):
        chi2 = (((val - s*mc)/err)**2).sum(); chis.append(round(chi2, 1))
        ax.errorbar(0.5*(lo+hi), val, yerr=err, xerr=0.5*(hi-lo), label="ALICE", **DATA)
        ax.stairs(s*mc, edges_of(lo, hi), label=f"{m}, $\\chi^2$ = {chi2:.0f}/{len(val)}", **STEP, **LOOK[m])
        ax.set_yscale("log"); ax.set_title(f"class {cl[0]}", fontsize=10); ax.legend(fontsize=7.5, loc="lower left")
        if row == 1: ax.set_xlabel(r"$p_T$ [GeV]")
    axes[row, 0].set_ylabel(r"$\frac{1}{N_{ev}}\frac{d^2N_d}{dy\,dp_T}$ [GeV$^{-1}$]")
    MULT[m] = dict(cls=cls, s=s, inv=s*inv0, chi2=chis)
plt.tight_layout(); save(fig, "fig_mult_spectra")

t12 = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, "ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table12.csv"))[0] if r[3] is not None]
b2data = []
for ic in range(len(CLASSES)):
    r_ = [r for r in read_hepdata_blocks(os.path.join(HEPDATA, f"ALICE_pp7TeV_multiplicity_HEPData_ins1721729_Table{6+ic}.csv"))[0] if r[3] is not None]
    sel = [r for r in r_ if r[1] >= 0.7-1e-9 and r[2] <= 0.8+1e-9]
    if sel: b2data.append((CLASSES[ic][3], np.mean([r[3] for r in sel]), np.mean([np.hypot(r[4], r[6]) for r in sel])))
NUM["mult"] = {"alice_nch": [c[3] for c in CLASSES], "alice_dp": [r[3] for r in t12],
               "alice_B2": [round(b[1]*1e3, 2) for b in b2data]}
fig, axes = plt.subplots(2, 3, figsize=(15, 8))
for row, m in enumerate(SETTINGS):
    a1, a2, a3 = axes[row]
    a1.errorbar([r[0] for r in t12], [r[3] for r in t12], yerr=[np.hypot(r[4], r[6]) for r in t12], label="ALICE", **DATA)
    a2.plot([r[0] for r in t12], np.array([r[3] for r in t12])/t12[0][3], "o", color="black", label="ALICE")
    a3.errorbar([b[0] for b in b2data], [b[1] for b in b2data], yerr=[b[2] for b in b2data], label="ALICE", **DATA)
    cls, s = MULT[m]["cls"], MULT[m]["s"]; nch = np.array([c["nch"] for c in cls])
    dp = np.array([2*s*c["d_dndy"]/(c["p_dndy"]+c["pbar_dndy"]) for c in cls])
    md_, mp_ = (pt_c > 1.4) & (pt_c < 1.6), (pt_c > 0.7) & (pt_c < 0.8)
    b2 = [s*0.5*(c["d"][md_]+c["dbar"][md_]).mean()/(0.5*(c["p"][mp_]+c["pbar"][mp_]).mean())**2 for c in cls]
    a1.plot(nch, dp, marker="s", label=m, **LOOK[m]); a2.plot(nch, dp/dp[-1], marker="s", label=m, **LOOK[m])
    a3.plot(nch, b2, marker="s", label=m, **LOOK[m])
    MULT[m].update(nch=[round(v, 2) for v in nch], dp=[float(f"{v:.3g}") for v in dp], b2=[round(v*1e3, 2) for v in b2])
    for a, yl, t in [(a1, r"$2d/(p+\bar p)$", r"$d/p$ vs multiplicity"), (a2, "relative to the lowest class", r"trend only (no $\sigma_0$)"),
                     (a3, r"$B_2$ at $0.7<p_T/A<0.8$ GeV [GeV$^2$]", r"$B_2$ vs multiplicity")]:
        a.set_ylabel(yl); a.set_title(f"{t}: ALICE, {m}", fontsize=10); a.legend(fontsize=8)
        if row == 1: a.set_xlabel(r"$\langle dN_{ch}/d\eta\rangle_{|\eta|<0.5}$")
    a1.set_ylim(0, None); a3.set_ylim(0, None)
plt.tight_layout(); save(fig, "fig_mult_dp_B2")
NUM["mult"].update({short(m): {k: v for k, v in d.items() if k in ("inv", "chi2", "nch", "dp", "b2")} for m, d in MULT.items()})

# ------------------------------------------------------------------------------------------
# 1/sigma_0 across the experiments (cf. the paper's Tables VII, VIII)
# ------------------------------------------------------------------------------------------
PAPERV = {"Pythia 8.186 (paper)": {"LEP": 1.80, "BaBar": 1.13, "0.9 d": 3.58, "0.9 dbar": 3.63, "2.76 d": 2.93,
                                   "2.76 dbar": 2.88, "7 d": 2.63, "7 dbar": 2.58},
          "Herwig++ (paper)":     {"LEP": 1.93, "BaBar": 0.63, "0.9 d": 5.85, "0.9 dbar": 6.15, "2.76 d": 4.03,
                                   "2.76 dbar": 4.00, "7 d": 3.35, "7 dbar": 3.40}}
CATS = ["LEP", "BaBar", "0.9 d", "0.9 dbar", "2.76 d", "2.76 dbar", "7 d", "7 dbar", "7 mult."]
def inv_of(m, cat):
    if cat == "LEP": return LEP[m]["inv"], LEP[m]["err"]
    if cat == "BaBar": return BAB[m]["inv"], BAB[m]["err"]
    if cat == "7 mult.": return (MULT[m]["inv"], 0.) if m in MULT else (np.nan, 0.)
    rs, sp = cat.split(); f = FITS[(m, rs, sp)]; return f["inv"], f["err"]
fig, axes = plt.subplots(1, 2, figsize=(12, 4.3), sharey=True)
X = np.arange(len(CATS))
for ax, m in zip(axes, SETTINGS):
    for (name, vals), mk in zip(PAPERV.items(), ["v", "o"]):
        ax.plot([i for i, c in enumerate(CATS) if c in vals], [vals[c] for c in CATS if c in vals], marker=mk, ms=6, label=name, **PAPER)
    for k, mm in enumerate((PY, m)):
        v = [inv_of(mm, c) for c in CATS]
        ax.errorbar(X + 0.12*(2*k - 1), [a for a, b in v], yerr=[b for a, b in v], fmt="s", ms=5,
                    color=LOOK[mm]["color"], label=mm)
    ax.set_xticks(X); ax.set_xticklabels(["LEP", "BaBar", "0.9\n$d$", "0.9\n$\\bar d$", "2.76\n$d$", "2.76\n$\\bar d$",
                                          "7\n$d$", "7\n$\\bar d$", "7 TeV\nclasses"], fontsize=9)
    ax.minorticks_off(); ax.set_ylim(0, 7); ax.legend(fontsize=7.5, loc="upper right", ncol=2)
    ax.set_title(f"fitted $1/\\sigma_0$: Pythia, {m}", fontsize=10)
axes[0].set_ylabel(r"$1/\sigma_0$ [b$^{-1}$]")
plt.tight_layout(); save(fig, "fig_invsigma0")

json.dump(NUM, open(os.path.join(DEUTERON_TUNE_REPORT, "numbers.json"), "w"), indent=1, default=float)
print("wrote numbers.json")
