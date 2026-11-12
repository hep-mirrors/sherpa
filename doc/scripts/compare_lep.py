"""Compares the LEP runs (lep_observables.py output) with ALEPH's measurements at the Z pole,
read from Rivet's reference files: charged multiplicity (ALEPH_1991_I319520) and pi, K, p
momentum spectra (ALEPH_1995_I382179).  Normalised as the Rivet analyses: per event, per unit x_p.

    python3 compare_lep.py <figure dir> setup1 setup2 ...     (setups = directories with lep.json)

Writes lep_nch.{png,pdf}, lep_xp.{png,pdf} and prints <n_ch>, chi2 per point and MC/data."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, "/Users/ekremorhandemiray/code_software/rivetbootstrap/local/lib/python3.12/site-packages")
import yoda

REF = "/Users/ekremorhandemiray/code_software/rivetbootstrap/local/share/Rivet"
LABEL = {"lund_tune1": "Sherpa Lund tuned", "s6_008": "Sherpa s6_008 (AHADIC)", "sherpa_default": "Sherpa default (AHADIC)"}
STYLE = {"lund_tune1": dict(color="#E69F00", ls="-", lw=2.2), "s6_008": dict(color="#D55E00", ls=":", lw=2.2),
         "sherpa_default": dict(color="0.45", ls="--", lw=1.6)}

# final report (doc/final_report): the final setting, the LEP-tuned strings (Lund_tune4), Pythia
LABEL["lund_mbtune"] = "Sherpa final setting"
STYLE["lund_mbtune"] = dict(color="#CC0000", ls="-", lw=2.2)
LABEL["app_tune"] = "Sherpa, LEP-tuned strings (Lund_tune4)"
STYLE["app_tune"] = dict(color="#009E73", ls="--", lw=1.6)
LABEL["pythia"] = "Pythia 8.318 (Monash)"
STYLE["pythia"] = dict(color="#0072B2", ls="-.", lw=1.8)

FS_COLOURS = {"0.25": "#0072B2", "0.4": "#009E73", "0.6": "#CC79A7"}
for v, c in FS_COLOURS.items():
    LABEL[f"lund_fs{v}"] = f"Sherpa Lund, FS_PT2MIN {v}"
    STYLE[f"lund_fs{v}"] = dict(color=c, ls="-", lw=1.6)

AS_COLOURS = {"0.8": "#56B4E9", "0.66": "#009E73", "0.5": "#0072B2", "0.35": "#CC79A7"}
for v, c in AS_COLOURS.items():
    LABEL[f"lund_as{v}"] = f"Sherpa Lund, FS_AS_FAC {v}"
    STYLE[f"lund_as{v}"] = dict(color=c, ls="-", lw=1.6)

def ref(analysis, name):
    ao = yoda.read(os.path.join(REF, analysis + ".yoda.gz"))[f"/REF/{analysis}/{name}"]
    return ao

def main(figdir, setups):
    os.makedirs(figdir, exist_ok=True)
    runs = {s: json.load(open(os.path.join(s, "lep.json"))) for s in setups}
    # ---- charged multiplicity
    d01 = ref("ALEPH_1991_I319520", "d01-x01-y01")
    n_data = np.array(d01.xEdges(), dtype=float)
    p_data = np.array([b.val() for b in d01.bins()]); e_data = np.array([b.totalErrAvg() for b in d01.bins()])
    mean = ref("ALEPH_1991_I319520", "d02-x01-y01").bins()[0]
    print(f"<n_ch>: ALEPH {mean.val():.2f} +- {mean.totalErrAvg():.2f}")
    fig, (ax, rx) = plt.subplots(2, 1, figsize=(7, 6), sharex=True, gridspec_kw={"height_ratios": [3, 1]})
    ax.errorbar(n_data, p_data, yerr=e_data, fmt="o", color="black", ms=4, label="ALEPH")
    for s, r in runs.items():
        N = r["events"]; counts = {int(k): v for k, v in r["nch"].items()}
        mc = np.array([counts.get(int(n), 0)/N for n in n_data])
        avg = sum(n*c for n, c in counts.items())/N
        ok = e_data > 0
        chi2 = (((mc - p_data)/e_data)[ok]**2).mean()
        print(f"  {LABEL[s]:28s} <n_ch> {avg:.2f}   P(n_ch) chi2/point {chi2:.2f}   ({N} events)")
        ax.plot(n_data, mc, label=f"{LABEL[s]}, $\\langle n\\rangle$ = {avg:.1f}", **STYLE[s])
        rx.plot(n_data[p_data > 0], mc[p_data > 0]/p_data[p_data > 0], **STYLE[s])
    ax.set_yscale("log"); ax.set_ylabel(r"$P(n_{\rm ch})$"); ax.legend(fontsize=8)
    ax.set_title(r"Charged multiplicity in $Z\to$ hadrons (ALEPH)")
    rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0, 2); rx.set_xlabel(r"$n_{\rm ch}$"); rx.set_ylabel("MC / ALEPH")
    plt.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(figdir, f"lep_nch.{ext}"), dpi=130)
    plt.close(fig)
    # ---- pi, K, p spectra
    fig, axes = plt.subplots(2, 3, figsize=(16, 6), sharex="col", gridspec_kw={"height_ratios": [3, 1]})
    for col, (sp, name, title) in enumerate((("pi", "d01-x01-y01", r"$\pi^\pm$"), ("K", "d02-x01-y01", r"$K^\pm$"),
                                             ("p", "d03-x01-y01", r"$p,\bar p$"))):
        ao = ref("ALEPH_1995_I382179", name)
        lo = np.array([b.xMin() for b in ao.bins()]); hi = np.array([b.xMax() for b in ao.bins()])
        val = np.array([b.val() for b in ao.bins()]); err = np.array([b.totalErrAvg() for b in ao.bins()])
        x = 0.5*(lo + hi)
        ax, rx = axes[0, col], axes[1, col]
        ax.errorbar(x, val, yerr=err, xerr=0.5*(hi - lo), fmt="o", color="black", ms=3, label="ALEPH")
        rx.fill_between(x, 1 - err/val, 1 + err/val, color="0.85", step="mid")
        print(f"{title} x_p spectrum:")
        for s, r in runs.items():
            xp = np.array(r["xp"][sp]); N = r["events"]
            mc = np.array([((xp >= l) & (xp < h)).sum()/(N*(h - l)) for l, h in zip(lo, hi)])
            chi2 = (((mc - val)/err)**2).mean()
            ratio = mc/val
            def rng(a, b): m = (x >= a) & (x < b); return ratio[m].mean() if m.any() else float("nan")
            print(f"  {LABEL[s]:28s} per event {len(xp)/N:.3f}   chi2/point {chi2:6.1f}   MC/ALEPH at x_p<0.05 {rng(0, .05):.2f}, "
                  f"0.05-0.2 {rng(.05, .2):.2f}, 0.2-0.5 {rng(.2, .5):.2f}, >0.5 {rng(.5, 1.01):.2f}")
            ax.plot(x, mc, label=LABEL[s], **STYLE[s]); rx.plot(x, ratio, **STYLE[s])
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_title(f"{title} in $Z\\to$ hadrons (ALEPH)"); ax.legend(fontsize=7)
        rx.axhline(1, color="k", lw=0.8); rx.set_ylim(0.5, 1.5); rx.set_xlabel(r"$x_p=|p|/E_{\rm beam}$")
    axes[0, 0].set_ylabel(r"$\frac{1}{N}\frac{dN}{dx_p}$"); axes[1, 0].set_ylabel("MC / ALEPH")
    plt.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(os.path.join(figdir, f"lep_xp.{ext}"), dpi=130)
    plt.close(fig)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2:])
