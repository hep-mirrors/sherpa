"""Energy dependence of the multiple scattering for Sherpa with Lund strings: scan of AMISIC's Eta.

    python3 eta_scan.py make  <settings.json>    cards for every Eta in ETAS at 0.9 and 2.76 TeV
    python3 eta_scan.py run   <settings.json>    run them (PARALLEL at a time)
    python3 eta_scan.py score <settings.json>    <dNch/deta> against ALICE -> table
    python3 eta_scan.py plot  <settings.json> <figure stem>   dNch ratios and fitted 1/sigma_0 against Eta

In AMISIC (MI_Parameters.C:141-161) the regularisation scale and the lower cut of the
scatterings scale with energy as pT0^2(s) = PT_0(ref)^2 (s/s_ref)^(2 Eta), the same for PT_Min,
with s_ref = (7 TeV)^2.  Eta therefore leaves 7 TeV unchanged and only moves the other energies.
A smaller Eta means a larger pT0 below 7 TeV, i.e. fewer scatterings there.

settings.json: {"name": ..., "base": "Lund_tune1", "AMISIC": {...}, "COLOUR_RECONNECTIONS": {...},
                "ETAS": [...], "EVENTS": "100k", "PARALLEL": 10}
Any other card parameter can be scanned instead of Eta, with
                "PARAM": ["REMNANTS", "2212", "ENERGY_SCALING_EXPO"], "TAG": "kexpo", "VALUES": [...],
                "LABEL": "REMNANTS ENERGY_SCALING_EXPO", "DEFAULT": 0.08
"base" is a set of run cards (paths.card_set).  The settings files are in SimpleLEP/TuneScan.
Runs go to SimpleLEP/TuneScan/runs_<TAG>_<name>/<energy>_<TAG><value>/ (TAG "eta" by default).  Score: <dNch/deta> at |eta|<0.5 in
INEL>0 events against ALICE (HEPData ins1394854, Table 15), and its ratio to the 7 TeV value of
the same setting (which removes the common offset of non-diffractive against INEL>0 events)."""
import json, os, subprocess, sys, time
import numpy as np
import yaml

from paths import TUNE_SCAN, SHERPA_LUND, card_set
SHERPA = SHERPA_LUND
ALICE_NCH = {"0.9": (3.75, 0.06), "2.76": (4.76, 0.08), "7": (5.98, 0.09)}   # ins1394854 Table 15, INEL>0
ENERGIES = ["0.9", "2.76"]

def load(settings_file):
    cfg = json.load(open(settings_file))
    cfg["param"] = cfg.get("PARAM", ["AMISIC", "Eta"])
    cfg["tag"] = cfg.get("TAG", "eta")
    cfg["values"] = cfg.get("VALUES", cfg.get("ETAS"))
    cfg["runs"] = os.path.join(TUNE_SCAN, f"runs_{cfg['tag']}_{cfg['name']}")
    return cfg

def points(cfg):
    return [(rs, v, point_dir(cfg, rs, v)) for v in cfg["values"] for rs in ENERGIES]

def point_dir(cfg, rs, value): return os.path.join(cfg["runs"], f"{rs}TeV_{cfg['tag']}{value:+.3f}")

def make(cfg):
    for rs, value, d in points(cfg):
        card = yaml.safe_load(open(os.path.join(card_set(cfg["base"]), f"Sherpa.Deuterons.ALICE_{rs}TeV.Lund.yaml")))
        for block in ("AMISIC", "COLOUR_RECONNECTIONS"):
            card.setdefault(block, {}).update(cfg.get(block, {}))
        block = card
        for key in cfg["param"][:-1]:          # YAML reads REMNANTS: 2212: as an integer key
            key = next((k for k in block if str(k) == key), key)
            block = block.setdefault(key, {})
        block[cfg["param"][-1]] = value
        card["EVENTS"] = cfg.get("EVENTS", "100k")
        card["DEUTERON_FORMATION"]["OUTPUT"] = os.path.join(d, "out")
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, "card.yaml"), "w") as f:
            f.write(f"# eta_scan.py {cfg['name']}: {rs} TeV, {':'.join(cfg['param'])} = {value}\n")
            yaml.safe_dump(card, f, sort_keys=False)
    print(len(points(cfg)), "cards in", cfg["runs"])

def done(d): return os.path.exists(os.path.join(d, "out", "summary.dat"))

def run(cfg):
    todo = [d for _, _, d in points(cfg) if not done(d)]
    running = {}
    print(len(todo), "runs to do")
    while todo or running:
        while todo and len(running) < cfg.get("PARALLEL", 10):
            d = todo.pop(0)
            running[d] = subprocess.Popen(["caffeinate", "-i", SHERPA, "-f", "card.yaml"], cwd=d,
                                          stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT)
            print(time.strftime("%H:%M"), "started", os.path.basename(d), flush=True)
        for d, p in list(running.items()):
            if p.poll() is not None:
                print(time.strftime("%H:%M"), "finished", os.path.basename(d), "rc", p.returncode, flush=True); del running[d]
        time.sleep(10)

def nch(outdir):
    m = np.atleast_2d(np.loadtxt(os.path.join(outdir, "multiplicity_events.dat")))
    return m[:, 2].sum()/m[:, 1].sum()

def score(cfg):
    ref7 = cfg.get("NCH_7TEV")          # the setting's own <dNch/deta> at 7 TeV (from its 7 TeV run)
    print(f"{cfg['tag']:>7s} " + "".join(f"{rs + ' TeV: MC/ALICE':>20s}{'ratio to 7 TeV':>16s}" for rs in ENERGIES))
    print(f"{'ALICE':>7s} " + "".join(f"{'1':>20s}{ALICE_NCH[rs][0]/ALICE_NCH['7'][0]:16.3f}" for rs in ENERGIES)
          + (f"     (7 TeV MC/ALICE {ref7/ALICE_NCH['7'][0]:.3f})" if ref7 else ""))
    for value in cfg["values"]:
        row = f"{value:7.3f} "
        for rs in ENERGIES:
            d = os.path.join(point_dir(cfg, rs, value), "out")
            if not os.path.exists(os.path.join(d, "multiplicity_events.dat")): row += f"{'-':>20s}{'-':>16s}"; continue
            v = nch(d)
            row += f"{v/ALICE_NCH[rs][0]:20.3f}{(v/ref7 if ref7 else float('nan')):16.3f}"
        print(row)

def plot(cfg, stem):
    "left: <dNch/deta> relative to 7 TeV against Eta; right: the deuteron 1/sigma_0 fitted at each energy"
    import contextlib, io
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    with contextlib.redirect_stdout(io.StringIO()):
        import scorecard
    ref7, inv7 = cfg["NCH_7TEV"], cfg["INV_SIGMA0_7TEV"]
    sqrts = {"0.9": 900., "2.76": 2760.}; table = {"0.9": 1, "2.76": 2}
    colour = {"0.9": "#0072B2", "2.76": "#D55E00"}
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.2))
    for rs in ENERGIES:
        etas, ratio, inv_d, inv_dbar = [], [], [], []
        for eta in cfg["values"]:
            d = os.path.join(point_dir(cfg, rs, eta), "out")
            if not os.path.exists(os.path.join(d, "summary.dat")): continue
            fit = scorecard.deuterons(d, sqrts[rs], table[rs])
            etas.append(eta); ratio.append(nch(d)/ref7); inv_d.append(fit["d"][0]); inv_dbar.append(fit["dbar"][0])
        ax.plot(etas, ratio, "o-", color=colour[rs], label=f"{rs} TeV, MC")
        ax.axhline(ALICE_NCH[rs][0]/ALICE_NCH["7"][0], color=colour[rs], ls="--", label=f"{rs} TeV, ALICE")
        bx.plot(etas, inv_d, "o-", color=colour[rs], label=f"{rs} TeV, $d$")
        bx.plot(etas, inv_dbar, "s:", color=colour[rs], label=rf"{rs} TeV, $\bar d$")
    bx.axhline(inv7[0], color="k", ls="-", lw=0.8, label=f"7 TeV, $d$ ({inv7[0]:.2f})")
    bx.axhline(inv7[1], color="k", ls=":", lw=0.8, label=rf"7 TeV, $\bar d$ ({inv7[1]:.2f})")
    for a in (ax, bx):
        a.axvline(cfg.get("DEFAULT", 0.026), color="0.6", lw=0.8); a.set_xlabel(cfg.get("LABEL", "AMISIC Eta")); a.legend(fontsize=7)
    ax.set_ylabel(r"$\langle dN_{\rm ch}/d\eta\rangle$ / its 7 TeV value"); ax.set_title("Charged-particle density")
    bx.set_ylabel(r"fitted $1/\sigma_0$ [b$^{-1}$]"); bx.set_title(r"Deuterons: $1/\sigma_0$ needed at each energy")
    bx.set_yscale("log")
    plt.tight_layout()
    for ext in ("png", "pdf"): fig.savefig(f"{stem}.{ext}", dpi=130)
    print("wrote", stem + ".png/.pdf")

if __name__ == "__main__":
    if len(sys.argv) < 3: sys.exit(__doc__)
    if sys.argv[1] == "plot": plot(load(sys.argv[2]), sys.argv[3])
    else: {"make": make, "run": run, "score": score}[sys.argv[1]](load(sys.argv[2]))
