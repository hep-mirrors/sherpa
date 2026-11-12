"""Step 1 of the Lund retune: a few Lund string parameters, fitted to ALEPH at the Z pole with
Sherpa's own shower (the Monash 2013 values were fitted with Pythia's shower).

    python3 lund_scan.py make   <settings.json>    cards: the Monash point, a Latin hypercube, EXTRA_POINTS
    python3 lund_scan.py run    <settings.json>    run them, PARALLEL at a time (each ~1.5 min for 200k)
    python3 lund_scan.py score  <settings.json>    chi2 of every finished point against ALEPH, best first
    python3 lund_scan.py inputs <settings.json>    apprentice input: <APPRENTICE>/runs/<point>/{params.dat,
                                                   sherpa.yoda}, data.json, weights.txt

settings.json: {"name": ..., "base": "lund_tune1/card.yaml", "EVENTS": "200k", "PARALLEL": 10,
                "POINTS": 60, "SEED": 2026, "RANGES": {"StringZ:aLund": [lo, hi], ...},
                "EXTRA_POINTS": {"app_tune": {"StringZ:aLund": ..., ...}}}
Each point is the base card (lund_tune1: shower, MPI and colour reconnections unchanged) plus the
Pythia settings under PYTHIA8: PARAMETERS; Pythia does only the string fragmentation.
Runs go to SimpleLEP/LEPRivet/runs_<name>/<point>/, events through events.fifo to lep_observables.py.

Observables (ALEPH, as compare_lep.py): P(n_ch) (ALEPH_1991_I319520 d01), <n_ch> (d02), and the
x_p spectra of pi, K, p (ALEPH_1995_I382179 d01-d03).  Errors: data (+) 5 % of the value, as the
pp tune.  chi2 = sum over observables of W x chi2/n, with W = 1 for n_ch, <n_ch>, pi, K and 2 for p."""
import json, os, subprocess, sys, time
import numpy as np
import yaml

from paths import SCRIPTS, LEP_RIVET, APPRENTICE, SHERPA_LUND
SHERPA = SHERPA_LUND
REF = "/Users/ekremorhandemiray/code_software/rivetbootstrap/local/share/Rivet"
MODEL_UNCERTAINTY = 0.05
MONASH = {"StringZ:aLund": 0.68, "StringZ:bLund": 0.98, "StringFlav:probQQtoQ": 0.081, "StringZ:aExtraDiquark": 0.97}
# observable -> (analysis, histogram, weight); the pi/K/p entries are the species key in lep.json
OBSERVABLES = {"nch": ("ALEPH_1991_I319520", "d01-x01-y01", 1.0),
               "nch_mean": ("ALEPH_1991_I319520", "d02-x01-y01", 1.0),
               "pi": ("ALEPH_1995_I382179", "d01-x01-y01", 1.0),
               "K": ("ALEPH_1995_I382179", "d02-x01-y01", 1.0),
               "p": ("ALEPH_1995_I382179", "d03-x01-y01", 2.0)}

def load(settings_file):
    cfg = json.load(open(settings_file))
    cfg["runs"] = os.path.join(LEP_RIVET, f"runs_{cfg['name']}")
    return cfg

def points(cfg):
    "name -> {Pythia setting: value}: Monash, the Latin hypercube in RANGES, EXTRA_POINTS"
    names = list(cfg["RANGES"])
    rng = np.random.default_rng(cfg.get("SEED", 2026))
    n = cfg["POINTS"]
    result = {"ref_monash": {k: MONASH[k] for k in names}}
    columns = []
    for k in names:                                  # one stratum per point in every parameter
        lo, hi = cfg["RANGES"][k]
        columns.append(lo + (hi - lo)*(rng.permutation(n) + rng.random(n))/n)
    for i in range(n):
        result[f"l{i+1:03d}"] = {k: round(float(columns[j][i]), 4) for j, k in enumerate(names)}
    result.update(cfg.get("EXTRA_POINTS", {}))
    return result

def make(cfg):
    base = yaml.safe_load(open(os.path.join(LEP_RIVET, cfg["base"])))
    for name, values in points(cfg).items():
        d = os.path.join(cfg["runs"], name)
        os.makedirs(d, exist_ok=True)
        card = dict(base)
        card["EVENTS"] = cfg.get("EVENTS", "200k")
        card["PYTHIA8"] = dict(base.get("PYTHIA8", {}))
        card["PYTHIA8"]["PARAMETERS"] = [{k: v} for k, v in values.items()]
        with open(os.path.join(d, "card.yaml"), "w") as f:
            f.write(f"# lund_scan.py {cfg['name']} point {name}, from {cfg['base']}\n")
            for k, v in values.items(): f.write(f"#   {k} = {v}\n")
            yaml.safe_dump(card, f, sort_keys=False)
    print(len(points(cfg)), "cards in", cfg["runs"])

def done(d): return os.path.exists(os.path.join(d, "lep.json"))

def start(d):
    "the reader on the pipe, then Sherpa writing into it"
    fifo = os.path.join(d, "events.fifo")
    if os.path.exists(fifo): os.remove(fifo)
    os.mkfifo(fifo)
    reader = subprocess.Popen([sys.executable, os.path.join(SCRIPTS, "lep_observables.py"), "events.fifo", "lep.json"], cwd=d,
                              stdout=open(os.path.join(d, "reader.log"), "w"), stderr=subprocess.STDOUT)
    sherpa = subprocess.Popen(["caffeinate", "-i", SHERPA, "-f", "card.yaml"], cwd=d,
                              stdout=open(os.path.join(d, "run.log"), "w"), stderr=subprocess.STDOUT)
    return reader, sherpa

def run(cfg):
    todo = [os.path.join(cfg["runs"], name) for name in points(cfg)]
    todo = [d for d in todo if not done(d)]
    running = {}
    print(len(todo), "runs to do", flush=True)
    while todo or running:
        while todo and len(running) < cfg.get("PARALLEL", 10):
            d = todo.pop(0)
            running[d] = start(d)
            print(time.strftime("%H:%M"), "started", os.path.basename(d), flush=True)
        for d, (reader, sherpa) in list(running.items()):
            if sherpa.poll() is not None and reader.poll() is not None:
                print(time.strftime("%H:%M"), "finished", os.path.basename(d), "rc", sherpa.returncode, reader.returncode, flush=True)
                del running[d]
            elif sherpa.poll() is not None and sherpa.returncode != 0:   # Sherpa failed: the reader waits forever
                reader.kill(); print(time.strftime("%H:%M"), "FAILED", os.path.basename(d), flush=True); del running[d]
        time.sleep(5)

def reference():
    "ALEPH: name -> (bin lower edges, upper edges, values, errors incl. model uncertainty)"
    import yoda
    ref = {}
    for name, (analysis, hist, _) in OBSERVABLES.items():
        ao = yoda.read(os.path.join(REF, analysis + ".yoda.gz"))[f"/REF/{analysis}/{hist}"]
        val = np.array([b.val() for b in ao.bins()]); err = np.array([b.totalErrAvg() for b in ao.bins()])
        if name == "nch":
            n = np.array(ao.xEdges(), dtype=float); lo, hi = n - 1, n + 1      # discrete: even n_ch
        elif name == "nch_mean":
            lo, hi = np.array([0.]), np.array([1.])
        else:
            lo = np.array([b.xMin() for b in ao.bins()]); hi = np.array([b.xMax() for b in ao.bins()])
        ref[name] = (lo, hi, val, np.hypot(err, MODEL_UNCERTAINTY*val))
    return ref

def predictions(d, ref):
    "a finished point on ALEPH's bins: name -> (values, statistical errors)"
    r = json.load(open(os.path.join(d, "lep.json")))
    N = r["events"]
    counts = {int(k): v for k, v in r["nch"].items()}
    out = {}
    n = 0.5*(ref["nch"][0] + ref["nch"][1])
    c = np.array([counts.get(int(round(x)), 0) for x in n], dtype=float)
    out["nch"] = (c/N, np.sqrt(np.maximum(c, 1))/N)
    mean = sum(k*v for k, v in counts.items())/N
    second = sum(k*k*v for k, v in counts.items())/N
    out["nch_mean"] = (np.array([mean]), np.array([np.sqrt((second - mean**2)/N)]))
    for sp in ("pi", "K", "p"):
        lo, hi = ref[sp][0], ref[sp][1]
        xp = np.array(r["xp"][sp])
        c = np.array([((xp >= l) & (xp < h)).sum() for l, h in zip(lo, hi)], dtype=float)
        out[sp] = (c/(N*(hi - lo)), np.sqrt(np.maximum(c, 1))/(N*(hi - lo)))
    return out

def chi2(pred, ref):
    "name -> chi2/n, and the weighted total"
    parts = {}
    for name, (_, _, val, err) in ref.items():
        mc, mce = pred[name]
        ok = err > 0
        parts[name] = float((((mc - val)**2/(err**2 + mce**2))[ok]).mean())
    return parts, sum(OBSERVABLES[k][2]*v for k, v in parts.items())

def score(cfg):
    ref = reference()
    rows = []
    for name, values in points(cfg).items():
        d = os.path.join(cfg["runs"], name)
        if not done(d): continue
        pred = predictions(d, ref)
        parts, total = chi2(pred, ref)
        rows.append((total, name, values, parts, pred["nch_mean"][0][0]))
    names = list(cfg["RANGES"])
    print(f"{'point':12s}" + "".join(f"{k.split(':')[1]:>14s}" for k in names)
          + f"{'<n_ch>':>8s}" + "".join(f"{k:>9s}" for k in OBSERVABLES) + f"{'TOTAL':>9s}")
    print(f"{'ALEPH':12s}" + " "*14*len(names) + f"{ref['nch_mean'][2][0]:8.2f}")
    for total, name, values, parts, mean in sorted(rows, key=lambda row: row[0]):
        print(f"{name:12s}" + "".join(f"{values[k]:14.4f}" for k in names) + f"{mean:8.2f}"
              + "".join(f"{parts[k]:9.2f}" for k in OBSERVABLES) + f"{total:9.2f}")

def inputs(cfg):
    "apprentice input in <APPRENTICE> (default SimpleLEP/TuneScan/apprentice/<name>)"
    import yoda
    out = cfg.get("APPRENTICE", os.path.join(APPRENTICE, cfg["name"]))
    ref = reference()
    used = 0
    for name, values in points(cfg).items():
        d = os.path.join(cfg["runs"], name)
        if not done(d): continue
        target = os.path.join(out, "runs", name)
        os.makedirs(target, exist_ok=True)
        with open(os.path.join(target, "params.dat"), "w") as f:
            for k, v in values.items(): f.write(f"{k.split(':')[1]} {v}\n")
        objects = []
        for obs, (mc, mce) in predictions(d, ref).items():
            edges = list(np.arange(len(mc) + 1, dtype=float))    # bin index: the edges do not matter to apprentice
            est = yoda.Estimate1D(edges, f"/ALEPH/{obs}")
            for i, (v, e) in enumerate(zip(mc, mce)):
                est.bin(i+1).setVal(float(v)); est.bin(i+1).setErr((-float(e), float(e)))
            objects.append(est)
        yoda.write(objects, os.path.join(target, "sherpa.yoda"))
        used += 1
    data = {}
    with open(os.path.join(out, "weights.txt"), "w") as f:
        for obs, (_, _, val, err) in ref.items():
            for i, (v, e) in enumerate(zip(val, err)):
                data[f"/ALEPH/{obs}#{i}"] = [float(v), float(e)]
            f.write(f"/ALEPH/{obs} {np.sqrt(OBSERVABLES[obs][2]/len(val)):.6f}\n")
    json.dump(data, open(os.path.join(out, "data.json"), "w"), indent=1)
    print(f"{used} points -> {out}/runs, data.json, weights.txt")

if __name__ == "__main__":
    if len(sys.argv) != 3: sys.exit(__doc__)
    {"make": make, "run": run, "score": score, "inputs": inputs}[sys.argv[1]](load(sys.argv[2]))
