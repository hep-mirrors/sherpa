"""Coalescence runs for the lund_tune report (tunedLund.tex): Sherpa with Lund strings (Lund_tune1,
Lund_tune4) and Pythia 8.318, with the coalescence models of Deuteron_Production.C instead of the
cross-section model.

    python3 run_coalescence.py list               the runs and their output directories
    python3 run_coalescence.py run [parallel]     run every run that has not finished (summary.dat)
    python3 run_coalescence.py seeds [parallel]   the same, plus the extra seeds of the Sherpa pp runs
    python3 run_coalescence.py merge              merge each pp run with its seeds -> <setup>_merged

Models (DEUTERON_FORMATION: MODEL), as the earlier coalescence study (deuteron_report/deuteron_report.tex):
  HardDisk     P = 1 if k < p_coal (and b < b_coal)
  Exponential  P = 0.75 exp(-b^2/C_b) exp(-k^2/C_p), C = coal^2/(0.75 pi)
  _k   momentum condition only (B_COAL: 0), as Dal & Raklev's coalescence (their Eq. 1)
  _kb  with the impact parameter b < 3 fm; Sherpa needs PYTHIA8: SET_VERTICES: true for the hadrons
       to have production points (Pythia's string-break vertices), Pythia sets Fragmentation:setVertices
p_coal = 0.2 GeV, b_coal = 3 fm.  Scans at LEP and 7 TeV (HardDisk): p_coal 0.15, 0.3 (k only) and
b_coal 2, 5 fm.
Output (in SimpleLEP): DeuteronFormationFromXSecs/coalescence_lund/<tune>/<model>/<setup> and
        DeuteronFormationFromXSecs/coalescence_lund/pythia/<model>/<setup>;
work directories (Sherpa's Process/, Results/, logs): LundRuns/coalescence/..."""
import glob, os, subprocess, sys, time

from paths import SCRIPTS, DEUTERON_DATA, LUND_RUNS, PYTHIA_REFERENCE, SHERPA_LUND, card_set
SHERPA = SHERPA_LUND
PYTHIA = os.path.join(PYTHIA_REFERENCE, "pythia_deuterons")
OUT = os.path.join(DEUTERON_DATA, "coalescence_lund")
WORK = os.path.join(LUND_RUNS, "coalescence")
TUNES = ["tune1", "tune4"]
# setup -> (card, Pythia sqrt(s), Pythia events, Pythia beams, seed)
SETUPS = {"LEP":            ("LEP",           91.2,   1000000, "ee",   91),
          "BaBar_10.58GeV": ("BaBar",         10.58,  2000000, "ee",   10),
          "ALICE_0.9TeV":   ("ALICE_0.9TeV",  900.,   500000,  "pp",  911),
          "ALICE_2.76TeV":  ("ALICE_2.76TeV", 2760.,  500000,  "pp", 2771),
          "ALICE_7TeV":     ("ALICE_7TeV",    7000.,  500000,  "pp", 7011)}
P_COAL, B_COAL = 0.2, 3.

def variants():
    "name -> (model, p_coal, b_coal, setups)"
    out = {}
    for model in ("HardDisk", "Exponential"):
        out[f"{model}_k"] = (model, P_COAL, 0., list(SETUPS))
        out[f"{model}_kb"] = (model, P_COAL, B_COAL, list(SETUPS))
    for p in (0.15, 0.3):
        out[f"HardDisk_k_p{p}"] = ("HardDisk", p, 0., ["LEP", "ALICE_7TeV"])
    for b in (2., 5.):
        out[f"HardDisk_kb_b{b:g}"] = ("HardDisk", P_COAL, b, ["LEP", "ALICE_7TeV"])
    return out

def sherpa_card(tune, setup, model, p, b, outdir):
    card = os.path.join(card_set(f"Lund_{tune}"), f"Sherpa.Deuterons.{SETUPS[setup][0]}.Lund.yaml")
    lines = []
    for line in open(card).read().splitlines():
        stripped = line.strip()
        if stripped.startswith("MODEL:") and line.startswith("  "):
            lines += [f"  MODEL:      {model}", f"  P_COAL:     {p}", f"  B_COAL:     {b}"]
            continue
        if stripped.startswith("OUTPUT:") and line.startswith("  "):
            lines.append(f"  OUTPUT:     {outdir}"); continue
        if stripped.startswith("TEST:") and line.startswith("  "):
            lines.append("  TEST:       false"); continue
        lines.append(line)
        if line == "  DECAYS: false" and b > 0.:
            lines.append("  SET_VERTICES: true   # hadron production points from Pythia (needed for b)")
    return f"# run_coalescence.py: Lund_{tune}, {model}, p_coal {p} GeV, b_coal {b} fm\n" + "\n".join(lines) + "\n"

def jobs():
    "(name, work directory, command, output directory)"
    out = []
    for name, (model, p, b, setups) in variants().items():
        for setup in setups:
            for tune in TUNES:
                outdir = os.path.join(OUT, tune, name, setup)
                work = os.path.join(WORK, tune, name, setup)
                out.append((f"{tune}/{name}/{setup}", work, ("sherpa", sherpa_card(tune, setup, model, p, b, outdir)), outdir))
            _, sqrts, nev, beams, seed = SETUPS[setup]
            outdir = os.path.join(OUT, "pythia", name, setup)
            work = os.path.join(WORK, "pythia", name, setup)
            cmd = [PYTHIA, str(sqrts), str(nev), str(seed), "1.8", outdir, beams, f"{model}:{p}:{b}"]
            out.append((f"pythia/{name}/{setup}", work, ("pythia", cmd), outdir))
    return out

# extra seeds for the Sherpa pp runs: coalescence uses only the few pairs at small k, so 100k
# events leave the spectra dominated by statistical fluctuations.  Seeds 1..n go to
# <setup>_seeds/seed_<n>; "merge" averages them with the first run into <setup>_merged.
EXTRA_SEEDS = {"HardDisk_k": {"ALICE_0.9TeV": 4, "ALICE_2.76TeV": 4, "ALICE_7TeV": 9},
               "Exponential_k": {"ALICE_0.9TeV": 4, "ALICE_2.76TeV": 4, "ALICE_7TeV": 9},
               "HardDisk_kb": {"ALICE_0.9TeV": 4, "ALICE_2.76TeV": 4, "ALICE_7TeV": 4},
               "Exponential_kb": {"ALICE_0.9TeV": 4, "ALICE_2.76TeV": 4, "ALICE_7TeV": 4}}

# Pythia is fast (500k pp events in ~12 min): momentum-only models to 5M (7 TeV) and 2.5M events
PYTHIA_EXTRA_SEEDS = {"HardDisk_k": {"ALICE_0.9TeV": 4, "ALICE_2.76TeV": 4, "ALICE_7TeV": 9},
                      "Exponential_k": {"ALICE_0.9TeV": 4, "ALICE_2.76TeV": 4, "ALICE_7TeV": 9}}

def pythia_seed_jobs():
    out = []
    for name, setups in PYTHIA_EXTRA_SEEDS.items():
        model, p, b, _ = variants()[name]
        for setup, n in setups.items():
            _, sqrts, nev, beams, seed0 = SETUPS[setup]
            for i in range(1, n + 1):
                seed = seed0 + 1000*i
                outdir = os.path.join(OUT, "pythia", name, setup + "_seeds", f"seed_{seed}")
                work = os.path.join(WORK, "pythia", name, setup + f"_seed{seed}")
                cmd = [PYTHIA, str(sqrts), str(nev), str(seed), "1.8", outdir, beams, f"{model}:{p}:{b}"]
                out.append((f"pythia/{name}/{setup}/seed{seed}", work, ("pythia", cmd), outdir))
    return out

def seed_jobs():
    out = []
    for name, setups in EXTRA_SEEDS.items():
        model, p, b, _ = variants()[name]
        for setup, n in setups.items():
            for tune in TUNES:
                for i in range(1, n + 1):
                    seed = 1000*i + 7
                    outdir = os.path.join(OUT, tune, name, setup + "_seeds", f"seed_{seed}")
                    work = os.path.join(WORK, tune, name, setup + f"_seed{seed}")
                    card = sherpa_card(tune, setup, model, p, b, outdir).replace("EVENTS:", f"RANDOM_SEED:    {seed}\nEVENTS:", 1)
                    out.append((f"{tune}/{name}/{setup}/seed{seed}", work, ("sherpa", card), outdir))
    return out

def merge():
    for name, setups in EXTRA_SEEDS.items():
        for setup in setups:
            for tune in TUNES + ["pythia"]:
                first = os.path.join(OUT, tune, name, setup)
                seeds = sorted(glob.glob(os.path.join(OUT, tune, name, setup + "_seeds", "seed_*")))
                runs = [first] + [d for d in seeds if done(d)]
                if not done(first) or len(runs) < 2: continue
                subprocess.run([sys.executable, os.path.join(SCRIPTS, "merge_seeds.py"),
                                os.path.join(OUT, tune, name, setup + "_merged")] + runs, check=True)

def done(outdir): return os.path.exists(os.path.join(outdir, "summary.dat"))

def start(work, what):
    os.makedirs(work, exist_ok=True)
    log = open(os.path.join(work, "run.log"), "w")
    if what[0] == "sherpa":
        with open(os.path.join(work, "card.yaml"), "w") as f: f.write(what[1])
        return subprocess.Popen(["caffeinate", "-i", SHERPA, "-f", "card.yaml"], cwd=work, stdout=log, stderr=subprocess.STDOUT)
    return subprocess.Popen(["caffeinate", "-i"] + what[1], cwd=work, stdout=log, stderr=subprocess.STDOUT)

def run(parallel, with_seeds=False):
    # the long pp runs first
    todo = sorted([j for j in jobs() + (seed_jobs() if with_seeds else []) if not done(j[3])],
                  key=lambda j: ("ALICE" not in j[0], "7TeV" not in j[0], j[0]))
    running = {}
    print(len(todo), "runs to do,", parallel, "at a time", flush=True)
    while todo or running:
        while todo and len(running) < parallel:
            name, work, what, outdir = todo.pop(0)
            os.makedirs(outdir, exist_ok=True)       # pythia_deuterons makes only the last directory
            running[name] = start(work, what)
            print(time.strftime("%H:%M"), "started", name, flush=True)
        for name, p in list(running.items()):
            if p.poll() is not None:
                print(time.strftime("%H:%M"), "finished", name, "rc", p.returncode, flush=True); del running[name]
        time.sleep(10)

if __name__ == "__main__":
    if len(sys.argv) < 2: sys.exit(__doc__)
    if sys.argv[1] == "list":
        for name, work, what, outdir in jobs(): print(("done  " if done(outdir) else "todo  ") + name)
    elif sys.argv[1] == "run": run(int(sys.argv[2]) if len(sys.argv) > 2 else 14)
    elif sys.argv[1] == "seeds": run(int(sys.argv[2]) if len(sys.argv) > 2 else 14, with_seeds=True)
    elif sys.argv[1] == "merge": merge()
