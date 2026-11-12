"""Runs every Lund_<setting> card (make_lund_cards.py) with the Pythia-8-enabled Sherpa, as the
ScanTune runs of that setting were made, then merges the 7 TeV seeds.

    python3 run_lund.py s6_008 [parallel jobs, default 8]

The cards are in <report>/runcards/Lund_<setting> (paths.card_set).
Each Sherpa job runs in its own directory SimpleLEP/LundRuns/<setting>/<card>/ (Sherpa writes Results/ and
Process/ into its working directory, which parallel jobs must not share); the OUTPUT of the
deuteron phase is made absolute and stays DeuteronFormationFromXSecs/lund_<setting>/<setup>.
ALICE 7 TeV, as for ScanTune s6_008: the card's own run plus seeds 7001-7003 (100k events each),
merged with merge_seeds.py into lund_<setting>/ALICE_7TeV (400k events).
Finished runs (summary.dat present) are not run again."""
import glob, os, subprocess, sys, time

from paths import SCRIPTS, SIMPLE_LEP, DEUTERON_DATA, LUND_RUNS, SHERPA_LUND, card_set
SIMPLELEP = SIMPLE_LEP          # the card OUTPUT paths are relative to it
SHERPA = SHERPA_LUND
EXTRA_7TEV_SEEDS = [7001, 7002, 7003]

def jobs(setting):
    "(run directory, card text, OUTPUT directory) of every run"
    out = []
    for card in sorted(glob.glob(os.path.join(card_set(f"Lund_{setting}"), "*.Lund.yaml"))):
        name = os.path.basename(card).replace(".Lund.yaml", "")
        text = open(card).read()
        output = [l.split("OUTPUT:")[1].split("#")[0].strip() for l in text.splitlines() if "OUTPUT:" in l][0]
        variants = [(name, text, output)]
        if name.endswith("ALICE_7TeV"):
            seeds_dir = output.replace("/ALICE_7TeV", "/ALICE_7TeV_seeds")
            variants = [(name, text, os.path.join(seeds_dir, "card_seed"))]
            for seed in EXTRA_7TEV_SEEDS:
                variants.append((f"{name}.seed{seed}", text.replace("EVENTS:", f"RANDOM_SEED:    {seed}\nEVENTS:", 1),
                                 os.path.join(seeds_dir, f"seed{seed}")))
        for run, card_text, outdir in variants:
            absolute = os.path.join(SIMPLELEP, outdir)
            card_text = "\n".join(l.replace(output, absolute) if "OUTPUT:" in l else l for l in card_text.splitlines()) + "\n"
            out.append((os.path.join(LUND_RUNS, setting, run), card_text, absolute))
    return out

def main(setting, parallel):
    todo = [j for j in jobs(setting) if not os.path.exists(os.path.join(j[2], "summary.dat"))]
    print(f"{len(todo)} runs to do, {parallel} at a time")
    running = {}
    while todo or running:
        while todo and len(running) < parallel:
            rundir, text, outdir = todo.pop(0)
            os.makedirs(rundir, exist_ok=True)
            with open(os.path.join(rundir, "card.yaml"), "w") as f: f.write(text)
            log = open(os.path.join(rundir, "run.log"), "w")
            running[rundir] = subprocess.Popen(["caffeinate", "-i", SHERPA, "-f", "card.yaml"], cwd=rundir,
                                               stdout=log, stderr=subprocess.STDOUT)
            print(time.strftime("%H:%M"), "started", os.path.basename(rundir), flush=True)
        for rundir, proc in list(running.items()):
            if proc.poll() is not None:
                print(time.strftime("%H:%M"), "finished", os.path.basename(rundir), "rc", proc.returncode, flush=True)
                del running[rundir]
        time.sleep(10)
    seeds = sorted(glob.glob(os.path.join(DEUTERON_DATA, f"lund_{setting}", "ALICE_7TeV_seeds", "*")))
    if seeds:
        target = os.path.join(DEUTERON_DATA, f"lund_{setting}", "ALICE_7TeV")
        subprocess.run([sys.executable, os.path.join(SCRIPTS, "merge_seeds.py"), target] + seeds, check=True)
        print("merged", len(seeds), "7 TeV runs into", target)

if __name__ == "__main__":
    if len(sys.argv) not in (2, 3): sys.exit(__doc__)
    main(sys.argv[1], int(sys.argv[2]) if len(sys.argv) == 3 else 8)
