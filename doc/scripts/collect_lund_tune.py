"""Copies everything of the final setting (Lund_mbtune) into doc/lund_tune/, one folder per
part, with the run cards next to the output they made and the Pythia counterpart:

    lund_tune/atlas_tune/   minimum-bias tune: ATLAS (and ALICE) data, Sherpa runs at 0.9, 7, 13 TeV, Pythia
    lund_tune/lep_tune/     ALEPH at the Z pole: Sherpa run, Pythia, the ALEPH reference data
    lund_tune/deuterons/    ALICE, BaBar, LEP deuterons: Sherpa runs, Pythia runs, HEPData tables

    python3 collect_lund_tune.py

The copied deuteron cards get OUTPUT: sherpa (or sherpa/<seed>), so that
    cd lund_tune/deuterons/<setup>;  Sherpa -f Sherpa.yaml
writes next to the card.  Nothing else in the cards changes.  The originals stay where they are."""
import glob
import os
import shutil

from paths import DOC, TUNING, NON_DIFFRACTIVE, LUND_TUNE, DEUTERON_DATA, LEP_RIVET, PYTHIA_REFERENCE, card_set

RIVET_REF = "/Users/ekremorhandemiray/code_software/rivetbootstrap/local/share/Rivet"
TARGET = LUND_TUNE


def copy_file(source, target):
    os.makedirs(os.path.dirname(target), exist_ok=True)
    shutil.copy2(source, target)


def copy_dir(source, target):
    if os.path.exists(target):
        shutil.rmtree(target)
    shutil.copytree(source, target)


def card_with_output(text, output, seed):
    "the card text with OUTPUT replaced, and a RANDOM_SEED line before EVENTS if seed is given"
    out = []
    for line in text.rstrip("\n").splitlines():
        # the original says where to run it from (SimpleLEP); the copy runs in its own folder
        if line.startswith("# Run from "):
            out.append("# Run from this folder:  Sherpa -f <this card>")
            continue
        if line.strip().startswith("OUTPUT:"):
            indent = line[:len(line) - len(line.lstrip())]
            out.append(indent + "OUTPUT:     " + output)
            continue
        if seed is not None and line.startswith("EVENTS:"):
            out.append("RANDOM_SEED:    %d" % seed)
        out.append(line)
    return "\n".join(out) + "\n"


def write(path, text):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write(text)


def atlas_tune():
    base = os.path.join(TARGET, "atlas_tune")
    for name in sorted(glob.glob(os.path.join(TUNING, "data", "*.dat"))):
        copy_file(name, os.path.join(base, "data", os.path.basename(name)))
    copy_file(os.path.join(TUNING, "pythia_minbias.cc"), os.path.join(base, "pythia_minbias.cc"))
    copy_file(os.path.join(TUNING, "cards", "base.yaml"), os.path.join(base, "base.yaml"))
    for name in ["fit_nd_round6.txt", "fit_nd_round6.json"]:
        copy_file(os.path.join(NON_DIFFRACTIVE, name), os.path.join(base, "fit", name))
    for energy in [900, 7000, 13000]:
        run = os.path.join(TUNING, "runs", "validate6", "p018", str(energy))
        folder = os.path.join(base, str(energy))
        copy_file(os.path.join(run, "card.yaml"), os.path.join(folder, "Sherpa.yaml"))
        copy_file(os.path.join(run, "out.json"), os.path.join(folder, "sherpa", "out.json"))
        for sample in ["nondiffractive", "inelastic"]:
            source = os.path.join(TUNING, "runs", "pythia", "%s_%d.json" % (sample, energy))
            copy_file(source, os.path.join(folder, "pythia", sample + ".json"))


def lep_tune():
    base = os.path.join(TARGET, "lep_tune")
    # the card without the EVENT_OUTPUT line: it is given on the command line
    lines = []
    for line in open(os.path.join(LEP_RIVET, "lund_mbtune", "card.yaml")).read().rstrip("\n").splitlines():
        if line.startswith("EVENT_OUTPUT:"):
            continue
        lines.append(line)
    write(os.path.join(base, "Sherpa.yaml"), "\n".join(lines) + "\n")
    copy_file(os.path.join(LEP_RIVET, "lund_mbtune", "lep.json"), os.path.join(base, "sherpa", "lep.json"))
    copy_file(os.path.join(LEP_RIVET, "pythia", "lep.json"), os.path.join(base, "pythia", "lep.json"))
    copy_file(os.path.join(TUNING, "pythia_lep.cc"), os.path.join(base, "pythia_lep.cc"))
    for analysis in ["ALEPH_1991_I319520", "ALEPH_1995_I382179"]:
        name = analysis + ".yoda.gz"
        copy_file(os.path.join(RIVET_REF, name), os.path.join(base, "data", name))


def deuterons():
    base = os.path.join(TARGET, "deuterons")
    cards = card_set("Lund_mbtune")
    sherpa_runs = os.path.join(DEUTERON_DATA, "lund_mbtune")
    pythia_runs = os.path.join(DEUTERON_DATA, "pythia")
    copy_dir(os.path.join(DEUTERON_DATA, "data"), os.path.join(base, "data"))
    copy_file(os.path.join(PYTHIA_REFERENCE, "pythia_deuterons.cc"), os.path.join(base, "pythia_deuterons.cc"))

    # one card per setup; ALICE 7 TeV = the card's own run + seeds 7001-7003 (run_lund.py), merged
    setups = []
    setups.append(("ALICE_0.9TeV", "Sherpa.Deuterons.ALICE_0.9TeV.Lund.yaml"))
    setups.append(("ALICE_2.76TeV", "Sherpa.Deuterons.ALICE_2.76TeV.Lund.yaml"))
    setups.append(("BaBar_10.58GeV", "Sherpa.Deuterons.BaBar.Lund.yaml"))
    setups.append(("LEP", "Sherpa.Deuterons.LEP.Lund.yaml"))
    for (setup, card) in setups:
        folder = os.path.join(base, setup)
        text = open(os.path.join(cards, card)).read()
        write(os.path.join(folder, "Sherpa.yaml"), card_with_output(text, "sherpa", None))
        copy_dir(os.path.join(sherpa_runs, setup), os.path.join(folder, "sherpa"))
        copy_dir(os.path.join(pythia_runs, setup), os.path.join(folder, "pythia"))

    folder = os.path.join(base, "ALICE_7TeV")
    text = open(os.path.join(cards, "Sherpa.Deuterons.ALICE_7TeV.Lund.yaml")).read()
    write(os.path.join(folder, "Sherpa.yaml"), card_with_output(text, "sherpa_seeds/card_seed", None))
    for seed in [7001, 7002, 7003]:
        write(os.path.join(folder, "Sherpa.seed%d.yaml" % seed), card_with_output(text, "sherpa_seeds/seed%d" % seed, seed))
    copy_dir(os.path.join(sherpa_runs, "ALICE_7TeV_seeds"), os.path.join(folder, "sherpa_seeds"))
    copy_dir(os.path.join(sherpa_runs, "ALICE_7TeV"), os.path.join(folder, "sherpa"))
    copy_dir(os.path.join(pythia_runs, "ALICE_7TeV"), os.path.join(folder, "pythia"))

    folder = os.path.join(base, "ALICE_7TeV_multiplicity")
    for seed in [1001, 1002, 1003, 1004, 1005]:
        card = "Sherpa.Deuterons.ALICE_7TeV_multiplicity.seed%d.Lund.yaml" % seed
        text = open(os.path.join(cards, card)).read()
        write(os.path.join(folder, "Sherpa.seed%d.yaml" % seed), card_with_output(text, "sherpa/seed_%d" % seed, None))
    copy_dir(os.path.join(sherpa_runs, "ALICE_7TeV_multiplicity"), os.path.join(folder, "sherpa"))
    copy_dir(os.path.join(pythia_runs, "ALICE_7TeV_multiplicity"), os.path.join(folder, "pythia"))


if __name__ == "__main__":
    atlas_tune()
    lep_tune()
    deuterons()
    print("copied into", TARGET)
