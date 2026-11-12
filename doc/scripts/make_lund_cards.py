"""Writes the Lund-string run cards: the frozen cards of a ScanTune setting with Sherpa's cluster
hadronisation (AHADIC) replaced by Pythia 8's Lund string fragmentation.

    python3 make_lund_cards.py s6_008      # ScanTune_s6_008/*.yaml -> Lund_s6_008/*.yaml
                                           # (folders of the card sets: paths.card_set)

Changes, nothing else:
  * the AHADIC block is removed (it only steers the cluster model);
  * FRAGMENTATION: Pythia8, with Pythia's default Lund parameters (Monash 2013, as our Pythia runs);
  * PYTHIA8: DECAYS: false, so that Sherpa (HADRONS++) still decays the hadrons, as in the
    ScanTune runs: only the hadronisation differs;
  * the OUTPUT directory scanTune_<setting> -> lund_<setting>.
Shower, multiple scatterings (AMISIC), colour reconnections, beams and the deuteron model stay as
in the ScanTune card.  Needs a Sherpa built with SHERPA_ENABLE_PYTHIA8=ON (build_lund at the top of the repository)."""
import glob, os, sys

from paths import card_set

def lund_card(text, setting):
    out, skip = [], False
    for line in text.splitlines():
        top_level = line[:1] not in (" ", "#", "") and ":" in line
        if top_level: skip = line.startswith("AHADIC:")
        if skip and not line.strip(): skip = False          # the block ends at the next blank line
        if skip: continue
        out.append(line.replace(f"scanTune_{setting}", f"lund_{setting}")
                       .replace(f"ScanTune_{setting}", f"Lund_{setting}")
                       .replace(".ScanTune.yaml", ".Lund.yaml"))
    out += ["# Lund string fragmentation (Pythia 8.318, default parameters) instead of AHADIC;",
            "# Sherpa's HADRONS++ decays the hadrons (make_lund_cards.py)",
            "FRAGMENTATION: Pythia8",
            "PYTHIA8:",
            "  DECAYS: false"]
    return "\n".join(out) + "\n"

if __name__ == "__main__":
    if len(sys.argv) != 2: sys.exit(__doc__)
    setting = sys.argv[1]
    source, target = card_set(f"ScanTune_{setting}"), card_set(f"Lund_{setting}")
    os.makedirs(target, exist_ok=True)
    for card in sorted(glob.glob(os.path.join(source, "*.ScanTune.yaml"))):
        name = os.path.basename(card).replace(".ScanTune.yaml", ".Lund.yaml")
        with open(os.path.join(target, name), "w") as f: f.write(lund_card(open(card).read(), setting))
        print(os.path.join(target, name))
