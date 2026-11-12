"""Copies the Lund_<from> cards to Lund_<to> with Pythia string-fragmentation settings added under
PYTHIA8: PARAMETERS (e.g. the LEP Lund retune of lund_scan.py).  Nothing else changes.

    python3 make_lund_variant.py tune1 tune3 StringZ:aLund=0.75 StringZ:bLund=0.9 ...

The new set goes to paths.card_set("Lund_<to>").  The OUTPUT directory lund_<from> -> lund_<to>; run them with  python3 run_lund.py <to>."""
import glob, os, sys

from paths import card_set

def variant(text, old, new, settings):
    lines = [line.replace(f"lund_{old}", f"lund_{new}").replace(f"Lund_{old}/", f"Lund_{new}/")
             for line in text.rstrip("\n").splitlines()]
    assert lines[-2:] == ["PYTHIA8:", "  DECAYS: false"], "expected the cards to end with the PYTHIA8 block"
    lines.insert(-2, f"# Lund_{new}: Lund_{old} with Lund string parameters retuned to LEP (make_lund_variant.py)")
    lines.append("  PARAMETERS:")
    lines += [f"  - {name}: {value}" for name, value in settings]
    return "\n".join(lines) + "\n"

if __name__ == "__main__":
    if len(sys.argv) < 4: sys.exit(__doc__)
    old, new = sys.argv[1], sys.argv[2]
    settings = [arg.split("=") for arg in sys.argv[3:]]
    target = card_set(f"Lund_{new}")
    os.makedirs(target, exist_ok=True)
    cards = sorted(glob.glob(os.path.join(card_set(f"Lund_{old}"), "*.Lund.yaml")))
    for card in cards:
        with open(os.path.join(target, os.path.basename(card)), "w") as f:
            f.write(variant(open(card).read(), old, new, settings))
    print(len(cards), "cards ->", target)
