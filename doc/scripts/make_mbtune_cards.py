"""Copies the Lund_tune1 cards to Lund_mbtune with the final minimum-bias tune put in
(doc/minimum_bias_tuning, run_mb.py round 6, best run p018, checked with twice the events in
validate6).  Everything of that tune goes in, not only the tuned values:
  AMISIC     SIGMA_ND_NORM, PT_0(ref), PT_Min(ref), Eta, MU_R_FACTOR      (tuned)
             PT_Min(IR) 0.3, nB_bins 400                                  (fixed in the tune)
  REMNANTS   MATTER_RADIUS_1, MATTER_FRACTION_1                           (tuned)
  COLOUR_RECONNECTIONS  ETA_Q, RESHUFFLE                                  (tuned)
  MPI_PDF_SET NNPDF23_lo_as_0130_qed, SHOWER: MI_KFACTOR_SCHEME 1         (fixed in the tune)
The Lund strings and the deuteron model stay as in Lund_tune1.  The e+e- cards (LEP, BaBar)
get the same values, so that all cards are one setting.  Only the colour reconnections act
there (ETA_Q 0.568 -> 0.625, RESHUFFLE 0.333 -> 0.330); AMISIC, the matter profile, the MPI PDF
and the MPI shower do not, as there are no multiple scatterings in e+e- events.

    python3 make_mbtune_cards.py <run_mb.py round> <point>        e.g.  validate6 p018

The OUTPUT directory lund_tune1 -> lund_mbtune; run them with  python3 run_lund.py mbtune."""
import glob
import os
import sys

from paths import card_set
from run_mb import ROUNDS


def tune_values(round_name, point):
    "{run-card path: value} of the point, the round's fixed settings included"
    values = {}
    fixed = ROUNDS[round_name]["fixed"]
    for path in fixed:
        values[path] = fixed[path]
    point_values = ROUNDS[round_name]["extra_points"][point]
    for path in point_values:
        values[path] = point_values[path]
    return values


def new_card(text, values, label):
    "the Lund_tune1 card text with the minimum-bias tune; and whether it has an AMISIC block"
    out = []
    in_block = ""
    is_pp = False
    has_shower = False
    for line in text.rstrip("\n").splitlines():
        line = line.replace("lund_tune1", "lund_mbtune").replace("Lund_tune1/", "Lund_mbtune/")
        stripped = line.strip()
        if not line.startswith(" ") and stripped.endswith(":"):
            in_block = stripped[:-1]
        if in_block == "AMISIC":
            is_pp = True
        # a card with its own SHOWER block (LEP, BaBar) gets the K-factor inside it
        if line == "SHOWER:":
            has_shower = True
            out.append(line)
            out.append("  MI_KFACTOR_SCHEME: %d   # fixed in the minimum-bias tune" % values["SHOWER:MI_KFACTOR_SCHEME"])
            continue
        if in_block == "AMISIC" and stripped.startswith("SIGMA_ND_NORM:"):
            out.append("  SIGMA_ND_NORM:  %.4f         # minimum-bias tune (Lund_tune1: 0.614)" % values["AMISIC:SIGMA_ND_NORM"])
            continue
        if in_block == "AMISIC" and stripped.startswith('"PT_0(ref)":'):
            out.append('  "PT_0(ref)":    %.4f         # minimum-bias tune (Lund_tune1: 1.21)' % values["AMISIC:PT_0(ref)"])
            continue
        if in_block == "AMISIC" and stripped.startswith('"PT_Min(ref)":'):
            out.append('  "PT_Min(ref)":  %.4f         # minimum-bias tune (Lund_tune1: 1.1)' % values["AMISIC:PT_Min(ref)"])
            out.append('  Eta:            %.4f         # minimum-bias tune (Lund_tune1: default 0.026)' % values["AMISIC:Eta"])
            out.append('  MU_R_FACTOR:    %.4f         # minimum-bias tune (Lund_tune1: default 0.5)' % values["AMISIC:MU_R_FACTOR"])
            out.append('  "PT_Min(IR)":   %.1f            # fixed in the minimum-bias tune' % values["AMISIC:PT_Min(IR)"])
            out.append('  nB_bins:        %d            # b table: integrates to 1.001 (default 20: 1.255)' % values["AMISIC:nB_bins"])
            continue
        if in_block == "REMNANTS" and stripped.startswith("MATTER_RADIUS_1:"):
            out.append("    MATTER_RADIUS_1:   %.4f    # minimum-bias tune (Lund_tune1: 0.85)" % values["REMNANTS:2212:MATTER_RADIUS_1"])
            continue
        if in_block == "REMNANTS" and stripped.startswith("MATTER_FRACTION_1:"):
            out.append("    MATTER_FRACTION_1: %.4f    # minimum-bias tune (Lund_tune1: 0.65)" % values["REMNANTS:2212:MATTER_FRACTION_1"])
            continue
        if in_block == "COLOUR_RECONNECTIONS" and stripped.startswith("ETA_Q:"):
            out.append("  ETA_Q:     %.4f              # minimum-bias tune (Lund_tune1: 0.568)" % values["COLOUR_RECONNECTIONS:ETA_Q"])
            continue
        if in_block == "COLOUR_RECONNECTIONS" and stripped.startswith("RESHUFFLE:"):
            out.append("  RESHUFFLE: %.4f              # minimum-bias tune (Lund_tune1: 0.333)" % values["COLOUR_RECONNECTIONS:RESHUFFLE"])
            continue
        out.append(line)
    if is_pp:
        out.append("# fixed in the minimum-bias tune: Pythia's LO PDF (and alpha_s) for the scatterings,")
        out.append("# and the CMW K-factor in the shower of the scatterings, as in the signal shower")
        out.append("MPI_PDF_SET: %s" % values["MPI_PDF_SET"])
        if not has_shower:
            out.append("SHOWER:")
            out.append("  MI_KFACTOR_SCHEME: %d" % values["SHOWER:MI_KFACTOR_SCHEME"])
    out.insert(0, "# Lund_mbtune: Lund_tune1 with the minimum-bias tune " + label + " (make_mbtune_cards.py)")
    return "\n".join(out) + "\n", is_pp


def main(round_name, point):
    values = tune_values(round_name, point)
    target = card_set("Lund_mbtune")
    os.makedirs(target, exist_ok=True)
    cards = sorted(glob.glob(os.path.join(card_set("Lund_tune1"), "*.Lund.yaml")))
    for card in cards:
        (text, is_pp) = new_card(open(card).read(), values, round_name + " " + point)
        # in a pp card every value must have been replaced
        if is_pp:
            checks = []
            checks.append("SIGMA_ND_NORM:  %.4f" % values["AMISIC:SIGMA_ND_NORM"])
            checks.append("MU_R_FACTOR:    %.4f" % values["AMISIC:MU_R_FACTOR"])
            checks.append("MATTER_RADIUS_1:   %.4f" % values["REMNANTS:2212:MATTER_RADIUS_1"])
            checks.append("ETA_Q:     %.4f" % values["COLOUR_RECONNECTIONS:ETA_Q"])
            checks.append("RESHUFFLE: %.4f" % values["COLOUR_RECONNECTIONS:RESHUFFLE"])
            for key in checks:
                assert key in text, (card, key)
        with open(os.path.join(target, os.path.basename(card)), "w") as f:
            f.write(text)
    print(len(cards), "cards ->", target)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
