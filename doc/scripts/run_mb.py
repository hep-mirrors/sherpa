"""Make and run the Sherpa runs of the minimum-bias tune.

    python3 run_mb.py make <round>        write runs/<round>/<point>/<energy>/card.yaml for every point
    python3 run_mb.py run <round>         run every job that has not finished, PARALLEL at a time
    python3 run_mb.py status <round>      list the jobs and whether they are done
    python3 run_mb.py diffraction         make and run Sherpa's diffractive samples (one set for all points)

runs/ and cards/ are in doc/minimum_bias_tuning.  A job is one point at one energy: Sherpa writes its events into a pipe (events.fifo), and
mb_analysis.py reads them and writes out.json.  The rounds (parameter ranges, number of
points, events) are defined in ROUNDS below.
"""
import os
import random
import signal
import subprocess
import sys
import time

import yaml

from paths import SCRIPTS, TUNING
SHERPA = "/Users/ekremorhandemiray/code_software/sherpaMainBranch/build_lund/bin/Sherpa"
BASE_CARD = os.path.join(TUNING, "cards", "base.yaml")
PARALLEL = 7          # about half of the 15 cores, to keep the laptop cool (was 8)
# A job whose Sherpa has not made its first event after this many seconds is stopped.
# This happens when the hard cross section above PT_Min is smaller than the (normalised)
# non-diffractive one: AMISIC (Amisic.C, FirstMinBias: XSratio < 1) then rejects every
# event and Sherpa retries forever.  Such a point is marked with the file no_events.
NO_EVENT_TIMEOUT = 600

# the parameters of the tune: card path, lowest and highest value of round 1
#   SIGMA_ND_NORM  scales the non-diffractive cross section that the number of scatterings is
#                  divided by: higher = fewer scatterings per event
#   PT_0(ref)      damping of the scattering cross section below PT_0 at 7 TeV
#   PT_Min(ref)    no scatterings below PT_Min at 7 TeV
#   Eta            energy dependence: PT_0 and PT_Min ~ (E/7 TeV)^(2 Eta); Eta >= 0 only
#   MATTER_RADIUS_1, MATTER_FRACTION_1: proton matter = Gaussian of radius 1 (fraction 1)
#                  + Gaussian of radius 1.0 (the rest); sets how much more central collisions
#                  scatter than peripheral ones
ROUNDS = {}
ROUNDS["round1"] = {
    "points": 80,
    "seed": 1,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "parameters": [
        ["AMISIC:SIGMA_ND_NORM", 0.25, 1.0],
        ["AMISIC:PT_0(ref)", 1.0, 3.0],
        ["AMISIC:PT_Min(ref)", 0.8, 1.6],
        ["AMISIC:Eta", 0.0, 0.15],
        ["REMNANTS:2212:MATTER_RADIUS_1", 0.3, 1.0],
        ["REMNANTS:2212:MATTER_FRACTION_1", 0.2, 0.9],
    ],
    # named points, run as well: Lund_tune1, and Sherpa's default (TuneSets/MPI-CRon) values
    "extra_points": {
        "lund_tune1": {"AMISIC:SIGMA_ND_NORM": 0.614, "AMISIC:PT_0(ref)": 1.21, "AMISIC:PT_Min(ref)": 1.1,
                       "AMISIC:Eta": 0.026, "REMNANTS:2212:MATTER_RADIUS_1": 0.85,
                       "REMNANTS:2212:MATTER_FRACTION_1": 0.65},
        "sherpa_default": {"AMISIC:SIGMA_ND_NORM": 0.44, "AMISIC:PT_0(ref)": 2.05, "AMISIC:PT_Min(ref)": 1.1,
                           "AMISIC:Eta": 0.026, "REMNANTS:2212:MATTER_RADIUS_1": 0.85,
                           "REMNANTS:2212:MATTER_FRACTION_1": 0.65},
    },
}

# round 2: around the best region of round 1 (minimum_bias_tuning/non_diffractive_events/fit_nd_round1.txt:
# SIGMA_ND_NORM 0.72, PT_0 1.50, PT_Min 0.80 (lower edge), Eta 0.014, MATTER_RADIUS_1 1.0
# (upper edge, = MATTER_RADIUS_2: one Gaussian), MATTER_FRACTION_1 0.52).
# PT_Min goes below the round-1 edge; MATTER_RADIUS_1 stays <= MATTER_RADIUS_2 = 1.0.
# Fit with  python3 fit_nd.py round2 round1  (the round-1 points
# inside these ranges are used as well; lund_tune1 and sherpa_default come from round 1).
ROUNDS["round2"] = {
    "points": 60,
    "seed": 2,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "parameters": [
        ["AMISIC:SIGMA_ND_NORM", 0.4, 0.95],
        ["AMISIC:PT_0(ref)", 1.2, 2.2],
        ["AMISIC:PT_Min(ref)", 0.5, 1.1],
        ["AMISIC:Eta", 0.0, 0.08],
        ["REMNANTS:2212:MATTER_RADIUS_1", 0.6, 1.0],
        ["REMNANTS:2212:MATTER_FRACTION_1", 0.2, 0.9],
    ],
    "extra_points": {},
}

# validation: the polynomial minimum of rounds 2+1 (minimum_bias_tuning/non_diffractive_events/fit_nd_round2_round1.json)
# and the best run point (round2 p004), with twice the events
ROUNDS["validate"] = {
    "points": 0,
    "seed": 3,
    "events": {900: "200k", 7000: "100k", 13000: "100k"},
    "parameters": ROUNDS["round2"]["parameters"],
    "extra_points": {
        "poly_best": {"AMISIC:SIGMA_ND_NORM": 0.8853, "AMISIC:PT_0(ref)": 1.2, "AMISIC:PT_Min(ref)": 0.7628,
                      "AMISIC:Eta": 0.0, "REMNANTS:2212:MATTER_RADIUS_1": 0.7730,
                      "REMNANTS:2212:MATTER_FRACTION_1": 0.2475},
        "p004": {"AMISIC:SIGMA_ND_NORM": 0.6447, "AMISIC:PT_0(ref)": 1.4682, "AMISIC:PT_Min(ref)": 0.5216,
                 "AMISIC:Eta": 0.0022, "REMNANTS:2212:MATTER_RADIUS_1": 0.786,
                 "REMNANTS:2212:MATTER_FRACTION_1": 0.4229},
    },
}

# round 3: MPI with Pythia's LO PDF (NNPDF2.3 LO, alpha_s(MZ) = 0.130), instead of the default
# PDF4LHC21 (NNLO).  Tests on top of p004 (scratch, not saved) showed that with the LO PDF a
# positive Eta works (PT_0 2.3, Eta 0.1: 0.9 and 13 TeV underlying event right), which the NNLO
# PDF did not allow.  New compared with round 2:
#   PT_Min(IR) = 0.3 (default 1.0): the floor of PT_Min(E) = max(PT_Min(IR), PT_Min(ref)(E/7TeV)^(2 Eta));
#                at 1.0 it made the PT_Min(ref) < 1 values of round 2 do nothing
#   MU_R_FACTOR scanned: mu_R^2 = MU_R_FACTOR^2 (pT^2 + PT_0^2) in alpha_s of the scatterings;
#                the only setting that softened the high-pT tail in the tests (MU_R_FACTOR 1)
#   PT_0 and Eta higher, around the LO test values (Pythia Monash: pT0 2.28 at 7 TeV, Eta = 0.215/2)
# The points are not comparable with rounds 1 and 2 (other PDF): fit round 3 on its own.
ROUNDS["round3"] = {
    "points": 80,
    "seed": 4,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "fixed": {"MPI_PDF_SET": "NNPDF23_lo_as_0130_qed", "AMISIC:PT_Min(IR)": 0.3},
    "parameters": [
        ["AMISIC:SIGMA_ND_NORM", 0.4, 1.0],
        ["AMISIC:PT_0(ref)", 1.8, 3.0],
        ["AMISIC:PT_Min(ref)", 0.4, 1.4],
        ["AMISIC:Eta", 0.03, 0.18],
        ["AMISIC:MU_R_FACTOR", 0.5, 1.0],
        ["REMNANTS:2212:MATTER_RADIUS_1", 0.6, 1.0],
        ["REMNANTS:2212:MATTER_FRACTION_1", 0.2, 0.9],
    ],
    "extra_points": {},
}

# validation of round 3: the polynomial minimum of round 3 (minimum_bias_tuning/non_diffractive_events/fit_nd_round3.json)
# and the best run point (round3 p039), with twice the events
ROUNDS["validate3"] = {
    "points": 0,
    "seed": 5,
    "events": {900: "200k", 7000: "100k", 13000: "100k"},
    "fixed": ROUNDS["round3"]["fixed"],
    "parameters": ROUNDS["round3"]["parameters"],
    "extra_points": {
        "poly_best": {"AMISIC:SIGMA_ND_NORM": 0.9021, "AMISIC:PT_0(ref)": 1.8, "AMISIC:PT_Min(ref)": 0.7404,
                      "AMISIC:Eta": 0.0709, "AMISIC:MU_R_FACTOR": 0.9048,
                      "REMNANTS:2212:MATTER_RADIUS_1": 0.6051,
                      "REMNANTS:2212:MATTER_FRACTION_1": 0.6431},
        "p039": {"AMISIC:SIGMA_ND_NORM": 0.7973, "AMISIC:PT_0(ref)": 1.9599, "AMISIC:PT_Min(ref)": 0.4826,
                 "AMISIC:Eta": 0.0516, "AMISIC:MU_R_FACTOR": 0.9045,
                 "REMNANTS:2212:MATTER_RADIUS_1": 0.6711,
                 "REMNANTS:2212:MATTER_FRACTION_1": 0.8313},
    },
}

# round 4: same settings as round 3 (LO PDF, PT_Min(IR) = 0.3), around the round-3 minimum.
# The round-3 minimum is at two lower edges: PT_0(ref) = 1.80 and MATTER_RADIUS_1 = 0.61, so
# both ranges go lower.  MU_R_FACTOR (minimum 0.90) goes above 1: mu_R up to 1.2 times
# sqrt(pT^2 + PT_0^2), within the usual factor-2 scale variation.
# Fit with  python3 fit_nd.py round4 round3  (the round-3 points inside these ranges are used as well).
ROUNDS["round4"] = {
    "points": 60,
    "seed": 6,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "fixed": ROUNDS["round3"]["fixed"],
    "parameters": [
        ["AMISIC:SIGMA_ND_NORM", 0.6, 1.0],
        ["AMISIC:PT_0(ref)", 1.3, 2.3],
        ["AMISIC:PT_Min(ref)", 0.4, 1.1],
        ["AMISIC:Eta", 0.03, 0.15],
        ["AMISIC:MU_R_FACTOR", 0.6, 1.2],
        ["REMNANTS:2212:MATTER_RADIUS_1", 0.4, 1.0],
        ["REMNANTS:2212:MATTER_FRACTION_1", 0.3, 0.9],
    ],
    "extra_points": {},
}

# validation of round 4: the polynomial minimum of rounds 4+3 (minimum_bias_tuning/non_diffractive_events/fit_nd_round4_round3.json)
# and the best run point (round4 p027), with twice the events
ROUNDS["validate4"] = {
    "points": 0,
    "seed": 7,
    "events": {900: "200k", 7000: "100k", 13000: "100k"},
    "fixed": ROUNDS["round4"]["fixed"],
    "parameters": ROUNDS["round4"]["parameters"],
    "extra_points": {
        "poly_best": {"AMISIC:SIGMA_ND_NORM": 1.0, "AMISIC:PT_0(ref)": 1.3911, "AMISIC:PT_Min(ref)": 1.0367,
                      "AMISIC:Eta": 0.0696, "AMISIC:MU_R_FACTOR": 1.1326,
                      "REMNANTS:2212:MATTER_RADIUS_1": 0.4,
                      "REMNANTS:2212:MATTER_FRACTION_1": 0.3018},
        "p027": {"AMISIC:SIGMA_ND_NORM": 0.9432, "AMISIC:PT_0(ref)": 1.5798, "AMISIC:PT_Min(ref)": 1.0764,
                 "AMISIC:Eta": 0.0817, "AMISIC:MU_R_FACTOR": 1.1533,
                 "REMNANTS:2212:MATTER_RADIUS_1": 0.4191,
                 "REMNANTS:2212:MATTER_FRACTION_1": 0.3514},
    },
}

# shower test: in minimum-bias events every scattering, the hardest one too, is showered by the MPI
# shower (Initialization_Handler.C: Shower_Handler(..., id-1), type 1 = MPI), which has its own
# defaults (Initialization_Handler.C:281-287): no CMW K-factor (MI_KFACTOR_SCHEME 0) and alpha_s
# at 0.66 kT^2 (MI_FS_AS_FAC, MI_IS_AS_FAC), against K-factor on, FS_AS_FAC 1.0, IS_AS_FAC 0.5 for
# the shower of a signal process.  Less radiation leaves more of a mini-jet's pT in its leading
# hadrons.  Both points are the round-4 fit minimum (validate4 poly_best):
#   mi_kfactor   only the K-factor switched on in the MPI shower
#   mi_as_main   the MPI shower with the signal shower's settings
MI_KFACTOR = {}
for path in ROUNDS["validate4"]["extra_points"]["poly_best"]:
    MI_KFACTOR[path] = ROUNDS["validate4"]["extra_points"]["poly_best"][path]
MI_KFACTOR["SHOWER:MI_KFACTOR_SCHEME"] = 1
MI_AS_MAIN = {}
for path in MI_KFACTOR:
    MI_AS_MAIN[path] = MI_KFACTOR[path]
MI_AS_MAIN["SHOWER:MI_FS_AS_FAC"] = 1.0
MI_AS_MAIN["SHOWER:MI_IS_AS_FAC"] = 0.5
ROUNDS["shower_test"] = {
    "points": 0,
    "seed": 8,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "fixed": ROUNDS["round4"]["fixed"],
    "parameters": ROUNDS["round4"]["parameters"],
    "extra_points": {
        "mi_kfactor": MI_KFACTOR,
        "mi_as_main": MI_AS_MAIN,
    },
}

# b-table test: AMISIC tabulates sigma_hard O(b) at nB_bins + 1 log-spaced b (Matter_Overlap.C:217-219,
# default 20, MI_Parameters.C:122) and interpolates linearly between them (Lookup_Tables.C:232).
# With 20 bins the table integrates to 1.25 sigma_hard instead of sigma_hard (the CheckTables warning
# of every run), so the proton radius is fixed (FixKandSmin) on a wrong table.  200 bins: a
# numerical setting, no change of the model.  Point: mi_kfactor of the shower test.
NB_BINS_200 = {}
for path in MI_KFACTOR:
    NB_BINS_200[path] = MI_KFACTOR[path]
NB_BINS_200["AMISIC:nB_bins"] = 200
ROUNDS["bins_test"] = {
    "points": 0,
    "seed": 9,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "fixed": ROUNDS["round4"]["fixed"],
    "parameters": ROUNDS["round4"]["parameters"],
    "extra_points": {
        "nb200": NB_BINS_200,
    },
}

# round 5: the re-tune after the two fixes, both fixed in every card:
#   AMISIC:nB_bins 400           the b table integrates to 1.001 instead of 1.255 (bins_test: <n> = sigma_hard/sigma_ND,
#                                hard scatterings as Pythia's); a numerical setting, not tuned
#   SHOWER:MI_KFACTOR_SCHEME 1   CMW K-factor in the shower of the scatterings, as in the signal shower (shower_test)
# The values of rounds 1-4 partly made up for the 1.25 of the b table, so the ranges are broad again;
# with the correct table there are about a third fewer scatterings, so SIGMA_ND_NORM and PT_0 go lower.
# Fit with  python3 fit_nd.py round5  (the earlier rounds are not comparable: other b table).
ROUNDS["round5"] = {
    "points": 80,
    "seed": 10,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "fixed": {"MPI_PDF_SET": "NNPDF23_lo_as_0130_qed", "AMISIC:PT_Min(IR)": 0.3,
              "AMISIC:nB_bins": 400, "SHOWER:MI_KFACTOR_SCHEME": 1},
    "parameters": [
        ["AMISIC:SIGMA_ND_NORM", 0.3, 1.0],
        ["AMISIC:PT_0(ref)", 1.2, 2.6],
        ["AMISIC:PT_Min(ref)", 0.4, 1.4],
        ["AMISIC:Eta", 0.0, 0.15],
        ["AMISIC:MU_R_FACTOR", 0.5, 1.2],
        ["REMNANTS:2212:MATTER_RADIUS_1", 0.4, 1.0],
        ["REMNANTS:2212:MATTER_FRACTION_1", 0.2, 0.9],
    ],
    "extra_points": {},
}

# validation of round 5: the polynomial minimum (minimum_bias_tuning/non_diffractive_events/fit_nd_round5.txt)
# and the best run point (round5 p077), with twice the events
ROUNDS["validate5"] = {
    "points": 0,
    "seed": 11,
    "events": {900: "200k", 7000: "100k", 13000: "100k"},
    "fixed": ROUNDS["round5"]["fixed"],
    "parameters": ROUNDS["round5"]["parameters"],
    "extra_points": {
        "poly_best": {"AMISIC:SIGMA_ND_NORM": 0.8694, "AMISIC:PT_0(ref)": 1.4111, "AMISIC:PT_Min(ref)": 0.6672,
                      "AMISIC:Eta": 0.0602, "AMISIC:MU_R_FACTOR": 0.7458,
                      "REMNANTS:2212:MATTER_RADIUS_1": 0.5426,
                      "REMNANTS:2212:MATTER_FRACTION_1": 0.7075},
        "p077": {"AMISIC:SIGMA_ND_NORM": 0.8708, "AMISIC:PT_0(ref)": 1.6248, "AMISIC:PT_Min(ref)": 0.6344,
                 "AMISIC:Eta": 0.0561, "AMISIC:MU_R_FACTOR": 0.6124,
                 "REMNANTS:2212:MATTER_RADIUS_1": 0.8676,
                 "REMNANTS:2212:MATTER_FRACTION_1": 0.219},
    },
}

# round 6: round 5 with the colour reconnections tuned as well.  validate5 showed that the
# tune reaches the measured multiplicity only with too many scatterings (p077: 6.6 per event
# against Pythia's 4.2, 1.3-1.5x Pythia's above 4 GeV, the pT tail), i.e. too few charged particles
# per scattering.  With Lund strings (FRAGMENTATION: Pythia8) that number is set by which partons
# the colour reconnections join into strings.  ETA_Q 0.568 / RESHUFFLE 0.333 come from the early
# s6_008 ALICE scan and were never fitted to these data; LEP hardly constrains them (one scattering).
# Sherpa's own MPI + CR tune (arXiv:2606.24702, underlying-event data, no minimum bias) has
# ETA_Q 0.63, RESHUFFLE 1/9, Q_0 1 (kept fixed here, as there).
# PT_0(ref) <= 2.2: above 2.28 every round-5 point made no events at 900 GeV (the no_events files).
ROUNDS["round6"] = {
    "points": 100,
    "seed": 12,
    "events": {900: "100k", 7000: "50k", 13000: "50k"},
    "fixed": ROUNDS["round5"]["fixed"],
    "parameters": [
        ["AMISIC:SIGMA_ND_NORM", 0.3, 1.0],
        ["AMISIC:PT_0(ref)", 1.1, 2.2],
        ["AMISIC:PT_Min(ref)", 0.4, 1.4],
        ["AMISIC:Eta", 0.0, 0.15],
        ["AMISIC:MU_R_FACTOR", 0.5, 1.2],
        ["REMNANTS:2212:MATTER_RADIUS_1", 0.4, 1.0],
        ["REMNANTS:2212:MATTER_FRACTION_1", 0.2, 0.9],
        ["COLOUR_RECONNECTIONS:ETA_Q", 0.3, 0.9],
        ["COLOUR_RECONNECTIONS:RESHUFFLE", 0.05, 0.5],
    ],
    "extra_points": {},
}

# validation of round 6: the best run point (round6 p018) and the polynomial minimum
# (non_diffractive_events/fit_nd_round6.txt; MATTER_FRACTION_1 at the edge, 0.9), with twice the events
ROUNDS["validate6"] = {
    "points": 0,
    "seed": 13,
    "events": {900: "200k", 7000: "100k", 13000: "100k"},
    "fixed": ROUNDS["round6"]["fixed"],
    "parameters": ROUNDS["round6"]["parameters"],
    "extra_points": {
        "p018": {"AMISIC:SIGMA_ND_NORM": 0.8524, "AMISIC:PT_0(ref)": 1.5174, "AMISIC:PT_Min(ref)": 0.4983,
                 "AMISIC:Eta": 0.053, "AMISIC:MU_R_FACTOR": 0.7185,
                 "REMNANTS:2212:MATTER_RADIUS_1": 0.9604,
                 "REMNANTS:2212:MATTER_FRACTION_1": 0.6824,
                 "COLOUR_RECONNECTIONS:ETA_Q": 0.625, "COLOUR_RECONNECTIONS:RESHUFFLE": 0.3299},
        "poly_best": {"AMISIC:SIGMA_ND_NORM": 0.8406, "AMISIC:PT_0(ref)": 1.4867, "AMISIC:PT_Min(ref)": 0.6947,
                      "AMISIC:Eta": 0.0617, "AMISIC:MU_R_FACTOR": 0.7054,
                      "REMNANTS:2212:MATTER_RADIUS_1": 0.7004,
                      "REMNANTS:2212:MATTER_FRACTION_1": 0.9,
                      "COLOUR_RECONNECTIONS:ETA_Q": 0.6676, "COLOUR_RECONNECTIONS:RESHUFFLE": 0.2869},
    },
}

ENERGIES = [900, 7000, 13000]
DIFFRACTION_PARTS = ["SDA", "SDB", "DD"]
DIFFRACTION_EVENT_TYPES = {"SDA": "DiffractiveA", "SDB": "DiffractiveB", "DD": "DiffractiveAB"}
DIFFRACTION_EVENTS = "100k"


def points(round_name):
    """[(point name, {parameter: value})]: the extra points, then the random points.
    The random points are uniform in the ranges, from a fixed seed, so 'make' always
    gives the same points."""
    config = ROUNDS[round_name]
    result = []
    for name in config["extra_points"]:
        result.append((name, config["extra_points"][name]))
    generator = random.Random(config["seed"])
    for i in range(config["points"]):
        values = {}
        for (path, low, high) in config["parameters"]:
            value = low + (high - low) * generator.random()
            values[path] = round(value, 4)
        result.append(("p%03d" % (i + 1), values))
    return result


def set_value(card, path, value):
    "card['A']['B']... = value for path 'A:B:...'; the key 2212 is an integer in the card"
    keys = path.split(":")
    node = card
    for key in keys[:-1]:
        if key.isdigit():
            key = int(key)
        if key not in node:
            node[key] = {}
        node = node[key]
    node[keys[-1]] = value


def write_card(directory, energy, events, settings, comment):
    "a copy of the base card with the beam energy, number of events and settings"
    card = yaml.safe_load(open(BASE_CARD))
    card["BEAM_ENERGIES"] = energy // 2
    card["EVENTS"] = events
    for path in settings:
        set_value(card, path, settings[path])
    os.makedirs(directory, exist_ok=True)
    with open(os.path.join(directory, "card.yaml"), "w") as f:
        f.write("# " + comment + ", made by run_mb.py from cards/base.yaml\n")
        for path in settings:
            f.write("#   " + path + " = " + str(settings[path]) + "\n")
        yaml.safe_dump(card, f, sort_keys=False)


def make(round_name):
    config = ROUNDS[round_name]
    n = 0
    for (name, values) in points(round_name):
        # the settings of the card: the round's fixed settings (if any), then the point's values
        settings = {}
        if "fixed" in config:
            for path in config["fixed"]:
                settings[path] = config["fixed"][path]
        for path in values:
            settings[path] = values[path]
        for energy in ENERGIES:
            directory = os.path.join(TUNING, "runs", round_name, name, str(energy))
            write_card(directory, energy, config["events"][energy], settings,
                       "minimum-bias tune " + round_name + ", point " + name + ", sqrt(s) = " + str(energy) + " GeV")
            n += 1
    print(n, "cards in runs/" + round_name)


def make_diffraction():
    for energy in ENERGIES:
        for part in DIFFRACTION_PARTS:
            directory = os.path.join(TUNING, "runs", "diffraction", str(energy), part)
            settings = {"AMISIC:EVENT_TYPE": DIFFRACTION_EVENT_TYPES[part]}
            write_card(directory, energy, DIFFRACTION_EVENTS, settings,
                       "Sherpa diffractive sample " + part + ", sqrt(s) = " + str(energy) + " GeV")


def jobs(round_name):
    "[(directory, energy)] of every run of the round"
    result = []
    if round_name == "diffraction":
        for energy in ENERGIES:
            for part in DIFFRACTION_PARTS:
                result.append((os.path.join(TUNING, "runs", "diffraction", str(energy), part), energy))
        return result
    for (name, values) in points(round_name):
        for energy in ENERGIES:
            result.append((os.path.join(TUNING, "runs", round_name, name, str(energy)), energy))
    return result


def done(directory):
    return os.path.exists(os.path.join(directory, "out.json"))


def has_no_events(directory):
    "the job was stopped before, because Sherpa made no event"
    return os.path.exists(os.path.join(directory, "no_events"))


def first_event_made(directory):
    "Sherpa prints 'Event 1 (' into run.log when it has made its first event"
    log = os.path.join(directory, "run.log")
    if not os.path.exists(log):
        return False
    return "Event 1 (" in open(log, errors="replace").read()


def start(directory, energy):
    "the reader on the pipe, then Sherpa writing into it"
    fifo = os.path.join(directory, "events.fifo")
    if os.path.exists(fifo):
        os.remove(fifo)
    os.mkfifo(fifo)
    reader_log = open(os.path.join(directory, "reader.log"), "w")
    reader = subprocess.Popen([sys.executable, os.path.join(SCRIPTS, "mb_analysis.py"), "events.fifo", str(energy), "out.json"],
                              cwd=directory, stdout=reader_log, stderr=subprocess.STDOUT)
    sherpa_log = open(os.path.join(directory, "run.log"), "w")
    # own process group, so that caffeinate and Sherpa can be stopped together (stop_job)
    sherpa = subprocess.Popen(["caffeinate", "-i", SHERPA, "-f", "card.yaml", "EVENT_OUTPUT: HepMC3_GenEvent[events.fifo]"],
                              cwd=directory, stdout=sherpa_log, stderr=subprocess.STDOUT, start_new_session=True)
    return reader, sherpa


def stop_job(reader, sherpa):
    "kill caffeinate with Sherpa under it, and the reader"
    os.killpg(sherpa.pid, signal.SIGKILL)
    sherpa.wait()
    reader.kill()
    reader.wait()


def run(round_name):
    # the slow 13 TeV jobs first, so that the last jobs to finish are short ones
    todo = []
    for first_energy in [13000, 7000, 900]:
        for (directory, energy) in jobs(round_name):
            if energy == first_energy and not done(directory) and not has_no_events(directory):
                todo.append((directory, energy))
    print(len(todo), "jobs to run,", PARALLEL, "at a time", flush=True)
    running = {}
    start_time = {}
    while len(todo) > 0 or len(running) > 0:
        while len(todo) > 0 and len(running) < PARALLEL:
            (directory, energy) = todo.pop(0)
            running[directory] = start(directory, energy)
            start_time[directory] = time.time()
            print(time.strftime("%H:%M"), "started", os.path.relpath(directory, TUNING), flush=True)
        for directory in list(running.keys()):
            (reader, sherpa) = running[directory]
            waited = time.time() - start_time[directory]
            if sherpa.poll() is None and waited > NO_EVENT_TIMEOUT and not first_event_made(directory):
                # no event after NO_EVENT_TIMEOUT seconds: AMISIC rejects every event
                stop_job(reader, sherpa)
                open(os.path.join(directory, "no_events"), "w").write(
                    "no event after %d s: AMISIC rejects every event (XSratio < 1)\n" % NO_EVENT_TIMEOUT)
                os.remove(os.path.join(directory, "events.fifo"))
                print(time.strftime("%H:%M"), "NO EVENTS, stopped", os.path.relpath(directory, TUNING), flush=True)
                del running[directory]
                continue
            if sherpa.poll() is not None and reader.poll() is not None:
                if done(directory):
                    print(time.strftime("%H:%M"), "finished", os.path.relpath(directory, TUNING), flush=True)
                else:
                    print(time.strftime("%H:%M"), "FAILED", os.path.relpath(directory, TUNING), flush=True)
                os.remove(os.path.join(directory, "events.fifo"))
                del running[directory]
            elif sherpa.poll() is not None and sherpa.returncode != 0:
                # Sherpa failed: the reader would wait for events forever
                reader.kill()
                print(time.strftime("%H:%M"), "FAILED (Sherpa exit", sherpa.returncode, ")",
                      os.path.relpath(directory, TUNING), flush=True)
                del running[directory]
        time.sleep(10)
    print(time.strftime("%H:%M"), "all done", flush=True)


def status(round_name):
    n_done = 0
    all_jobs = jobs(round_name)
    for (directory, energy) in all_jobs:
        if done(directory):
            n_done += 1
        else:
            print("not done:", os.path.relpath(directory, TUNING))
    print(n_done, "of", len(all_jobs), "done")


if __name__ == "__main__":
    command = sys.argv[1]
    if command == "make":
        make(sys.argv[2])
    elif command == "run":
        run(sys.argv[2])
    elif command == "status":
        status(sys.argv[2])
    elif command == "diffraction":
        make_diffraction()
        run("diffraction")
