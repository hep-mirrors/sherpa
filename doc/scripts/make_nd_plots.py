"""Plots of the non-diffractive minimum-bias tune, one PDF per measurement, into
minimum_bias_tuning/non_diffractive_events/figures/.

    python3 make_nd_plots.py            first candidates  -> figures/
    python3 make_nd_plots.py final      the final tune    -> figures_final/

Every plot shows the data and, all with non-diffractive events only (runs in minimum_bias_tuning/):
  first candidates:
    Sherpa Lund_tune1                     runs/round1/lund_tune1
    Sherpa Lund_tune1 with default MPI    runs/round1/sherpa_default
    candidate A: polynomial minimum       runs/validate/poly_best
    candidate B: best run point p004      runs/validate/p004
    Pythia 8 Monash                       runs/pythia/nondiffractive_<E>.json
  final (one Sherpa setting only, the one used at LEP and in pp):
    Sherpa final minimum-bias tune        runs/validate6/p018 (round 6, twice the events)
    Pythia 8 Monash                       runs/pythia/nondiffractive_<E>.json
with MC/data below, and chi^2/bin (in the fit range of combine.py) in the legend.

The measurements not used in the tune with n_ch >= 1, INEL>0 or per inelastic event
contain diffractive events in the data; the curves here do not.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from paths import TUNING, NON_DIFFRACTIVE

import combine
from observables import OBSERVABLES, ENERGIES
from make_mb_plots import KIND_TITLE, KIND_X, title_of

FIGURES = os.path.join(NON_DIFFRACTIVE, "figures")

# (label, run directory or None for Pythia, colour, line style)
SOURCES = []
SOURCES.append(("Sherpa Lund_tune1", os.path.join(TUNING, "runs", "round1", "lund_tune1"), "#E69F00", "-"))
SOURCES.append(("Sherpa Lund_tune1 + default MPI", os.path.join(TUNING, "runs", "round1", "sherpa_default"), "0.5", ":"))
SOURCES.append(("Sherpa candidate A (fit minimum)", os.path.join(TUNING, "runs", "validate", "poly_best"), "#CC0000", "-"))
SOURCES.append(("Sherpa candidate B (p004)", os.path.join(TUNING, "runs", "validate", "p004"), "#009E73", "--"))
SOURCES.append(("Pythia 8 Monash", None, "#0072B2", "-."))

FINAL_SOURCES = []
FINAL_SOURCES.append(("Sherpa final tune (round 6, p018)", os.path.join(TUNING, "runs", "validate6", "p018"), "#CC0000", "-"))
FINAL_SOURCES.append(("Pythia 8 Monash", None, "#0072B2", "-."))


def predictions(sources):
    """[(label, {observable: (lows, highs, values, errors)}, colour, style)]"""
    result = []
    for (label, directory, color, style) in sources:
        values = {}
        for energy in ENERGIES:
            if directory is None:
                out = os.path.join(TUNING, "runs", "pythia", "nondiffractive_%d.json" % energy)
            else:
                out = os.path.join(directory, str(energy), "out.json")
            one_energy = combine.combine([(combine.load(out), 1.0)])
            for name in one_energy:
                values[name] = one_energy[name]
        result.append((label, values, color, style))
    return result


def plot(entry, generators, figures):
    name = entry["name"]
    lows, highs, data, data_errors = combine.read_data(name)
    centres = []
    for i in range(len(lows)):
        centres.append(0.5 * (lows[i] + highs[i]))

    figure = plt.figure(figsize=(7.0, 6.0))
    top = figure.add_axes([0.13, 0.36, 0.83, 0.56])
    bottom = figure.add_axes([0.13, 0.10, 0.83, 0.24], sharex=top)

    # data, and its relative error band in the ratio panel
    top.errorbar(centres, data, yerr=data_errors, fmt="o", color="black", ms=3.5, label="data", zorder=5)
    band_low = []
    band_high = []
    for i in range(len(data)):
        if data[i] > 0:
            band_low.append(1 - data_errors[i] / data[i])
            band_high.append(1 + data_errors[i] / data[i])
        else:
            band_low.append(1.0)
            band_high.append(1.0)
    edges = list(lows) + [highs[-1]]
    band_low_steps = band_low + [band_low[-1]]
    band_high_steps = band_high + [band_high[-1]]
    bottom.fill_between(edges, band_low_steps, band_high_steps, step="post", color="0.85")
    bottom.axhline(1.0, color="black", lw=0.8)

    for (label, values, color, style) in generators:
        mc_lows, mc_highs, mc, mc_errors = values[name]
        c2, n = combine.chi2(name, values)
        text = label
        if n > 0:
            text = label + r"  ($\chi^2$/bin %.1f)" % (c2 / n)
        steps = list(mc) + [mc[-1]]
        top.step(edges, steps, where="post", color=color, ls=style, lw=1.6, label=text)
        ratio = []
        for i in range(len(mc)):
            if data[i] > 0 and mc[i] > 0:
                ratio.append(mc[i] / data[i])
            else:
                ratio.append(float("nan"))
        ratio_steps = ratio + [ratio[-1]]
        bottom.step(edges, ratio_steps, where="post", color=color, ls=style, lw=1.6)

    if entry["kind"] in ("pt", "nch", "ident"):
        top.set_yscale("log")
    if entry["kind"] in ("pt", "ue_nch", "ue_ptsum"):
        top.set_xscale("log")
    top.set_ylabel(KIND_TITLE[entry["kind"]])
    bottom.set_xlabel(KIND_X[entry["kind"]])
    bottom.set_ylabel("MC / data")
    # the pT spectra need more room: the final tune is 2x the data at 10 GeV
    if entry["kind"] == "pt":
        bottom.set_ylim(0.4, 2.6)
    else:
        bottom.set_ylim(0.5, 1.5)
    if entry["tune"]:
        use = "used in the tune"
    else:
        use = "not used in the tune (" + entry["reason"] + ")"
    top.set_title(title_of(entry) + "\n" + use + "; MC: non-diffractive events", fontsize=8.5)
    top.legend(fontsize=7, frameon=False)
    plt.setp(top.get_xticklabels(), visible=False)
    figure.savefig(os.path.join(figures, name + ".pdf"))
    plt.close(figure)


def main(which):
    if which == "final":
        figures = os.path.join(NON_DIFFRACTIVE, "figures_final")
        sources = FINAL_SOURCES
    else:
        figures = FIGURES
        sources = SOURCES
    os.makedirs(figures, exist_ok=True)
    generators = predictions(sources)
    for entry in OBSERVABLES:
        plot(entry, generators, figures)
    print(len(OBSERVABLES), "plots in", figures)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main("first")
