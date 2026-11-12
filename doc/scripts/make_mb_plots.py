"""Plots of the minimum-bias tune, one PDF per measurement, into minimum_bias_tuning/figures/.

    python3 make_mb_plots.py <run directory of the new tune> [<label of the new tune>]

Every plot shows the data, Sherpa with Lund_tune1, Sherpa's default multiple-scattering
values, the new tune and Pythia 8 (Monash, all inelastic events), with MC/data below.
Sherpa = non-diffractive + diffractive runs, weighted with Sherpa's cross sections (combine.py).
Also: chi2_vs_parameters.pdf, the chi^2 of every scan point against each parameter.
"""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import combine
from observables import OBSERVABLES, ENERGIES

from paths import TUNING
FIGURES = os.path.join(TUNING, "figures")
DIFFRACTION_DIR = os.path.join(TUNING, "runs", "diffraction")

KIND_TITLE = {
    "eta": r"$1/N_{\rm ev}\;{\rm d}N_{\rm ch}/{\rm d}\eta$",
    "pt": r"$1/N_{\rm ev}\;1/(2\pi p_T)\;{\rm d}^2N_{\rm ch}/{\rm d}\eta\,{\rm d}p_T$ [GeV$^{-2}$]",
    "nch": r"$1/N_{\rm ev}\;{\rm d}N_{\rm ev}/{\rm d}n_{\rm ch}$",
    "ptnch": r"$\langle p_T\rangle$ [GeV]",
    "ue_nch": r"transverse region: $\langle {\rm d}^2N_{\rm ch}/{\rm d}\eta\,{\rm d}\phi\rangle$",
    "ue_ptsum": r"transverse region: $\langle {\rm d}^2\sum p_T/{\rm d}\eta\,{\rm d}\phi\rangle$ [GeV]",
    "ident": r"$1/N_{\rm ev}\;{\rm d}^2N/{\rm d}y\,{\rm d}p_T$ [GeV$^{-1}$]",
}
KIND_X = {
    "eta": r"$\eta$",
    "pt": r"$p_T$ [GeV]",
    "nch": r"$n_{\rm ch}$",
    "ptnch": r"$n_{\rm ch}$",
    "ue_nch": r"$p_T$ of the leading track [GeV]",
    "ue_ptsum": r"$p_T$ of the leading track [GeV]",
    "ident": r"$p_T$ [GeV]",
}
SPECIES_NAME = {211: r"$\pi^+ + \pi^-$", 321: r"$K^+ + K^-$", 2212: r"$p + \bar p$"}


def title_of(entry):
    "a one-line description of the measurement"
    energy = entry["energy"] / 1000.0
    if entry["kind"] == "ident":
        return "%s, %s at |y| < 0.5, %g TeV (per inelastic event)" % (
            entry["analysis"], SPECIES_NAME[entry["species"]], energy)
    if entry["kind"] in ("ue_nch", "ue_ptsum"):
        return "%s, underlying event, %g TeV, $p_T > %g$ GeV, $|\\eta| < %g$" % (
            entry["analysis"], energy, entry["ptmin"], entry["etamax"])
    selection = "$n_{\\rm ch} \\geq %d$" % entry["nchmin"]
    if entry["analysis"].startswith("ALICE"):
        selection = "INEL>0"
    return "%s, %g TeV, $p_T > %g$ GeV, $|\\eta| < %g$, %s" % (
        entry["analysis"], energy, entry["ptmin"], entry["etamax"], selection)


def predictions(new_tune_dir):
    """{generator label: {observable: (lows, highs, values, errors)}} for all energies"""
    sources = []
    sources.append(("Sherpa Lund_tune1", os.path.join(TUNING, "runs", "round1", "lund_tune1"), "sherpa", "#E69F00", "-"))
    sources.append(("Sherpa default MPI", os.path.join(TUNING, "runs", "round1", "sherpa_default"), "sherpa", "0.5", ":"))
    sources.append((NEW_LABEL, new_tune_dir, "sherpa", "#CC0000", "-"))
    sources.append(("Pythia 8 Monash", os.path.join(TUNING, "runs", "pythia"), "pythia", "#0072B2", "-."))
    result = []
    for (label, directory, generator, color, style) in sources:
        values = {}
        for energy in ENERGIES:
            if generator == "sherpa":
                out = os.path.join(directory, str(energy), "out.json")
                if not os.path.exists(out):
                    continue
                one_energy = combine.sherpa_prediction(out, DIFFRACTION_DIR, energy)
            else:
                out = os.path.join(directory, "inelastic_%d.json" % energy)
                if not os.path.exists(out):
                    continue
                one_energy = combine.combine([(combine.load(out), 1.0)])
            for name in one_energy:
                values[name] = one_energy[name]
        result.append((label, values, color, style))
    return result


def plot(entry, generators):
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
    bottom.fill_between(edges, band_low_steps, band_high_steps, step="post", color="0.85", label="data error")
    bottom.axhline(1.0, color="black", lw=0.8)

    for (label, values, color, style) in generators:
        if name not in values:
            continue
        mc_lows, mc_highs, mc, mc_errors = values[name]
        c2, n = combine.chi2(name, values)
        text = label
        if n > 0:
            text = label + r"  ($\chi^2$/bin %.1f)" % (c2 / n)
        steps = list(mc) + [mc[-1]]
        top.step(edges, steps, where="post", color=color, ls=style, lw=1.8, label=text)
        ratio = []
        for i in range(len(mc)):
            if data[i] > 0 and mc[i] > 0:
                ratio.append(mc[i] / data[i])
            else:
                ratio.append(float("nan"))
        ratio_steps = ratio + [ratio[-1]]
        bottom.step(edges, ratio_steps, where="post", color=color, ls=style, lw=1.8)

    if entry["kind"] in ("pt", "nch", "ident"):
        top.set_yscale("log")
    if entry["kind"] in ("pt", "ue_nch", "ue_ptsum"):
        top.set_xscale("log")
    top.set_ylabel(KIND_TITLE[entry["kind"]])
    bottom.set_xlabel(KIND_X[entry["kind"]])
    bottom.set_ylabel("MC / data")
    bottom.set_ylim(0.5, 1.5)
    if entry["tune"]:
        use = "used in the tune"
    else:
        use = "not used in the tune (" + entry["reason"] + ")"
    top.set_title(title_of(entry) + "\n" + use, fontsize=8.5)
    top.legend(fontsize=7.5, frameon=False)
    plt.setp(top.get_xticklabels(), visible=False)
    figure.savefig(os.path.join(FIGURES, name + ".pdf"))
    plt.close(figure)


def main(new_tune_dir):
    os.makedirs(FIGURES, exist_ok=True)
    generators = predictions(new_tune_dir)
    for entry in OBSERVABLES:
        plot(entry, generators)
        print("figures/" + entry["name"] + ".pdf")


NEW_LABEL = "Sherpa MB tune"

if __name__ == "__main__":
    if len(sys.argv) > 2:
        NEW_LABEL = sys.argv[2]
    main(sys.argv[1])
