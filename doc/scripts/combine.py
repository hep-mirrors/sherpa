"""From the raw sums of mb_analysis.py to the measured quantities, and the chi^2 against data.

Sherpa's minimum-bias prediction is the sum of its event classes, each weighted with its
cross section:

    non-diffractive (AMISIC EVENT_TYPE Perturbative)   sigma_ND
    single diffraction, proton A dissociates           sigma_SD,A   (EVENT_TYPE DiffractiveA)
    single diffraction, proton B dissociates           sigma_SD,B   (EVENT_TYPE DiffractiveB)
    double diffraction                                 sigma_DD     (EVENT_TYPE DiffractiveAB)

Elastic events leave no particle in any of the acceptances, so they never enter.
The diffractive cross sections are the ones Sherpa reports for these event types
(tests/xsecs/xsecs.txt), and sigma_ND = sigma_tot - sigma_el - sigma_SD,A - sigma_SD,B - sigma_DD
with Sherpa's total cross section (Donnachie-Landshoff, Schuler-Sjostrand Eq. 4, as in
AMISIC++/Tools/Hadronic_XSec_Calculator.C):  sigma_tot = 21.70 s^0.0808 + 56.08 s^-0.4525 mb.

An event of a sample with N events and cross section sigma has the weight sigma / N.
"""
import json
import math
import os

from observables import OBSERVABLES

from paths import TUNING
MODEL_UNCERTAINTY = 0.05      # 5 % added to every data point in the chi^2 (as for pi/K/p in ../examples/.../TuneScan)

# Sherpa's elastic and diffractive cross sections [mb] (tests/xsecs/xsecs.txt, read off the
# "XS" of a 10-event run of each EVENT_TYPE)
SHERPA_XSECS = {}
SHERPA_XSECS[900] = {"el": 12.7826, "SDA": 5.84894, "SDB": 5.84894, "DD": 6.4138}
SHERPA_XSECS[2760] = {"el": 16.0627, "SDA": 6.41332, "SDB": 6.41332, "DD": 7.93544}
SHERPA_XSECS[7000] = {"el": 19.3527, "SDA": 6.85004, "SDB": 6.85004, "DD": 9.25827}
SHERPA_XSECS[13000] = {"el": 21.8631, "SDA": 7.12616, "SDB": 7.12616, "DD": 10.1662}


def sigma_total(sqrts):
    "Sherpa's total pp cross section [mb] (Donnachie-Landshoff)"
    s = sqrts * sqrts
    return 21.70 * s ** 0.0808 + 56.08 * s ** (-0.4525)


def sherpa_cross_sections(sqrts):
    "sigma_ND, sigma_SD,A, sigma_SD,B, sigma_DD [mb]"
    x = SHERPA_XSECS[sqrts]
    sigma_nd = sigma_total(sqrts) - x["el"] - x["SDA"] - x["SDB"] - x["DD"]
    return {"ND": sigma_nd, "SDA": x["SDA"], "SDB": x["SDB"], "DD": x["DD"]}


def load(json_file):
    return json.load(open(json_file))


def combine(samples):
    """samples: list of (raw sums from mb_analysis.py, cross section).
    Returns {observable name: (low edges, high edges, values, errors)}."""
    # the weight of one event of each sample
    event_weights = []
    for (raw, sigma) in samples:
        event_weights.append(sigma / raw["events"])

    result = {}
    names = list(samples[0][0]["observables"].keys())
    for name in names:
        first = samples[0][0]["observables"][name]
        kind = first["kind"]
        edges = first["edges"]
        nbins = len(edges) - 1
        entry = None
        for e in OBSERVABLES:
            if e["name"] == name:
                entry = e

        lows = []
        highs = []
        values = []
        errors = []
        for i in range(nbins):
            lows.append(edges[i])
            highs.append(edges[i + 1])
            width = edges[i + 1] - edges[i]

            if kind in ("eta", "pt", "nch", "ident"):
                # (weighted sum in the bin) / (weighted number of selected events) / bin width
                sum_w = 0.0
                sum_w2 = 0.0
                selected = 0.0
                for k in range(len(samples)):
                    obs = samples[k][0]["observables"][name]
                    w = event_weights[k]
                    sum_w += w * obs["sum_w"][i]
                    sum_w2 += w * w * obs["sum_w2"][i]
                    selected += w * obs["selected"]
                if selected <= 0:
                    values.append(0.0)
                    errors.append(0.0)
                    continue
                value = sum_w / selected / width
                error = math.sqrt(sum_w2) / selected / width
                if kind == "pt":
                    # 1/(2 pi pT) d2N/deta dpT: the 1/pT is in the fill weight; per unit eta
                    eta_range = 2 * entry["etamax"]
                    value = value / (2 * math.pi * eta_range)
                    error = error / (2 * math.pi * eta_range)
                values.append(value)
                errors.append(error)
            else:
                # profiles: weighted average of y in the bin, error of the average
                sum_w = 0.0
                sum_wy = 0.0
                sum_wy2 = 0.0
                entries = 0.0
                for k in range(len(samples)):
                    obs = samples[k][0]["observables"][name]
                    w = event_weights[k]
                    sum_w += w * obs["sum_w"][i]
                    sum_wy += w * obs["sum_wy"][i]
                    sum_wy2 += w * obs["sum_wy2"][i]
                    entries += obs["sum_w"][i]
                if sum_w <= 0:
                    values.append(0.0)
                    errors.append(0.0)
                    continue
                mean = sum_wy / sum_w
                variance = sum_wy2 / sum_w - mean * mean
                if variance < 0:
                    variance = 0.0
                if entries > 1:
                    error = math.sqrt(variance / (entries - 1))
                else:
                    error = mean
                values.append(mean)
                errors.append(error)
        result[name] = (lows, highs, values, errors)
    return result


def sherpa_prediction(nd_json, diffraction_dir, sqrts):
    """Sherpa's minimum-bias prediction: the non-diffractive run nd_json plus the diffractive
    runs in diffraction_dir/<sqrts>/{SDA,SDB,DD}/out.json, weighted with their cross sections."""
    sigmas = sherpa_cross_sections(sqrts)
    samples = []
    samples.append((load(nd_json), sigmas["ND"]))
    for part in ("SDA", "SDB", "DD"):
        f = os.path.join(diffraction_dir, str(sqrts), part, "out.json")
        samples.append((load(f), sigmas[part]))
    return combine(samples)


def read_data(name):
    "data/<name>.dat: low edges, high edges, values, total errors"
    lows = []
    highs = []
    values = []
    errors = []
    for line in open(os.path.join(TUNING, "data", name + ".dat")):
        if line.startswith("#"):
            continue
        columns = line.split()
        lows.append(float(columns[0]))
        highs.append(float(columns[1]))
        values.append(float(columns[2]))
        errors.append(float(columns[3]))
    return lows, highs, values, errors


def in_fit_range(entry, low, high, data_value):
    """the bins used in the chi^2.  The measurements reach far into tails that 100k
    minimum-bias events do not populate; those bins are left out (the same for every
    generator setting)."""
    kind = entry["kind"]
    energy = entry["energy"]
    if data_value <= 0:
        return False
    if kind == "pt":
        return high <= 10.0
    if kind == "ident":
        return low >= 0.3 and high <= 3.0
    if kind == "ue_nch" or kind == "ue_ptsum":
        if energy == 900:
            return high <= 6.0
        return high <= 10.0
    if kind == "nch":
        # n_ch up to where the measured fraction falls below 1e-4 per unit n_ch
        return data_value >= 1e-4
    if kind == "ptnch":
        return high <= 60.0
    return True


def chi2(name, prediction):
    """chi^2 per bin of one measurement: (MC - data)^2 / (data error^2 + MC error^2 + (5% data)^2),
    over the bins in the fit range.  Returns (chi^2 sum, number of bins)."""
    entry = None
    for e in OBSERVABLES:
        if e["name"] == name:
            entry = e
    lows, highs, data, data_errors = read_data(name)
    mc_lows, mc_highs, mc, mc_errors = prediction[name]
    total = 0.0
    n = 0
    for i in range(len(data)):
        if not in_fit_range(entry, lows[i], highs[i], data[i]):
            continue
        difference = mc[i] - data[i]
        variance = data_errors[i] ** 2 + mc_errors[i] ** 2 + (MODEL_UNCERTAINTY * data[i]) ** 2
        total += difference * difference / variance
        n += 1
    return total, n
