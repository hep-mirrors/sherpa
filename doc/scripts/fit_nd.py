"""The minimum-bias tune with Sherpa's non-diffractive events only.

    python3 fit_nd.py <round> [<round> ...]

Same method as fit_tune.py (Professor: a quadratic polynomial per bin, chi^2 per bin of
each measurement, summed over the measurements), with one difference: the Sherpa prediction
of a run point is its non-diffractive run alone (AMISIC EVENT_TYPE Perturbative), without
Sherpa's diffractive samples.

Why: Sherpa's diffractive events are far too hard (the diffractive system is split into a
quark and a diquark in a random direction in its rest frame, so they get up to M/2 of
transverse momentum; see minimum_bias_tuning/non_diffractive_events/SUMMARY.md).  No parameter changes this, so the diffractive
samples are left out.  The measurements used in the fit (tune = True in observables.py)
were chosen such that diffraction changes them by 0-7 % (Pythia, minimum_bias_tuning/diffraction_check.txt).

Reads the runs of minimum_bias_tuning/runs/<round>/; writes fit_nd_<rounds>.txt and
fit_nd_<rounds>.json into minimum_bias_tuning/non_diffractive_events/.
"""
import json
import os
import sys

import numpy as np
from scipy.optimize import minimize

from paths import TUNING, NON_DIFFRACTIVE

import combine
from observables import OBSERVABLES, ENERGIES
from run_mb import ROUNDS, points


def tune_observables():
    "the entries of observables.py used in the fit"
    result = []
    for entry in OBSERVABLES:
        if entry["tune"]:
            result.append(entry)
    return result


def load_point_nd(round_name, point_name):
    """MC values of one point at all energies, non-diffractive events only:
    {observable name: (lows, highs, values, errors)}, or None if a run is missing"""
    prediction = {}
    for energy in ENERGIES:
        out = os.path.join(TUNING, "runs", round_name, point_name, str(energy), "out.json")
        if not os.path.exists(out):
            return None
        # one sample, so its cross section does not matter: weight 1
        one_energy = combine.combine([(combine.load(out), 1.0)])
        for name in one_energy:
            prediction[name] = one_energy[name]
    return prediction


def scale_to_unit(values, parameters):
    "parameter values -> x in [-1, 1] over the round-1 range"
    x = []
    for k in range(len(parameters)):
        (path, low, high) = parameters[k]
        x.append(2.0 * (values[path] - low) / (high - low) - 1.0)
    return x


def scale_from_unit(x, parameters):
    values = {}
    for k in range(len(parameters)):
        (path, low, high) = parameters[k]
        values[path] = low + (x[k] + 1.0) * (high - low) / 2.0
    return values


def polynomial_terms(x):
    "1, x_i, x_i x_j (i <= j): the terms of the quadratic polynomial"
    terms = [1.0]
    for i in range(len(x)):
        terms.append(x[i])
    for i in range(len(x)):
        for j in range(i, len(x)):
            terms.append(x[i] * x[j])
    return terms


def sort_key_first(item):
    return item[0]


def main(round_names):
    parameters = ROUNDS[round_names[0]]["parameters"]
    observables = tune_observables()

    # --- the MC values of all finished points ------------------------------------------
    run_points = []          # (round, point name, parameter values, prediction)
    for round_name in round_names:
        for (point_name, values) in points(round_name):
            prediction = load_point_nd(round_name, point_name)
            if prediction is None:
                continue
            run_points.append((round_name, point_name, values, prediction))
    print(len(run_points), "finished points")

    # --- the bins of the fit ------------------------------------------------------------
    bins = []                # (observable index, bin index, data value, data error)
    n_bins_of = []           # number of bins of each observable
    for k in range(len(observables)):
        entry = observables[k]
        lows, highs, data, data_errors = combine.read_data(entry["name"])
        n = 0
        for i in range(len(data)):
            if combine.in_fit_range(entry, lows[i], highs[i], data[i]):
                bins.append((k, i, data[i], data_errors[i]))
                n += 1
        n_bins_of.append(n)
    print(len(bins), "bins in", len(observables), "measurements")

    # --- direct chi^2 of every run point --------------------------------------------------
    direct = []
    for (round_name, point_name, values, prediction) in run_points:
        total = 0.0
        per_observable = []
        for k in range(len(observables)):
            c2, n = combine.chi2(observables[k]["name"], prediction)
            per_observable.append(c2 / n)
            total += c2 / n
        direct.append((total, round_name, point_name, values, per_observable))
    direct.sort(key=sort_key_first)

    # --- Pythia 8 Monash, non-diffractive events, for comparison --------------------------
    pythia = {}
    for energy in ENERGIES:
        out = os.path.join(TUNING, "runs", "pythia", "nondiffractive_%d.json" % energy)
        one_energy = combine.combine([(combine.load(out), 1.0)])
        for name in one_energy:
            pythia[name] = one_energy[name]
    pythia_per_observable = []
    pythia_total = 0.0
    for k in range(len(observables)):
        c2, n = combine.chi2(observables[k]["name"], pythia)
        pythia_per_observable.append(c2 / n)
        pythia_total += c2 / n

    # --- the polynomial of every bin -------------------------------------------------------
    # only the random points: the named points (lund_tune1, sherpa_default) are compared directly
    design = []
    for (round_name, point_name, values, prediction) in run_points:
        if not point_name.startswith("p"):
            continue
        in_range = True
        x = scale_to_unit(values, parameters)
        for xi in x:
            if xi < -1.0001 or xi > 1.0001:
                in_range = False
        if in_range:
            design.append((x, prediction))
    matrix = []
    for (x, prediction) in design:
        matrix.append(polynomial_terms(x))
    matrix = np.array(matrix)
    print(len(design), "points for", matrix.shape[1], "coefficients per bin")

    coefficients = []        # per bin
    mc_error = []            # per bin: average MC error over the points
    residual_ratio = []      # per bin: rms(polynomial - MC) / average MC error
    for (k, i, data_value, data_error) in bins:
        name = observables[k]["name"]
        mc_values = []
        mc_errors = []
        for (x, prediction) in design:
            mc_values.append(prediction[name][2][i])
            mc_errors.append(prediction[name][3][i])
        mc_values = np.array(mc_values)
        solution = np.linalg.lstsq(matrix, mc_values, rcond=None)[0]
        coefficients.append(solution)
        average_error = float(np.mean(mc_errors))
        mc_error.append(average_error)
        fitted = matrix.dot(solution)
        rms = float(np.sqrt(np.mean((fitted - mc_values) ** 2)))
        if average_error > 0:
            residual_ratio.append(rms / average_error)
        else:
            residual_ratio.append(0.0)
    print("polynomial vs MC: median rms/(MC error) = %.2f, 90%% quantile = %.2f"
          % (np.median(residual_ratio), np.quantile(residual_ratio, 0.9)))

    # --- minimise the chi^2 of the polynomials ---------------------------------------------
    def chi2_of(x):
        terms = np.array(polynomial_terms(x))
        per_observable = [0.0] * len(observables)
        for b in range(len(bins)):
            (k, i, data_value, data_error) = bins[b]
            prediction = float(terms.dot(coefficients[b]))
            variance = data_error ** 2 + mc_error[b] ** 2 + (combine.MODEL_UNCERTAINTY * data_value) ** 2
            per_observable[k] += (prediction - data_value) ** 2 / variance
        total = 0.0
        for k in range(len(observables)):
            per_observable[k] = per_observable[k] / n_bins_of[k]
            total += per_observable[k]
        return total, per_observable

    def objective(x):
        return chi2_of(x)[0]

    # start from the 10 best run points, and from the centre
    starts = []
    for (total, round_name, point_name, values, per_observable) in direct[:10]:
        starts.append(scale_to_unit(values, parameters))
    starts.append([0.0] * len(parameters))
    bounds = []
    for k in range(len(parameters)):
        bounds.append((-1.0, 1.0))
    best = None
    for x0 in starts:
        clipped = []
        for xi in x0:
            clipped.append(min(1.0, max(-1.0, xi)))
        result = minimize(objective, clipped, method="L-BFGS-B", bounds=bounds)
        if best is None or result.fun < best.fun:
            best = result
    best_values = scale_from_unit(list(best.x), parameters)
    best_total, best_per_observable = chi2_of(list(best.x))

    # --- output ------------------------------------------------------------------------------
    lines = []
    lines.append("minimum-bias tune, non-diffractive events only, rounds " + " ".join(round_names))
    lines.append("%d points, %d bins in %d measurements" % (len(run_points), len(bins), len(observables)))
    lines.append("polynomial vs MC: median rms/(MC error) = %.2f, 90%% quantile = %.2f"
                 % (np.median(residual_ratio), np.quantile(residual_ratio, 0.9)))
    lines.append("")
    lines.append("best values (polynomial minimum):")
    for (path, low, high) in parameters:
        at_edge = ""
        value = best_values[path]
        if abs(value - low) < 0.01 * (high - low) or abs(value - high) < 0.01 * (high - low):
            at_edge = "   <- at the edge of the scan range"
        lines.append("  %-34s %.4f   (range %.3f - %.3f)%s" % (path, value, low, high, at_edge))
    lines.append("predicted total chi^2 (sum over %d measurements of chi^2/bin) = %.2f" % (len(observables), best_total))
    lines.append("")
    lines.append("chi^2 per bin of each measurement:")
    lines.append("  %-24s %7s | %7s | %10s | %14s | %7s" % ("", "poly", "best run", "lund_tune1", "default MPI", "Pythia"))
    reference = {}
    for (total, round_name, point_name, values, per_observable) in direct:
        if point_name in ("lund_tune1", "sherpa_default") and point_name not in reference:
            reference[point_name] = (total, per_observable)
    # rounds without these named points (round 3): take them from round 1
    for point_name in ("lund_tune1", "sherpa_default"):
        if point_name in reference:
            continue
        prediction = load_point_nd("round1", point_name)
        total = 0.0
        per_observable = []
        for k in range(len(observables)):
            c2, n = combine.chi2(observables[k]["name"], prediction)
            per_observable.append(c2 / n)
            total += c2 / n
        reference[point_name] = (total, per_observable)
    best_run = direct[0]
    for k in range(len(observables)):
        row = "  %-24s %7.2f | %8.2f | %10.2f | %14.2f | %7.2f" % (
            observables[k]["name"], best_per_observable[k], best_run[4][k],
            reference["lund_tune1"][1][k], reference["sherpa_default"][1][k], pythia_per_observable[k])
        lines.append(row)
    row = "  %-24s %7.2f | %8.2f | %10.2f | %14.2f | %7.2f" % (
        "TOTAL", best_total, best_run[0], reference["lund_tune1"][0], reference["sherpa_default"][0], pythia_total)
    lines.append(row)
    lines.append("  (best run = " + best_run[2] + ")")
    lines.append("")
    lines.append("the 15 best run points (direct chi^2, no polynomial):")
    for (total, round_name, point_name, values, per_observable) in direct[:15]:
        text = "  %-8s %-14s total %7.2f  " % (round_name, point_name, total)
        for (path, low, high) in parameters:
            text += " %s=%.3f" % (path.split(":")[-1], values[path])
        lines.append(text)
    output = "\n".join(lines)
    print(output)
    out_name = "fit_nd_" + "_".join(round_names)
    with open(os.path.join(NON_DIFFRACTIVE, out_name + ".txt"), "w") as f:
        f.write(output + "\n")
    with open(os.path.join(NON_DIFFRACTIVE, out_name + ".json"), "w") as f:
        json.dump({"best": best_values, "chi2": best_total}, f, indent=1)
    print("->", os.path.join(NON_DIFFRACTIVE, out_name + ".txt"))


if __name__ == "__main__":
    main(sys.argv[1:])
