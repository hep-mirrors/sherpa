"""Turns finished minimum-bias tune points into apprentice input, so that apprentice can be
compared with fit_nd.py on the same points, bins, errors and weights.

    python3 make_apprentice_mb_inputs.py <output dir> <round> [<more rounds>]
        e.g. ../minimum_bias_tuning/apprentice/round4_round3 round4 round3

As in fit_nd.py: the parameter ranges are those of the first round, and points of the other
rounds are used only if they are inside these ranges.  Run it with apprentice's python (it
needs yoda), see ../examples/RemnantTests/SimpleLEP/TuneScan/apprentice/README.md.

writes
    <output dir>/runs/<point>/params.dat    the point's parameters
    <output dir>/runs/<point>/sherpa.yoda   the bins of the fit (combine.in_fit_range), non-diffractive events
    <output dir>/data.json                  the ATLAS values and errors, {"/MB/<measurement>#<bin>": [value, error]}
    <output dir>/weights.txt                one weight per measurement

The bins of one measurement are numbered 0, 1, 2, ... in the order of the data (only the bins
in the fit range), so the YODA bin edges are 0, 1, 2, ... and not the real ones.

The chi^2 is the one of fit_nd.py:
    chi^2 = sum over measurements of (1/n) sum over its n bins of (MC - data)^2 / (data error^2 + MC error^2 + (0.05 data)^2)
    data.json error = sqrt(data error^2 + (0.05 data)^2); the MC error comes from the error
    surrogates (app-tune2 -e); weight = sqrt(1/n), because apprentice multiplies each
    measurement's chi^2 by weight^2.
"""
import json
import math
import os
import sys

import yoda

import combine
from fit_nd import tune_observables, load_point_nd
from run_mb import ROUNDS, points


def short_name(path):
    "AMISIC:PT_0(ref) -> PT_0(ref), REMNANTS:2212:MATTER_RADIUS_1 -> MATTER_RADIUS_1"
    return path.split(":")[-1]


def bins_in_fit(entry):
    "the indices of the data bins in the fit range"
    lows, highs, data, data_errors = combine.read_data(entry["name"])
    indices = []
    for i in range(len(data)):
        if combine.in_fit_range(entry, lows[i], highs[i], data[i]):
            indices.append(i)
    return indices


def inside(values, parameters):
    "the point is inside the ranges of the first round"
    for (path, low, high) in parameters:
        if values[path] < low - 1e-6 or values[path] > high + 1e-6:
            return False
    return True


def main(out, round_names):
    parameters = ROUNDS[round_names[0]]["parameters"]
    observables = tune_observables()
    os.makedirs(os.path.join(out, "runs"), exist_ok=True)

    used = 0
    for round_name in round_names:
        for (point_name, values) in points(round_name):
            # only the random points, as in fit_nd.py
            if not point_name.startswith("p"):
                continue
            if not inside(values, parameters):
                continue
            prediction = load_point_nd(round_name, point_name)
            if prediction is None:
                continue
            target = os.path.join(out, "runs", round_name + "_" + point_name)
            os.makedirs(target, exist_ok=True)
            with open(os.path.join(target, "params.dat"), "w") as f:
                for (path, low, high) in parameters:
                    f.write(short_name(path) + " " + str(values[path]) + "\n")
            objects = []
            for entry in observables:
                indices = bins_in_fit(entry)
                mc_lows, mc_highs, mc, mc_errors = prediction[entry["name"]]
                edges = []
                for b in range(len(indices) + 1):
                    edges.append(float(b))
                estimate = yoda.Estimate1D(edges, "/MB/" + entry["name"])
                for b in range(len(indices)):
                    i = indices[b]
                    estimate.bin(b + 1).setVal(float(mc[i]))
                    estimate.bin(b + 1).setErr((-float(mc_errors[i]), float(mc_errors[i])))
                objects.append(estimate)
            yoda.write(objects, os.path.join(target, "sherpa.yoda"))
            used += 1

    data = {}
    weight_lines = []
    for entry in observables:
        lows, highs, values, errors = combine.read_data(entry["name"])
        indices = bins_in_fit(entry)
        for b in range(len(indices)):
            i = indices[b]
            # error = sqrt(data error^2 + (0.05 data)^2)
            error = math.sqrt(errors[i] ** 2 + (combine.MODEL_UNCERTAINTY * values[i]) ** 2)
            data["/MB/" + entry["name"] + "#" + str(b)] = [values[i], error]
        # weight^2 = 1/n: each measurement counts as its chi^2 per bin
        weight_lines.append("/MB/" + entry["name"] + " %.6f" % math.sqrt(1.0 / len(indices)))
    with open(os.path.join(out, "data.json"), "w") as f:
        json.dump(data, f, indent=1)
    with open(os.path.join(out, "weights.txt"), "w") as f:
        f.write("\n".join(weight_lines) + "\n")
    print(used, "points ->", os.path.join(out, "runs") + ", data.json, weights.txt")


if __name__ == "__main__":
    main(os.path.abspath(sys.argv[1]), sys.argv[2:])
