"""How much does each measurement depend on diffraction?  Pythia 8 (Monash) is run twice at each
energy: all inelastic events, and non-diffractive events only.  For every measurement this
prints the average over the fit bins of (non-diffractive only) / (all inelastic), and the
largest deviation from 1.  Measurements used in the tune should be close to 1.

    python3 diffraction_check.py  > ../minimum_bias_tuning/diffraction_check.txt
"""
import os

import combine
from observables import OBSERVABLES
from paths import TUNING

print("%-24s %-6s %10s %10s %14s" % ("measurement", "tune", "ND/inel", "max |dev|", "Pythia chi2/bin"))
for entry in OBSERVABLES:
    name = entry["name"]
    energy = entry["energy"]
    inelastic = combine.combine([(combine.load(os.path.join(TUNING, "runs", "pythia", "inelastic_%d.json" % energy)), 1.0)])
    nondiffractive = combine.combine([(combine.load(os.path.join(TUNING, "runs", "pythia", "nondiffractive_%d.json" % energy)), 1.0)])
    lows, highs, data, data_errors = combine.read_data(name)
    ratios = []
    for i in range(len(data)):
        if not combine.in_fit_range(entry, lows[i], highs[i], data[i]):
            continue
        if inelastic[name][2][i] <= 0:
            continue
        ratios.append(nondiffractive[name][2][i] / inelastic[name][2][i])
    average = sum(ratios) / len(ratios)
    largest = 0.0
    for r in ratios:
        if abs(r - 1) > largest:
            largest = abs(r - 1)
    c2, n = combine.chi2(name, inelastic)
    print("%-24s %-6s %10.3f %10.3f %14.2f" % (name, entry["tune"], average, largest, c2 / n))
