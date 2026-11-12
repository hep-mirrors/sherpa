"""Accuracy of apprentice surrogates on held-back points, relative to the data uncertainty.

    python3 check_apprentice_surrogates.py <approx.json> <runs dir> <validation.txt> <data.json>

For every held-back point and bin: (surrogate - MC)/sigma, with sigma = the error used in the
tune (data.json: stat (+) syst (+) 5 %), and the relative deviation (surrogate - MC)/MC.
A surrogate is good enough for tuning if these are well below 1 sigma."""
import json, os, sys
import numpy as np
import apprentice
import apprentice.io as io

approx_file, runs_dir, validation_file, data_file = sys.argv[1:5]
binids, approximations = io.readApprox(approx_file)
data = json.load(open(data_file))
points = [line.split(",")[0].strip() for line in open(validation_file) if line.strip() and not line.startswith("#")]
pull, rel, obs = [], [], []
for point in points:
    d = os.path.join(runs_dir, point)
    params = io.read_paramsfile(os.path.join(d, "params.dat"))
    x = list(params.values())
    histos = io.read_histos(os.path.join(d, "sherpa.yoda"))
    for b, a in zip(binids, approximations):
        name, num = b.split("#")
        mc = histos[name][int(num)][2]
        pred = a.predict(x)
        pull.append((pred - mc)/data[b][1]); rel.append((pred - mc)/mc); obs.append(name)
pull, rel, obs = np.abs(pull), np.abs(rel), np.array(obs)
print(f"{len(points)} held-back points, {len(binids)} bins")
print(f"{'observable':22s} {'|dev|/sigma: median':>20s} {'90%':>6s} {'max':>6s}   {'|dev|/MC: median':>16s} {'90%':>6s}")
for name in sorted(set(obs)):
    m = obs == name
    print(f"{name:22s} {np.median(pull[m]):20.2f} {np.percentile(pull[m], 90):6.2f} {pull[m].max():6.2f}   "
          f"{100*np.median(rel[m]):15.1f}% {100*np.percentile(rel[m], 90):5.1f}%")
