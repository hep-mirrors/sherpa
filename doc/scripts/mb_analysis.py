"""Fill the minimum-bias observables of observables.py from a HepMC3 (ASCII) event stream.
Used instead of Rivet, because the local Rivet 4.1.4 / YODA installation crashes when it
writes its results (same reason as ../examples/RemnantTests/SimpleLEP/LEPRivet).

    mkfifo events.fifo
    python3 mb_analysis.py events.fifo <sqrt(s) [GeV]> out.json &
    Sherpa -f card.yaml 'EVENT_OUTPUT: HepMC3_GenEvent[events.fifo]'

Particles: all final-state (status 1) charged particles.  Sherpa (and Pythia with
ParticleDecays:limitTau0 = on, tau0Max = 10 mm) keep K0S, Lambda, Sigma, Xi and Omega
undecayed, so their decay products never appear; this is the "primary particle"
definition of ATLAS and ALICE (mean lifetime above 30 ps, i.e. c tau > 10 mm).

The output holds raw sums only (no normalisation); combine.py turns them into the
measured quantities, after mixing non-diffractive and diffractive samples.
For every observable:
  histograms (kind eta, pt, nch, ident):
     sum_w[i]      sum of fill weights in bin i
     sum_w2[i]     sum of squared fill weights in bin i
     selected      number of events that pass the event selection
  profiles (kind ptnch, ue_nch, ue_ptsum):
     sum_w[i], sum_wy[i], sum_wy2[i]   for the average of y in bin i of x
"""
import json
import math
import os
import sys

from observables import observables_at

from paths import TUNING

# 3 x electric charge of the quarks d u s c b t
QUARK_CHARGE3 = {1: -1, 2: 2, 3: -1, 4: 2, 5: -1, 6: 2}
STRANGE_BARYONS = [3112, 3222, 3312, 3334]       # Sigma-, Sigma+, Xi-, Omega-: not counted by ATLAS at 13 TeV


def charge3(pid):
    "three times the electric charge of a PDG code (same function as LEPRivet/lep_observables.py)"
    a = abs(pid)
    if pid > 0:
        sign = 1
    else:
        sign = -1
    if a in (11, 13, 15):
        return -3 * sign
    if a in (12, 14, 16, 22, 21, 23, 25, 130, 310):
        return 0
    if a == 24:
        return 3 * sign
    if a > 1000000000:
        return 0
    n = a % 10000
    q1 = (n // 1000) % 10
    q2 = (n // 100) % 10
    q3 = (n // 10) % 10
    if q1 == 0:
        # meson q2 q3bar; an even (up-type) q2 is the quark, an odd one the antiquark
        if q2 == 0:
            return 0
        c = QUARK_CHARGE3[q2] - QUARK_CHARGE3[q3]
        if q2 % 2 == 0:
            return c * sign
        else:
            return -c * sign
    return (QUARK_CHARGE3[q1] + QUARK_CHARGE3[q2] + QUARK_CHARGE3[q3]) * sign


def read_edges(name):
    "bin edges of the measurement, from data/<name>.dat"
    edges = []
    for line in open("%s/data/%s.dat" % (TUNING, name)):
        if line.startswith("#"):
            continue
        columns = line.split()
        low = float(columns[0])
        high = float(columns[1])
        if len(edges) == 0:
            edges.append(low)
        edges.append(high)
    return edges


def find_bin(edges, x):
    "index i with edges[i] <= x < edges[i+1], or -1 outside"
    if x < edges[0] or x >= edges[-1]:
        return -1
    low = 0
    high = len(edges) - 1
    while high - low > 1:
        middle = (low + high) // 2
        if x >= edges[middle]:
            low = middle
        else:
            high = middle
    return low


def delta_phi(phi1, phi2):
    "|phi1 - phi2| folded into [0, pi]"
    d = abs(phi1 - phi2)
    if d > math.pi:
        d = 2 * math.pi - d
    return d


class Observable:
    def __init__(self, entry):
        self.entry = entry
        self.name = entry["name"]
        self.kind = entry["kind"]
        self.edges = read_edges(self.name)
        nbins = len(self.edges) - 1
        self.sum_w = [0.0] * nbins
        self.sum_w2 = [0.0] * nbins
        self.sum_wy = [0.0] * nbins
        self.sum_wy2 = [0.0] * nbins
        self.selected = 0

    def fill_histogram(self, x, weight):
        i = find_bin(self.edges, x)
        if i < 0:
            return
        self.sum_w[i] += weight
        self.sum_w2[i] += weight * weight

    def fill_profile(self, x, y):
        i = find_bin(self.edges, x)
        if i < 0:
            return
        self.sum_w[i] += 1.0
        self.sum_wy[i] += y
        self.sum_wy2[i] += y * y

    def result(self):
        out = {}
        out["kind"] = self.kind
        out["edges"] = self.edges
        out["selected"] = self.selected
        out["sum_w"] = self.sum_w
        out["sum_w2"] = self.sum_w2
        out["sum_wy"] = self.sum_wy
        out["sum_wy2"] = self.sum_wy2
        return out


def analyse_event(particles, observables):
    """particles: list of (pid, px, py, pz, e) of the final-state particles of one event.
    Fills every observable."""
    # the charged particles, with the quantities the cuts need
    charged = []
    for (pid, px, py, pz, e) in particles:
        if charge3(pid) == 0:
            continue
        pt = math.sqrt(px * px + py * py)
        if pt <= 0.0:
            continue
        p = math.sqrt(pt * pt + pz * pz)
        eta = 0.5 * math.log((p + pz) / (p - pz))
        phi = math.atan2(py, px)
        rapidity = 0.5 * math.log((e + pz) / (e - pz))
        charged.append((abs(pid), pt, eta, phi, rapidity))

    for obs in observables:
        entry = obs.entry
        kind = obs.kind

        if kind == "ident":
            # pi, K or p at |y| < 0.5, per inelastic event: every event counts
            obs.selected += 1
            for (apid, pt, eta, phi, rapidity) in charged:
                if apid == entry["species"] and abs(rapidity) < 0.5:
                    obs.fill_histogram(pt, 1.0)
            continue

        # the particles in the acceptance of this measurement
        accepted = []
        for (apid, pt, eta, phi, rapidity) in charged:
            if entry["no_strange_baryons"] and apid in STRANGE_BARYONS:
                continue
            if pt < entry["ptmin"]:
                continue
            if abs(eta) >= entry["etamax"]:
                continue
            accepted.append((pt, eta, phi))

        # the number of particles that selects the event
        if entry["sel_etamax"] == entry["etamax"]:
            n_select = len(accepted)
        else:
            n_select = 0
            for (pt, eta, phi) in accepted:
                if abs(eta) < entry["sel_etamax"]:
                    n_select += 1

        if kind == "ue_nch" or kind == "ue_ptsum":
            # underlying event: leading track with pT > 1 GeV; transverse region
            # 60 < |dphi| < 120 degrees to it; densities per unit eta and phi
            if len(accepted) == 0:
                continue
            lead_pt = 0.0
            lead_phi = 0.0
            for (pt, eta, phi) in accepted:
                if pt > lead_pt:
                    lead_pt = pt
                    lead_phi = phi
            if lead_pt < 1.0:
                continue
            obs.selected += 1
            n_transverse = 0
            ptsum_transverse = 0.0
            for (pt, eta, phi) in accepted:
                dphi = delta_phi(phi, lead_phi)
                if dphi >= math.pi / 3 and dphi < 2 * math.pi / 3:
                    n_transverse += 1
                    ptsum_transverse += pt
            area = 2 * entry["etamax"] * 2 * math.pi / 3        # d eta d phi of the transverse region
            if kind == "ue_nch":
                obs.fill_profile(lead_pt, n_transverse / area)
            else:
                obs.fill_profile(lead_pt, ptsum_transverse / area)
            continue

        if n_select < entry["nchmin"]:
            continue
        obs.selected += 1

        if kind == "nch":
            obs.fill_histogram(n_select, 1.0)
        elif kind == "eta":
            for (pt, eta, phi) in accepted:
                obs.fill_histogram(eta, 1.0)
        elif kind == "pt":
            # weight 1/pT: the invariant yield 1/(2 pi pT) d2N/deta dpT (2 pi and eta range in combine.py)
            for (pt, eta, phi) in accepted:
                obs.fill_histogram(pt, 1.0 / pt)
        elif kind == "ptnch":
            # every particle adds its pT to the bin of the event's n_ch
            for (pt, eta, phi) in accepted:
                obs.fill_profile(n_select, pt)


def main(stream, energy, out_file):
    observables = []
    for entry in observables_at(energy):
        observables.append(Observable(entry))

    n_events = 0
    weight_min = None
    weight_max = None
    particles = []
    in_event = False
    with open(stream) as f:
        for line in f:
            tag = line[:2]
            if tag == "E ":
                if in_event:
                    analyse_event(particles, observables)
                    n_events += 1
                in_event = True
                particles = []
            elif tag == "W " and in_event:
                columns = line.split()
                # a line cut off by a stopped Sherpa (watch_jobs.py) is skipped
                if len(columns) < 2:
                    continue
                weight = float(columns[1])
                if weight_min is None or weight < weight_min:
                    weight_min = weight
                if weight_max is None or weight > weight_max:
                    weight_max = weight
            elif tag == "P " and in_event:
                # P id parent pdg px py pz e m status
                columns = line.split()
                # a line cut off by a stopped Sherpa (watch_jobs.py) is skipped
                if len(columns) < 10:
                    continue
                if columns[9] != "1":
                    continue
                pid = int(columns[3])
                px = float(columns[4])
                py = float(columns[5])
                pz = float(columns[6])
                e = float(columns[7])
                particles.append((pid, px, py, pz, e))
    if in_event:
        analyse_event(particles, observables)
        n_events += 1

    result = {}
    result["energy"] = energy
    result["events"] = n_events
    result["weight_min_max"] = [weight_min, weight_max]      # all events must have the same weight (unweighted)
    result["observables"] = {}
    for obs in observables:
        result["observables"][obs.name] = obs.result()
    with open(out_file, "w") as f:
        json.dump(result, f)
    print(n_events, "events ->", out_file)


if __name__ == "__main__":
    main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
