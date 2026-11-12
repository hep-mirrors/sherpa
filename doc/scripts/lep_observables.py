"""LEP observables from a HepMC3 (ASCII) event stream, as the Rivet analyses ALEPH_1991_I319520
and ALEPH_1995_I382179 define them.  Used instead of Rivet because the local Rivet 4.1.4 / YODA
installation crashes when it writes its results (and Sherpa's Rivet interface crashes on start).

    mkfifo events.fifo
    python3 lep_observables.py events.fifo out.json &
    Sherpa -f card.yaml 'EVENT_OUTPUT: HepMC3_GenEvent[events.fifo]'

Per event: all final-state (status 1) charged particles.  Filled:
  n_ch                     charged multiplicity (ALEPH_1991_I319520, d01 and the mean d02)
  xp_pi, xp_K, xp_p        x_p = |p| / E_beam of pi+-, K+-, p/pbar (ALEPH_1995_I382179, d01-d03)
The output json holds the raw sums; compare_lep.py normalises them like Rivet (1/N_events,
per unit x_p) and compares with the reference data."""
import json, sys
from collections import Counter, defaultdict

QUARK_CHARGE3 = {1: -1, 2: 2, 3: -1, 4: 2, 5: -1, 6: 2}     # 3 x charge of d u s c b t

def charge3(pid):
    "three times the electric charge of a PDG code (leptons, gauge bosons, hadrons)"
    a = abs(pid); sign = 1 if pid > 0 else -1
    if a in (11, 13, 15): return -3*sign
    if a in (12, 14, 16, 22, 21, 23, 25, 130, 310): return 0
    if a == 24: return 3*sign
    if a > 1000000000: return 0                              # nuclei (deuterons are not formed here)
    n = a % 10000
    q1, q2, q3 = (n // 1000) % 10, (n // 100) % 10, (n // 10) % 10
    if q1 == 0:                                              # meson q2 q3bar
        if q2 == 0: return 0
        # PDG convention: an even (up-type) q2 is the quark, an odd one the antiquark
        # (211 = u dbar, 321 = u sbar, 411 = c dbar, 521 = u bbar)
        c = QUARK_CHARGE3[q2] - QUARK_CHARGE3[q3]
        return c*sign if q2 % 2 == 0 else -c*sign
    return (QUARK_CHARGE3[q1] + QUARK_CHARGE3[q2] + QUARK_CHARGE3[q3])*sign

# the common charged final-state particles, checked explicitly
for pid, c in ((211, 3), (321, 3), (411, 3), (431, 3), (521, 3), (311, 0), (111, 0), (2212, 3), (3222, 3), (3112, -3), (3312, -3), (3334, -3), (11, -3), (13, -3)):
    assert charge3(pid) == c and charge3(-pid) == -c, pid

def main(stream, out, ebeam=45.6):
    nevents, nch_counts = 0, Counter()
    xp = {"pi": [], "K": [], "p": []}
    species = {211: "pi", 321: "K", 2212: "p"}
    nch = 0; in_event = False
    with open(stream) as f:
        for line in f:
            tag = line[:1]
            if tag == "E":
                if in_event: nch_counts[nch] += 1; nevents += 1
                in_event, nch = True, 0
            elif tag == "P" and in_event:
                parts = line.split()
                # P id parent pdg px py pz e m status
                if parts[9] != "1": continue
                pid = int(parts[3])
                if charge3(pid) == 0: continue
                nch += 1
                sp = species.get(abs(pid))
                if sp:
                    px, py, pz = float(parts[4]), float(parts[5]), float(parts[6])
                    xp[sp].append(round((px*px + py*py + pz*pz)**0.5/ebeam, 5))
    if in_event: nch_counts[nch] += 1; nevents += 1
    json.dump({"events": nevents, "ebeam": ebeam, "nch": dict(nch_counts), "xp": xp}, open(out, "w"))
    print(f"{nevents} events -> {out}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])
