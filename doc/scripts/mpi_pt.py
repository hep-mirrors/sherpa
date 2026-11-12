"""The pT of every parton scattering of the multiple interactions, per event: Sherpa against
Pythia, at parton level.  This separates the two possible sources of too many charged
particles at high pT: too many hard scatterings, or too hard jets from each scattering.

    python3 mpi_pt.py sherpa <HepMC3 file or pipe> <output .dat>
    python3 mpi_pt.py compare <sherpa .dat> <pythia .dat>

Sherpa writes each scattering (a Hard_Collision blob, the hardest one included in
minimum-bias events) as a HepMC vertex with status 2 (SHERPA/Tools/HepMC3_Interface.C:704);
the first two particles out of it are the two scattered partons, with equal and opposite pT.
The Pythia file comes from ../minimum_bias_tuning/pythia_mpi_pt.cc, which takes the pT from
pythia.info.pTMPI(i).  Both files have the same columns:
    low edge, high edge, scatterings per event (all), scatterings per event (hardest only)
"""
import math
import sys

N_BINS = 60
WIDTH = 0.5


def read_sherpa(stream, out_file):
    all_scatterings = [0.0] * N_BINS
    hardest = [0.0] * N_BINS
    n_events = 0
    n_scatterings = 0
    event_pts = []            # pT of the scatterings of the current event
    mpi_vertex = None         # id of the current status-2 vertex, until its first particle is read
    in_event = False

    def finish_event():
        "fill the histograms with the scatterings of one event; the hardest is the largest pT"
        largest = None
        for pt in event_pts:
            b = int(pt / WIDTH)
            if b < N_BINS:
                all_scatterings[b] += 1.0
            if largest is None or pt > largest:
                largest = pt
        if largest is not None:
            b = int(largest / WIDTH)
            if b < N_BINS:
                hardest[b] += 1.0

    with open(stream) as f:
        for line in f:
            tag = line[:2]
            if tag == "E ":
                if in_event:
                    finish_event()
                    n_events += 1
                    n_scatterings += len(event_pts)
                in_event = True
                event_pts = []
                mpi_vertex = None
            elif tag == "V ":
                # V <id> <status> [parents] @ position
                columns = line.split()
                if columns[2] == "2":
                    mpi_vertex = columns[1]
                else:
                    mpi_vertex = None
            elif tag == "P " and mpi_vertex is not None:
                # P <id> <production vertex> <pdg> px py pz e m status
                columns = line.split()
                if columns[2] != mpi_vertex:
                    continue
                px = float(columns[4])
                py = float(columns[5])
                # pT = sqrt(px^2 + py^2) of the first scattered parton; the second has the same pT
                event_pts.append(math.sqrt(px * px + py * py))
                mpi_vertex = None
    if in_event:
        finish_event()
        n_events += 1
        n_scatterings += len(event_pts)

    with open(out_file, "w") as out:
        out.write("# Sherpa, scatterings = HepMC vertices with status 2\n")
        out.write("# events " + str(n_events) + "\n")
        out.write("# mean scatterings per event " + str(n_scatterings / n_events) + "\n")
        out.write("# low high all_per_event hardest_per_event\n")
        for b in range(N_BINS):
            out.write("%g %g %g %g\n" % (b * WIDTH, (b + 1) * WIDTH,
                                          all_scatterings[b] / n_events, hardest[b] / n_events))
    print(n_events, "events,", n_scatterings / n_events, "scatterings per event ->", out_file)


def read_dat(name):
    "header lines, and the columns low, high, all, hardest"
    header = []
    lows = []
    highs = []
    all_values = []
    hardest_values = []
    for line in open(name):
        if line.startswith("#"):
            header.append(line.strip())
            continue
        columns = line.split()
        lows.append(float(columns[0]))
        highs.append(float(columns[1]))
        all_values.append(float(columns[2]))
        hardest_values.append(float(columns[3]))
    return header, lows, highs, all_values, hardest_values


def compare(sherpa_file, pythia_file):
    sherpa = read_dat(sherpa_file)
    pythia = read_dat(pythia_file)
    print("Sherpa:", sherpa[0][1], "|", sherpa[0][2])
    print("Pythia:", pythia[0][1], "|", pythia[0][2])
    print("scatterings per event in each pT bin (all scatterings), and Sherpa/Pythia")
    print("%-11s %12s %12s %8s   %12s %12s %8s" % ("pT [GeV]", "Sherpa all", "Pythia all", "ratio",
                                                 "Sherpa hard", "Pythia hard", "ratio"))
    for b in range(len(sherpa[1])):
        low = sherpa[1][b]
        if low not in (0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0):
            continue
        row = "%4.1f-%4.1f " % (low, sherpa[2][b])
        for column in (3, 4):
            s = sherpa[column][b]
            p = pythia[column][b]
            if p > 0:
                ratio = "%8.2f" % (s / p)
            else:
                ratio = "%8s" % "-"
            row += " %12.3e %12.3e %s  " % (s, p, ratio)
        print(row)


if __name__ == "__main__":
    if sys.argv[1] == "sherpa":
        read_sherpa(sys.argv[2], sys.argv[3])
    elif sys.argv[1] == "compare":
        compare(sys.argv[2], sys.argv[3])
