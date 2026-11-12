"""Copy the measured values of every histogram in observables.py from Rivet's reference
files into plain text files, so that they can be read (and checked) without YODA:

    data/<name>.dat     columns: bin low edge, bin high edge, value, total error

The total error is the one Rivet stores (statistical and systematic uncertainties added
in quadrature by HEPData), averaged over its down and up parts.

    python3 make_reference.py
"""
import os
import yoda
from observables import OBSERVABLES

RIVET_DATA = "/Users/ekremorhandemiray/code_software/rivetbootstrap/local/share/Rivet"
from paths import TUNING
OUT = os.path.join(TUNING, "data")


def main():
    os.makedirs(OUT, exist_ok=True)
    for entry in OBSERVABLES:
        reference_file = os.path.join(RIVET_DATA, entry["analysis"] + ".yoda.gz")
        all_histograms = yoda.read(reference_file)
        path = "/REF/" + entry["analysis"] + "/" + entry["histogram"]
        histogram = all_histograms[path]
        edges = histogram.xEdges()
        bins = histogram.bins()
        out_file = os.path.join(OUT, entry["name"] + ".dat")
        with open(out_file, "w") as f:
            f.write("# " + entry["analysis"] + " " + entry["histogram"] + "  sqrt(s) = " + str(entry["energy"]) + " GeV\n")
            f.write("# kind " + entry["kind"] + ", charged particles with pT > " + str(entry["ptmin"]) + " GeV, |eta| < "
                    + str(entry["etamax"]) + ", events with n_ch >= " + str(entry["nchmin"]) + "\n")
            f.write("# low_edge high_edge value total_error\n")
            for i in range(len(bins)):
                low = edges[i]
                high = edges[i + 1]
                value = bins[i].val()
                error = bins[i].totalErrAvg()
                f.write(f"{low:.6g} {high:.6g} {value:.6g} {error:.6g}\n")
        print(entry["name"], len(bins), "bins")


if __name__ == "__main__":
    main()
