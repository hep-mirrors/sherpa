#!/usr/bin/env python3
"""Assert a Sherpa+Rivet yoda from this campaign is what it claims to be.

Two things are checked, both of which fail silently otherwise:

1. REQUIRED WEIGHT COLUMNS. The budget reads L2 and L3 off the L4 sample as
   YFS.NoIFI and YFS.NoCoulomb. If Sherpa was built without YFS: Ladder_Weights
   it accepts the setting, drops it, emits no such columns, and budget.py falls
   back to the nominal -- reporting L2 = L3 = L4, a perfectly plausible-looking
   table in which three rungs are the same number.

2. ACCEPTANCE = integral(/yfs-ww/xs) / value(/_XSEC). For every card in this
   campaign the final state has exactly one negative and one positive charged
   lepton (the tau is Stable: 1), so yfsww::analyze vetoes nothing and this
   ratio is 1 by construction. Anything else means the sample is not what the
   card says -- and in a MERGED file it is the signature of rivet-merge having
   dropped objects whose binning did not match across inputs while still
   normalising by the total sumW: the survivor comes out low by exactly the
   fraction of seeds that contributed the other binning. That is the failure
   that once produced a histogram at 0.5 of its own xs.

   The binning DOES change with sqrt(s) here (yfsww books against sqrtS()), so
   merging two energies is a live hazard, not a hypothetical one.

Usage:
    check-yoda.py [--require "YFS.A YFS.B"] [--acceptance-tol 0.005]
                  [--no-acceptance] <file.yoda.gz> [...]
"""
import argparse, gzip, io, os, sys


def _open(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt")
    return io.open(path, "rt")


def read(path):
    """-> (columns, xsec{col:(v,e)}, xs_integral{col:value})

    The xs "histogram" is a single bin over [0,1] filled once per surviving
    event, so after finalize() its integral is crossSection() x acceptance.
    YODA v3 Estimate1D writes an underflow row FIRST and an overflow row LAST
    around the N in-range rows -- both `nan` here -- so the in-range rows are
    rows[1:-1] and nothing else. Reading row 0 would give nan and turn every
    check into a false alarm.
    """
    xsec, xsint = {}, {}
    cur = kind = None
    edges, rows, grab = [], [], False
    with _open(path) as f:
        for line in f:
            if line.startswith("BEGIN "):
                tag = line.split()[-1]
                col = tag.split("[")[1][:-1] if "[" in tag else ""
                base = tag.split("[")[0]
                edges, rows, grab = [], [], False
                if base == "/_XSEC":
                    cur, kind = col, "xsec"
                elif base == "/yfs-ww/xs":
                    cur, kind = col, "xs"
                else:
                    cur = kind = None
                continue
            if line.startswith("END "):
                if kind == "xs" and rows and len(edges) >= 2:
                    inner = rows[1:-1]
                    tot = 0.0
                    for i, v in enumerate(inner):
                        if i + 1 < len(edges):
                            tot += v * (edges[i + 1] - edges[i])
                    xsint[cur] = tot
                cur = kind = None
                grab = False
                continue
            if cur is None:
                continue
            if line.startswith("Edges("):
                edges = [float(x) for x in
                         line[line.index("[") + 1:line.rindex("]")].split(",")]
                continue
            if line.startswith("# value"):
                grab = True
                continue
            if not grab:
                continue
            p = line.split()
            if not p:
                continue
            if kind == "xsec":
                xsec[cur] = (float(p[0]), abs(float(p[1])))
                cur = kind = None
                grab = False
            else:
                try:
                    rows.append(float(p[0]))
                except ValueError:
                    rows.append(float("nan"))
    return sorted(xsec), xsec, xsint


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+")
    ap.add_argument("--require", default="",
                    help="space-separated weight column names that must exist")
    ap.add_argument("--acceptance-tol", type=float, default=0.005)
    ap.add_argument("--no-acceptance", action="store_true")
    a = ap.parse_args()

    need = [c for c in a.require.split() if c]
    rc = 0
    for path in a.files:
        if not os.path.exists(path):
            print("ERROR: no such file: %s" % path, file=sys.stderr)
            rc = 1
            continue
        cols, xsec, xshist = read(path)
        if "" not in xsec:
            print("ERROR: %s has no nominal /_XSEC" % path, file=sys.stderr)
            rc = 1
            continue
        v, e = xsec[""]
        print("  %s" % os.path.basename(path))
        print("    sigma = %.6g pb +- %.3g (%.2f%%)"
              % (v, e, 100 * e / v if v else float("nan")))
        named = [c for c in cols if c and not c.startswith("EXTRA__")]
        print("    columns (%d): %s" % (len(named), " ".join(named) or "<none>"))

        missing = [c for c in need if c not in xsec]
        if missing:
            print("ERROR: %s is missing required columns: %s"
                  % (path, " ".join(missing)), file=sys.stderr)
            print("       Sherpa dropped the YFS setting that produces them, or the"
                  " binary predates it.", file=sys.stderr)
            rc = 1

        if not a.no_acceptance:
            if "" not in xshist:
                print("ERROR: %s has no /yfs-ww/xs -- the analysis did not run"
                      % path, file=sys.stderr)
                rc = 1
            else:
                acc = xshist[""] / v if v else float("nan")
                flag = "" if abs(acc - 1.0) <= a.acceptance_tol else "   <-- WRONG"
                print("    acceptance = xs-histogram / _XSEC = %.6f%s" % (acc, flag))
                if flag:
                    print("ERROR: %s: acceptance %.6f is not 1. Either the analysis"
                          " vetoed events it should not have, or a merge dropped"
                          " objects while keeping the full sumW normalisation."
                          % (path, acc), file=sys.stderr)
                    rc = 1
    return rc


if __name__ == "__main__":
    sys.exit(main())
