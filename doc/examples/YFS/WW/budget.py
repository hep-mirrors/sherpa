#!/usr/bin/env python3
"""Correction budget for the FCC-ee threshold YFS WW study.

Reads the yodas produced by the cards in cards/generated and prints, per energy
and channel, the cross section at each rung of the YFS ladder, the shift each
rung is worth, and that shift expressed as an equivalent shift in m_W.

The ladder is six rungs in four runs. L2 and L3 are not separate productions:
they are exactly-correlated named weights on the L4 sample, so their ratios
carry almost no MC error of their own.

    L0  Born                       L0_born run
    L1  + ISR                      L1_isr run
    L2  + FSR                      L5_lo x NoCoulomb x NoIFI
    L3  + IFI                      L5_lo x NoCoulomb
    L4  + Coulomb                  L5_lo nominal   (= the LO baseline)

and then two ALTERNATIVE higher-order treatments of the same content, each
measured against that same LO baseline rather than against each other:

    EEX (leading-log)              L4_master nominal / L5_lo nominal
    NLO (exact O(alpha))           L6 nominal / L6 YFS.LO

The ladder proper runs at BETA: 0 throughout. At the default BETA: 2 a
Born-level run puts the EEX series into its own nominal -- nothing overwrites it
there -- so mixing BETA between rungs silently folds a +4% leading-log
correction into the comparison. In an NLO run BETA cannot reach the nominal at
all (the NLO block overwrites m_real), so L6 keeps the default and gets a valid
YFS.EEX column out of it.

dsigma/dm_W comes from the generator itself: three ISR-only runs at
m_W = 80.279 / 80.379 / 80.479, central-differenced. Nothing external enters the
conversion.

Usage:  budget.py <rundir> [<rundir> ...]
where each rundir holds the per-rung subdirectories (L0_born/, L1_isr/, ...)
and optionally an m_W scan in mw/<mass>/.
"""
import gzip, glob, math, os, sys


def combine(parts):
    """Inverse-variance mean of independent runs, per column.

    Seeds are combined here rather than with yodamerge because the budget only
    needs the cross sections, and a 1/err^2 mean is the right combination for
    independent samples of the same integral.
    """
    out = {}
    keys = set()
    for p in parts:
        keys |= set(p)
    for k in keys:
        num = den = 0.
        for p in parts:
            if k not in p:
                continue
            v, e = p[k]
            if e <= 0.:
                continue
            w = 1. / (e * e)
            num += w * v
            den += w
        if den > 0.:
            out[k] = (num / den, den ** -0.5)
    return out


def xsecs(path):
    """{column -> (value, error)} from a Sherpa+Rivet yoda. '' is the nominal.

    If path holds no yoda itself but its subdirectories do, they are treated as
    independent seeds and combined.
    """
    direct = _xsecs_one(path)
    if direct:
        return direct
    parts = [_xsecs_one(d) for d in sorted(glob.glob(os.path.join(path, "*")))
             if os.path.isdir(d)]
    parts = [p for p in parts if p]
    return combine(parts) if parts else {}


def _xsecs_one(path):
    out, cur, grab = {}, None, False
    for fn in sorted(glob.glob(os.path.join(path, "*.yoda.gz"))):
        with gzip.open(fn, "rt") as f:
            for line in f:
                if line.startswith("BEGIN") and "/_XSEC" in line and "/RAW/" not in line:
                    tag = line.strip().split()[-1]
                    cur = tag.split("[")[1][:-1] if "[" in tag else ""
                    grab = False
                    continue
                if cur is not None and line.startswith("# value"):
                    grab = True
                    continue
                if grab:
                    p = line.split()
                    out[cur] = (float(p[0]), abs(float(p[1])))
                    cur, grab = None, False
        if out:
            break
    return out


def rung_table(rundir):
    """The six rungs as (label, sigma, error). Missing runs are skipped."""
    rows = []
    def get(sub):
        d = os.path.join(rundir, sub)
        return xsecs(d) if os.path.isdir(d) else {}

    born, isr, master, nlo = get("L0_born"), get("L1_isr"), get("L4_master"), get("L6_nlo")
    lo = get("L5_lo")
    # The ladder is built on the BETA:0 baseline when it is there. Falling back
    # to L4_master keeps old run directories readable, but that mixes BETA
    # between rungs -- so say so rather than printing a table that looks fine.
    base = lo if lo else master
    if not lo and master:
        print("  WARNING: no L5_lo (BETA:0) baseline; using L4_master, which "
              "carries the EEX series in its nominal.")

    if born:   rows.append(("L0 Born",            born[""][0],   born[""][1]))
    if isr:    rows.append(("L1 + ISR",           isr[""][0],    isr[""][1]))
    if base:
        n, e = base[""]
        # The named weights are ratios on the same events, so their statistical
        # error on the RATIO is far below the sample's own error. The absolute
        # errors below are the sample error scaled by the ratio -- honest for the
        # cross section, pessimistic for the rung-to-rung shift.
        nocoul = base.get("YFS.NoCoulomb", (n, e))[0] / n
        noifi  = base.get("YFS.NoIFI",     (n, e))[0] / n
        rows.append(("L2 + FSR",     n * nocoul * noifi, e * nocoul * noifi))
        rows.append(("L3 + IFI",     n * nocoul,         e * nocoul))
        rows.append(("L4 + Coulomb", n,                  e))
    return rows


def higher_order(rundir):
    """The two alternative higher-order treatments, each against its own LO.

    Returned as (label, fractional shift, error) so the caller can print them
    apart from the ladder -- they are alternatives, and stacking them as
    successive rows is exactly the mistake that makes the EEX look like part of
    the NLO rung.
    """
    def get(sub):
        d = os.path.join(rundir, sub)
        return xsecs(d) if os.path.isdir(d) else {}
    out = []
    lo, master, nlo = get("L5_lo"), get("L4_master"), get("L6_nlo")
    if lo and master:
        a, ae = lo[""]
        b, be = master[""]
        r = b / a - 1.
        out.append(("EEX (leading-log)", r,
                    abs(b / a) * math.hypot(ae / a, be / b)))
    if nlo and "YFS.LO" in nlo:
        b, be = nlo[""]
        a, ae = nlo["YFS.LO"]
        # Same events in numerator and denominator, so this ratio is far better
        # determined than either cross section; the quoted error is an upper
        # bound, not the real one.
        out.append(("NLO (exact O(alpha))", b / a - 1.,
                    abs(b / a) * math.hypot(ae / a, be / b)))
    return out


def dlnsigma_dmw(rundir):
    """(dln(sigma)/dm_W [1/GeV], its error) from the three-point m_W scan."""
    pts = []
    for d in sorted(glob.glob(os.path.join(rundir, "mw", "*"))):
        v = xsecs(d)
        if not v:
            continue
        mw = float(os.path.basename(d).replace("p", "."))
        pts.append((mw, v[""][0], v[""][1]))
    if len(pts) < 3:
        return None, None
    pts.sort()
    (m0, s0, e0), (_, s1, _), (m2, s2, e2) = pts[0], pts[1], pts[2]
    d = (s2 - s0) / (m2 - m0) / s1
    # central difference: the two endpoint errors add in quadrature
    err = math.hypot(e0, e2) / (m2 - m0) / s1
    return d, err


def report(rundir):
    print("=" * 78)
    print(rundir)
    rows = rung_table(rundir)
    if not rows:
        print("  no runs found")
        return
    slope, slope_err = dlnsigma_dmw(rundir)
    if slope:
        print("  dln(sigma)/dm_W = %+.4f /GeV  (+-%.4f)   "
              "=> 1%% in sigma == %.1f MeV in m_W"
              % (slope, slope_err, 1000 * 0.01 / abs(slope)))
    else:
        print("  dln(sigma)/dm_W unavailable (need the three-point m_W scan)")

    hdr = "%-16s %12s %10s %12s" % ("rung", "sigma [pb]", "stat", "d(sigma)/sigma")
    if slope:
        hdr += " %12s" % "= dm_W [MeV]"
    print("  " + hdr)
    prev = None
    for label, s, e in rows:
        # A rung whose own error is comparable to the shift it is supposed to
        # measure is not a measurement. Marked rather than dropped, so an
        # unconverged run is visible instead of quietly absent.
        bad = " (!)" if 100 * e / s > 5. else ""
        line = "%-12s%-4s %12.6g %9.2f%% " % (label, bad, s, 100 * e / s)
        if prev is None:
            line += "%12s" % "-"
            if slope:
                line += " %12s" % "-"
        else:
            rel = s / prev - 1.
            line += "%+11.3f%%" % (100 * rel)
            if slope:
                line += " %+12.1f" % (1000 * rel / slope)
        print("  " + line)
        prev = s

    ho = higher_order(rundir)
    if ho:
        print("  " + "-" * 60)
        print("  higher-order treatments (ALTERNATIVES, each vs its own LO):")
        for label, r, err in ho:
            line = "%-22s %+11.3f%%" % (label, 100 * r)
            if slope:
                line += " %12s" % ("%+.1f MeV" % (1000 * r / slope))
            print("  " + line)
    if any(100 * e / sv > 5. for _, sv, e in rows):
        print("  (!) statistically unconverged - not a result")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for d in sys.argv[1:]:
        report(d)
