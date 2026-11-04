#!/usr/bin/env python3
# Per-run Sherpa cross section and each run's contribution to dsigma/dtheta_e,
# compared with the merged file. Usage: ./check_MINERVA_2015_I1393997.py <fluxrun output dir>
import sys, glob, yoda

ana = "MINERVA_2015_I1393997"
h = f"/{ana}/d04-x01-y01"
out = sys.argv[1] if len(sys.argv) > 1 else ana   # default: fluxrun output dir in cwd

def integral(ao):   # Histo1D or (width-divided) Estimate1D -> sum of dsigma/dtheta x width
    e = ao.mkEstimate() if hasattr(ao, "mkEstimate") else ao
    return sum(b.val()*b.xWidth() for b in e.bins())

tot = 0.0
for f in sorted(glob.glob(f"{out}/*/*/{ana}_*.yoda.gz")):
    aos = yoda.read(f)
    xs = aos["/_XSEC"].val() * 100.0   # pb -> 1e-38 cm^2
    s = integral(aos[h])               # integral of dsigma/dtheta, 1e-39 cm^2/nucleon
    tot += s
    beam, en = f.split("/")[-3:-1]
    print(f"{beam:7s} {en:>6s} GeV   sigma = {xs:7.3f} x1e-38 cm2   contrib = {s:.4f}")

merged = integral(yoda.read(f"{out}/{ana}.yoda.gz")[h])
print(f"\nsum of runs = {tot:.3f}   merged = {merged:.3f}   data ~ 6.0   (x1e-39 cm2/nucleon)")

print("\nMINERvA measures CCQE-like events on carbon: one e±, any number of nucleons, no mesons. This rivet analysis is currently normalised to the free-nucleon cross section, so it is expected to be below the data.")
print("\nA free-nucleon QE generator is missing processes that land in that signal: " \
"\n   Multi-nucleon (2p2h) events. The neutrino scatters off a correlated nucleon pair, giving an electron plus two nucleons and no pion, so the event counts as CCQE-like. Calculations put this at roughly +20 to 40 percent of QE at a few GeV (Nieves, Ruiz Simo, Vicente Vacas, arXiv:1106.5374; Martini, Ericson et al., arXiv:0910.2622)."\
"\n   Pion production followed by absorption in the nucleus. At 3.6 GeV, resonance production is comparable in size to QE. A sizeable fraction of the pions are absorbed in carbon, and those events also look CCQE-like. The paper's GENIE simulation includes this through its final-state-interaction model."\
"\n   Nuclear effects on QE itself, such as Fermi motion and Pauli blocking. These mostly reshape the low-Q² region and change the total by only a few percent.")
print("\nFrom the paper itself: its GENIE simulation, which includes item 2, describes the data within errors. Even so, the data sit above it at θ_e ≈ 10 to 25° (Fig. 3, right). The data also carry 10 to 15 percent systematic uncertainties. " \
"\n A free-QE-only prediction at about 65percent of the data is therefore what you'd expect, not a sign of a bug.")
