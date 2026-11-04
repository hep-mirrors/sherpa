#!/usr/bin/env python3
# Overlays the nu_e, nu_ebar and total flux written by MINERVA_2015_I1393997.
# Usage: ./plot_flux_MINERVA_2015_I1393997.py [merged yoda, default MINERVA_2015_I1393997/MINERVA_2015_I1393997.yoda.gz]
import sys, yoda
import matplotlib.pyplot as plt

ana = "MINERVA_2015_I1393997"
f = sys.argv[1] if len(sys.argv) > 1 else f"{ana}/{ana}.yoda.gz"
aos = yoda.read(f)

for suffix, label in [("", r"$\nu_e + \bar{\nu}_e$"), ("_nue", r"$\nu_e$"), ("_nuebar", r"$\bar{\nu}_e$")]:
    e = aos[f"/{ana}/d01-x01-y01{suffix}"]
    plt.stairs([b.val() for b in e.bins()], e.xEdges(), label=label)

plt.xlabel(r"$E_\nu$ [GeV]")
plt.ylabel(r"Flux [$10^{-10}$ neutrinos/cm$^2$/GeV/P.O.T.]")
plt.legend()
plt.savefig("flux_MINERVA_2015_I1393997.pdf")
