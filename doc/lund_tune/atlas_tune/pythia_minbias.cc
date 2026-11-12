// Pythia 8 (default Monash 2013 tune) minimum-bias pp events, written as HepMC3 for
// mb_analysis.py.  Pythia is a reference here: its minimum-bias model was tuned to
// LHC data, and it has a realistic diffraction model, so it shows how much the
// measurements with n_ch >= 1 or 2 depend on diffraction.
//
// Usage: pythia_minbias <sqrt(s) [GeV]> <events> <seed> <inelastic|nondiffractive> <hepmc output>
//   inelastic      : SoftQCD:inelastic (non-diffractive + single + double + central diffraction)
//   nondiffractive : SoftQCD:nonDiffractive only (what Sherpa's AMISIC Perturbative mode makes)
// Particles with c tau > 10 mm are not decayed (K0S, Lambda, Sigma, Xi, Omega), as in Sherpa
// and in the ATLAS/ALICE definition of primary particles.
//
// Build:
//   g++ -O2 -std=c++17 pythia_minbias.cc -o pythia_minbias \
//     $(/Users/ekremorhandemiray/Downloads/pythia8318/bin/pythia8-config --cxxflags --libs) \
//     -I/Users/ekremorhandemiray/code_software/rivetbootstrap/local/include \
//     -L/Users/ekremorhandemiray/code_software/rivetbootstrap/local/lib -lHepMC3 \
//     -Wl,-rpath,/Users/ekremorhandemiray/code_software/rivetbootstrap/local/lib \
//     -Wl,-rpath,/Users/ekremorhandemiray/Downloads/pythia8318/lib

#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC3.h"
#include <iostream>
#include <string>

using namespace Pythia8;

int main(int argc, char* argv[]) {
  if (argc != 6) {
    std::cerr << "usage: pythia_minbias <sqrt(s) [GeV]> <events> <seed> <inelastic|nondiffractive> <hepmc output>\n";
    return 1;
  }
  double sqrts = std::stod(argv[1]);
  long nevents = std::stol(argv[2]);
  int seed = std::stoi(argv[3]);
  std::string mode = argv[4];
  std::string output = argv[5];

  Pythia pythia;
  pythia.readString("Beams:idA = 2212");
  pythia.readString("Beams:idB = 2212");
  pythia.readString("Beams:eCM = " + std::to_string(sqrts));
  if (mode == "inelastic") {
    pythia.readString("SoftQCD:inelastic = on");
  } else if (mode == "nondiffractive") {
    pythia.readString("SoftQCD:nonDiffractive = on");
  } else {
    std::cerr << "unknown mode " << mode << "\n";
    return 1;
  }
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 10");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 100000");
  if (!pythia.init()) return 1;

  Pythia8ToHepMC hepmc(output);
  for (long i = 0; i < nevents; ++i) {
    if (!pythia.next()) continue;
    hepmc.writeNextEvent(pythia);
  }
  pythia.stat();
  return 0;
}
