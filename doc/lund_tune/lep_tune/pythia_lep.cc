// Pythia 8 at the Z pole, for the LEP plots of the combined report: the same quantities as
// doc/scripts/lep_observables.py fills from Sherpa's events, written in the same lep.json format,
// so that compare_lep.py plots both.
//   n_ch                  charged multiplicity (ALEPH_1991_I319520)
//   x_p of pi, K, p       x_p = |p| / E_beam (ALEPH_1995_I382179)
// e+e- -> gamma*/Z -> q qbar (u d s c b) at 91.2 GeV, no initial-state radiation (PDF:lepton off),
// Pythia's default (Monash) tune.  Particles with c tau > 1000 mm are stable, as in the Sherpa
// LEP runs (HADRON_DECAYS: Max_Proper_Lifetime 1000): K0S and Lambda decay, as ALEPH counts them.
//
//   g++ -O2 -std=c++17 pythia_lep.cc -o pythia_lep \
//     $(/Users/ekremorhandemiray/Downloads/pythia8318/bin/pythia8-config --cxxflags --libs) \
//     -Wl,-rpath,/Users/ekremorhandemiray/Downloads/pythia8318/lib
//   ./pythia_lep <events> <seed> <output lep.json>

#include "Pythia8/Pythia.h"
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <vector>

using namespace Pythia8;

int main(int argc, char* argv[]) {
  if (argc != 4) {
    std::cerr << "usage: pythia_lep <events> <seed> <output lep.json>" << std::endl;
    return 1;
  }
  int n_events = std::stoi(argv[1]);
  int seed = std::stoi(argv[2]);
  std::string out_name = argv[3];
  double e_beam = 45.6;

  Pythia pythia;
  pythia.readString("Beams:idA = 11");
  pythia.readString("Beams:idB = -11");
  pythia.readString("Beams:eCM = 91.2");
  pythia.readString("PDF:lepton = off");
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
  pythia.readString("23:onMode = off");
  pythia.readString("23:onIfAny = 1 2 3 4 5");
  pythia.readString("ParticleDecays:limitTau0 = on");
  pythia.readString("ParticleDecays:tau0Max = 1000");
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 100000");
  pythia.init();

  std::map<int, long> nch_counts;
  std::vector<double> xp_pi;
  std::vector<double> xp_K;
  std::vector<double> xp_p;
  long accepted = 0;
  for (int i = 0; i < n_events; ++i) {
    if (!pythia.next()) continue;
    accepted += 1;
    int nch = 0;
    for (int j = 0; j < pythia.event.size(); ++j) {
      const Particle& particle = pythia.event[j];
      if (!particle.isFinal()) continue;
      if (!particle.isCharged()) continue;
      nch += 1;
      int id = std::abs(particle.id());
      // x_p = |p| / E_beam
      double xp = particle.pAbs() / e_beam;
      if (id == 211) xp_pi.push_back(xp);
      if (id == 321) xp_K.push_back(xp);
      if (id == 2212) xp_p.push_back(xp);
    }
    nch_counts[nch] += 1;
  }
  pythia.stat();

  // lep.json: {"events": N, "ebeam": E, "nch": {"<n>": count}, "xp": {"pi": [...], "K": [...], "p": [...]}}
  std::ofstream out(out_name);
  out << "{\"events\": " << accepted << ", \"ebeam\": " << e_beam << ", \"nch\": {";
  bool first = true;
  for (std::map<int, long>::iterator it = nch_counts.begin(); it != nch_counts.end(); ++it) {
    if (!first) out << ", ";
    out << "\"" << it->first << "\": " << it->second;
    first = false;
  }
  out << "}, \"xp\": {";
  const char* names[3] = {"pi", "K", "p"};
  std::vector<double>* lists[3] = {&xp_pi, &xp_K, &xp_p};
  for (int s = 0; s < 3; ++s) {
    if (s > 0) out << ", ";
    out << "\"" << names[s] << "\": [";
    for (size_t k = 0; k < lists[s]->size(); ++k) {
      if (k > 0) out << ", ";
      char buffer[32];
      std::snprintf(buffer, sizeof(buffer), "%.5f", (*lists[s])[k]);
      out << buffer;
    }
    out << "]";
  }
  out << "}}" << std::endl;
  std::cout << accepted << " events -> " << out_name << std::endl;
  return 0;
}
