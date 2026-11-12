// Pythia 8 reference for the Sherpa Deuteron_Formation phase.
//
// Non-diffractive pp events (or e+e- -> gamma*/Z -> q qbar), particles with
// lifetimes above 100 fm/c kept
// stable (arXiv:1504.07242, Sec. IV), and Pythia's own implementation of the
// Dal & Raklev cross-section model (HadronLevel:DeuteronProduction).  Writes
// the same histograms as SHERPA/Single_Events/Deuteron_Formation.C, so that
// XSecs/Plotting.ipynb can read both:
//   nucleon_{p,pbar}_invyield_midy  1/(2 pi pT) d^2N/dy dpT, |y|<0.5, per event
//   pairs_k_{pn,pp,nn}              (anti)nucleon pairs per event per GeV of k
//   formed_{d,dbar}_invyield_midy   deuterons formed by Pythia's model
//   formed_{d,dbar}_p               deuterons per event per GeV of lab momentum
//   summary.dat                     per-event counts, including the ALEPH and
//                                   OPAL momentum windows (as Deuteron_Formation.C)
//   partner_*                       origin of the deuteron-forming pairs, as
//                                   FillPartners in Deuteron_Formation.C
//   string_pairs_*                  Pythia only: for the pn pairs within k < 0.5 GeV,
//                                   whether the two nucleons come from the same string
//                                   and how far apart along it (rank distance)
//   partner_*anti*                  antibaryons near the protons, as FillPartners
//   meson_{pi,K}_{plus,minus}_invyield_midy  charged pions and kaons, |y|<0.5, as FillMesons
//   multiplicity_*                  ALICE multiplicity classes (V0M, INEL>0), as
//                                   WriteMultiplicityClasses in Deuteron_Formation.C
//
// Usage: pythia_deuterons <sqrt(s) [GeV]> <events> <seed> <1/sigma_0 [1/b]> <outdir> [pp|ee] [coalescence]
// DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence option added (original: pythia_deuterons.cc.before_coalescence)
//   coalescence (optional): <HardDisk|Exponential>:<p_coal [GeV]>:<b_coal [fm]>, the coalescence
//                 models of Deuteron_Production.C (Deuteron_Coalescence::Probabilities), applied to
//                 all p n (pbar nbar) pairs with the same k and impact parameter b; b_coal <= 0:
//                 momentum condition only.  Pythia's own deuteron model is then switched off, its
//                 hadron production vertices on (Fragmentation:setVertices), and the expected
//                 deuterons (weight P per pair, d gamma kinematics) written as expected_* and
//                 pairs_b_pn, as Deuteron_Formation.C
//   pp (default): non-diffractive pp collisions (SoftQCD:nonDiffractive)
//   ee          : e+ e- -> gamma*/Z -> q qbar (no initial-state photons), as Sherpa's
//                 Sherpa.Deuterons.yaml (Z pole) and Sherpa.Deuterons.BaBar.yaml (10.58 GeV)

#include "Pythia8/Pythia.h"
#include <cmath>
#include <fstream>
#include <map>
#include <string>
#include <sys/stat.h>
#include <vector>

using namespace Pythia8;

struct Histo {
  double lo, hi; std::vector<double> w, w2;
  Histo(int n = 50, double l = 0., double h = 5.) : lo(l), hi(h), w(n+1, 0.), w2(n+1, 0.) {}
  void fill(double x, double v) {
    if (x < lo) return;
    size_t i = std::min(size_t((x-lo)/(hi-lo)*(w.size()-1)), w.size()-1);   // last: overflow
    w[i] += v; w2[i] += v*v;
  }
  // same layout as Sherpa's histogram files: lower bin edge, value, error;
  // values per event and per unit of x, last row overflow
  void write(const std::string & path, double nev) const {
    std::ofstream out(path);
    double bw = (hi-lo)/double(w.size()-1);
    out << "# Pythia " << PYTHIA_VERSION << ": bin low edge, per event per unit, error\n";
    for (size_t i = 0; i < w.size(); ++i)
      out << lo + i*bw << "  " << w[i]/(nev*bw) << "  " << std::sqrt(w2[i])/(nev*bw) << "\n";
  }
};

int main(int argc, char * argv[]) {
  if (argc < 6 || argc > 8) {
    std::cerr << "usage: " << argv[0] << " <sqrt(s) GeV> <events> <seed> <1/sigma_0 [1/b]> <outdir> [pp|ee] [HardDisk|Exponential:p_coal:b_coal]\n";
    return 1;
  }
  std::string beams = argc >= 7 ? argv[6] : "pp";
  // DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence model (optional 7th argument)
  std::string coal_model; double p_coal = 0., b_coal = 0.;
  if (argc == 8) {
    std::string spec = argv[7];
    size_t c1 = spec.find(':'), c2 = spec.rfind(':');
    if (c1 == std::string::npos || c2 == c1) { std::cerr << "coalescence: <HardDisk|Exponential>:<p_coal>:<b_coal>\n"; return 1; }
    coal_model = spec.substr(0, c1);
    p_coal = std::atof(spec.substr(c1+1, c2-c1-1).c_str());
    b_coal = std::atof(spec.substr(c2+1).c_str());
    if (coal_model != "HardDisk" && coal_model != "Exponential") { std::cerr << "unknown coalescence model\n"; return 1; }
  }
  const bool coal = !coal_model.empty();
  if (beams != "pp" && beams != "ee") { std::cerr << "beams must be pp or ee\n"; return 1; }
  double ecm = std::atof(argv[1]), invsig0 = std::atof(argv[4]);
  int nev = std::atoi(argv[2]), seed = std::atoi(argv[3]);
  std::string outdir = argv[5];
  mkdir(outdir.c_str(), 0755);

  Pythia pythia;
  pythia.readString("Beams:eCM = " + std::to_string(ecm));
  if (beams == "pp") {
    pythia.readString("SoftQCD:nonDiffractive = on");   // as Sherpa's Amisic and the paper
  } else {
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("PDF:lepton = off");                 // no initial-state photon radiation
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");   // e+ e- -> gamma*/Z
    pythia.readString("23:onMode = off");                  // ... -> q qbar only (d u s c b)
    pythia.readString("23:onIfAny = 1 2 3 4 5");
  }
  pythia.readString("ParticleDecays:limitTau0 = on"); // tau > 100 fm/c stable (paper, Sec. IV)
  pythia.readString("ParticleDecays:tau0Max = 1e-10");  // mm/c
  pythia.readString(std::string("HadronLevel:DeuteronProduction = ") + (coal ? "off" : "on"));
  // DEUTERON-STUDY (E. Demiray, 2026-11-04): Pythia's deuteron model off and vertices on in coalescence mode
  if (coal) pythia.readString("Fragmentation:setVertices = on");   // hadron production points
  // Pythia's normalisation: 1/sigma_0 = 2.63/b corresponds to norm = 119.6
  pythia.readString("DeuteronProduction:norm = " + std::to_string(119.6*2.63/invsig0));
  pythia.readString("Random:setSeed = on");
  pythia.readString("Random:seed = " + std::to_string(seed));
  pythia.readString("Next:numberCount = 50000");
  if (!pythia.init()) return 1;
  // DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence
  const double m_d = 1.875613;                        // deuteron mass [GeV]

  std::map<std::string, Histo> h;
  for (std::string s : {"nucleon_p", "nucleon_pbar", "formed_d", "formed_dbar"})
    h[s + "_invyield_midy"] = Histo(50, 0., 5.);
  for (std::string s : {"formed_d", "formed_dbar"}) h[s + "_p"] = Histo(100, 0., 5.);
  for (std::string s : {"pn", "pp", "nn"}) h["pairs_k_" + s] = Histo(100, 0., 5.);
  // DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence histograms
  if (coal) {
    for (std::string s : {"expected_d", "expected_dbar"}) {
      h[s + "_invyield_midy"] = Histo(50, 0., 5.);
      h[s + "_p"] = Histo(100, 0., 5.);
    }
    h["pairs_b_pn"] = Histo(100, 0., 100.);           // fm, as Deuteron_Formation.C
  }
  for (std::string o : {"direct", "Delta", "Nstar", "strangebaryon", "other"})
    h["proton_origin_" + o + "_invyield_midy"] = Histo(50, 0., 5.);
  // origin of the deuteron-forming pairs, identical to FillPartners in Deuteron_Formation.C
  for (std::string s : {"partner_protons_pT", "partner_protons_k05_pT", "partner_protons_npart05_pT",
                        "partner_protons_cone_pT", "partner_pairs_cone_pT", "partner_pairs_direct_direct_pT",
                        "partner_pairs_direct_decay_pT", "partner_pairs_decay_decay_pT"})
    h[s] = Histo(25, 0., 5.);
  // string origin of the pn pairs within k < 0.5 GeV, in proton pT: same string with rank
  // distance 1-2 (neighbours or one hadron between), 3-5, > 5; different strings; and the
  // same-string pairs with an antibaryon in between (two separate diquark breaks)
  for (std::string s : {"string_pairs_same_rank1-2_pT", "string_pairs_same_rank3-5_pT", "string_pairs_same_rank6+_pT",
                        "string_pairs_different_pT", "string_pairs_same_antibaryon_between_pT", "string_pairs_other_pT"})
    h[s] = Histo(25, 0., 5.);
  // antibaryons near the protons, identical to FillPartners in Deuteron_Formation.C: in proton
  // pT, the protons and the protons with a partner that have a baryon of opposite baryon number
  // within k < 1 GeV; and for proton pT 1.5-5 GeV the k of the nearest such baryon
  for (std::string s : {"partner_protons_anti1_pT", "partner_pairs_anti1_pT"}) h[s] = Histo(25, 0., 5.);
  for (std::string s : {"partner_protons_kanti_pT1.5-5", "partner_pairs_kanti_pT1.5-5"}) h[s] = Histo(50, 0., 5.);
  // charged pions and kaons at |y|<0.5, as FillMesons in Deuteron_Formation.C (ALICE
  // arXiv:1504.00024). Nothing decays weakly here (tau0Max), so all are counted.
  for (std::string s : {"pi_plus", "pi_minus", "K_plus", "K_minus"})
    h["meson_" + s + "_invyield_midy"] = Histo(100, 0., 10.);
  // ALICE multiplicity classes (arXiv:1902.09290), as Deuteron_Formation.C: INEL>0 events
  // (a charged particle at |eta|<1) in rows of nV0, the number of charged particles in
  // V0-A (2.8 <= eta <= 5.1) and V0-C (-3.7 <= eta <= -1.7); per row the number of events,
  // the summed charged particles at |eta|<0.5, the summed dN/dy (|y|<0.5) of d, dbar, p,
  // pbar and their summed spectra 1/(2 pi pT) dN/dpT dy (and squares) in 0.1 GeV bins
  const int mult_nmax = 600, mult_npt = 50;
  const double mult_binwidth = 0.1;
  const std::vector<std::string> mult_species = {"d", "dbar", "p", "pbar"};
  std::vector<double> mult_events(mult_nmax, 0.), mult_nch_sum(mult_nmax, 0.);
  std::map<std::string, std::vector<double>> mult_count_sum, mult_spectrum_sum, mult_spectrum_square_sum;
  for (const std::string & sp : mult_species) {
    mult_count_sum[sp].assign(mult_nmax, 0.);
    mult_spectrum_sum[sp].assign(mult_nmax*mult_npt, 0.);
    mult_spectrum_square_sum[sp].assign(mult_nmax*mult_npt, 0.);
  }
  std::map<std::string, double> sum, sum2;
  int nacc = 0;

  for (int iev = 0; iev < nev; ++iev) {
    if (!pythia.next()) continue;
    ++nacc;
    const Event & ev = pythia.event;
    std::map<std::string, double> c;
    // nucleons entering the formation: final ones, plus those that Pythia's
    // model has just bound into a deuteron (its mothers)
    std::vector<int> nuc;
    for (int i = 0; i < ev.size(); ++i) {
      int id = ev[i].idAbs();
      if (ev[i].isFinal() && (id == 2212 || id == 2112)) nuc.push_back(i);
      if (ev[i].isFinal() && id == 1000010020) {
        bool anti = ev[i].id() < 0;
        std::string sp = anti ? "dbar" : "d";
        c["formed_" + sp + "_total"] += 1.;
        double p = ev[i].pAbs(), cth = p > 0. ? ev[i].pz()/p : 0.;
        h["formed_" + sp + "_p"].fill(p, 1.);
        // ALEPH (hep-ex/0604023) and OPAL windows, as in Deuteron_Formation.C
        if (p > 0.62 && p < 1.03 && std::abs(cth) < 0.95) c["formed_" + sp + "_ALEPH_window"] += 1.;
        if (p > 0.35 && p < 1.1) c["formed_" + sp + "_OPAL_window"] += 1.;
        if (std::abs(ev[i].y()) < 0.5 && ev[i].pT() > 0.)
          h[anti ? "formed_dbar_invyield_midy" : "formed_d_invyield_midy"].fill(ev[i].pT(), 1./(2.*M_PI*ev[i].pT()));
        for (int m : ev[i].motherList()) {
          int idm = ev[m].idAbs();
          if ((idm == 2212 || idm == 2112) && !ev[m].isFinal()) nuc.push_back(m);
        }
      }
    }
    for (int i : nuc) {
      if (ev[i].idAbs() != 2212) continue;
      std::string sp = ev[i].id() > 0 ? "p" : "pbar";
      c["nucleon_" + sp + "_total"] += 1.;
      double pn = ev[i].pAbs();
      if (pn > 0.31 && pn < 0.515 && std::abs(ev[i].pz()) < 0.95*pn) c["nucleon_" + sp + "_ALEPH_halfwindow"] += 1.;
      if (std::abs(ev[i].y()) < 0.5 && ev[i].pT() > 0.) {
        h["nucleon_" + sp + "_invyield_midy"].fill(ev[i].pT(), 1./(2.*M_PI*ev[i].pT()));
        // origin: immediate parent, as ProtonOrigin in Deuteron_Formation.C
        std::string o = "direct";
        int m = ev[i].mother1();
        if (ev[i].statusAbs() == 91 && m > 0) {
          int kf = ev[m].idAbs();
          int q1 = (kf/1000)%10, q2 = (kf/100)%10, q3 = (kf/10)%10;
          if (!ev[m].isHadron() || q1 == 0) o = "other";            // mesons, leptons, ...
          else if (kf == 1114 || kf == 2114 || kf == 2214 || kf == 2224) o = "Delta";
          else if (q1 > 3 || q2 > 3 || q3 > 3) o = "other";
          else if (q1 == 3 || q2 == 3 || q3 == 3) o = "strangebaryon";
          else o = "Nstar";
          // parents of antiprotons charge-conjugated, so that p and pbar add up
          c["proton_midy_parent_" + pythia.particleData.name(ev[i].id() > 0 ? ev[m].id() : -ev[m].id())] += 1.;
        }
        h["proton_origin_" + o + "_invyield_midy"].fill(ev[i].pT(), 1./(2.*M_PI*ev[i].pT()));
        c["proton_midy_from_" + o] += 1.;
      }
    }
    // origin of the deuteron-forming pairs: for every proton at |y| < 0.5, the same-sign
    // neutrons within k < 0.5 GeV, the charged pT in a cone DeltaR < 0.5 around it, and
    // whether proton and nearest neutron are made directly (not in a hadron decay)
    // primary ancestor of a hadron (status 81-89: made in the string fragmentation itself)
    auto primary = [&ev](int i) {
      for (int n = 0; n < 50 && i > 0 && !(ev[i].statusAbs() >= 81 && ev[i].statusAbs() <= 89); ++n) i = ev[i].mother1();
      return i; };
    auto kcm = [](const Particle & A, const Particle & B) {
      double s = (A.p() + B.p()).m2Calc(), mA2 = A.m()*A.m(), mB2 = B.m()*B.m();
      return std::sqrt(std::max(0., (s-mA2-mB2)*(s-mA2-mB2) - 4.*mA2*mB2)/s); };
    std::vector<int> charged;   // as the nucleons before the formation: deuterons out, bound nucleons in
    for (int i = 0; i < ev.size(); ++i)
      if (ev[i].isFinal() && ev[i].isCharged() && ev[i].idAbs() != 1000010020 && ev[i].pT() > 0.15 && std::abs(ev[i].eta()) < 2.)
        charged.push_back(i);
    for (int i : nuc)
      if (!ev[i].isFinal() && ev[i].idAbs() == 2212 && ev[i].pT() > 0.15 && std::abs(ev[i].eta()) < 2.) charged.push_back(i);
    // baryons as the formation sees them: final-state ones (lifetimes above 100 fm/c kept
    // stable, e.g. Lambda) and the nucleons bound into deuterons
    std::vector<int> baryons;
    for (int i = 0; i < ev.size(); ++i)
      if (ev[i].isFinal() && ev[i].isHadron() && ev[i].particleDataEntry().isBaryon()) baryons.push_back(i);
    for (int i : nuc)
      if (!ev[i].isFinal()) baryons.push_back(i);
    for (int i : nuc) {
      if (ev[i].idAbs() != 2212 || std::abs(ev[i].y()) >= 0.5 || ev[i].pT() <= 0.) continue;
      double pT = ev[i].pT(), knear = 1.e12, npart = 0., cone = 0.; int near = -1;
      for (int j : nuc) {
        if (ev[j].idAbs() != 2112 || (ev[j].id() > 0) != (ev[i].id() > 0)) continue;
        double k = kcm(ev[i], ev[j]);
        if (k < 0.5) npart += 1.;
        if (k < knear) { knear = k; near = j; }
      }
      for (int j : charged) {
        if (j == i) continue;
        double dphi = std::abs(ev[j].phi() - ev[i].phi()); if (dphi > M_PI) dphi = 2.*M_PI - dphi;
        if (std::pow(ev[j].eta() - ev[i].eta(), 2) + dphi*dphi < 0.25) cone += ev[j].pT();
      }
      h["partner_protons_pT"].fill(pT, 1.);
      h["partner_protons_npart05_pT"].fill(pT, npart);
      h["partner_protons_cone_pT"].fill(pT, cone);
      // nearest baryon of opposite baryon number
      double kanti = 1.e12;
      for (int j : baryons)
        if ((ev[j].id() > 0) != (ev[i].id() > 0)) kanti = std::min(kanti, kcm(ev[i], ev[j]));
      bool highpT = pT >= 1.5 && pT < 5.;
      if (kanti < 1.) h["partner_protons_anti1_pT"].fill(pT, 1.);
      if (highpT)     h["partner_protons_kanti_pT1.5-5"].fill(kanti, 1.);
      if (near >= 0 && knear < 0.5) {
        h["partner_protons_k05_pT"].fill(pT, 1.);
        h["partner_pairs_cone_pT"].fill(pT, cone);
        if (kanti < 1.) h["partner_pairs_anti1_pT"].fill(pT, 1.);
        if (highpT)     h["partner_pairs_kanti_pT1.5-5"].fill(kanti, 1.);
        // string origin: primary hadrons of one string share its partons (mother1, mother2)
        // and are stored next to each other, ordered along the string (StringFragmentation.cc)
        int a = primary(i), b = primary(near);
        bool string83 = a > 0 && b > 0 && (ev[a].statusAbs() == 83 || ev[a].statusAbs() == 84)
                                         && (ev[b].statusAbs() == 83 || ev[b].statusAbs() == 84);
        if (!string83) h["string_pairs_other_pT"].fill(pT, 1.);
        else if (ev[a].mother1() != ev[b].mother1() || ev[a].mother2() != ev[b].mother2())
          h["string_pairs_different_pT"].fill(pT, 1.);
        else {
          int lo = std::min(a, b), hi = std::max(a, b), dr = hi - lo;
          h[dr <= 2 ? "string_pairs_same_rank1-2_pT" : (dr <= 5 ? "string_pairs_same_rank3-5_pT" : "string_pairs_same_rank6+_pT")].fill(pT, 1.);
          bool anti = false;
          for (int m = lo + 1; m < hi; ++m)
            if (ev[m].isHadron() && ev[m].particleDataEntry().isBaryon() && (ev[m].id() > 0) != (ev[i].id() > 0)) anti = true;
          if (anti) h["string_pairs_same_antibaryon_between_pT"].fill(pT, 1.);
        }
        bool pdir = ev[i].statusAbs() != 91, ndir = ev[near].statusAbs() != 91;
        h[(pdir && ndir) ? "partner_pairs_direct_direct_pT" : ((pdir || ndir) ? "partner_pairs_direct_decay_pT"
                                                                             : "partner_pairs_decay_decay_pT")].fill(pT, 1.);
      }
    }
    // all unique nucleon-nucleon and antinucleon-antinucleon pairs; coalescence: expected deuterons
    std::vector<std::pair<Vec4, double>> coal_d, coal_dbar;   // DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence
    for (size_t a = 0; a < nuc.size(); ++a) {
      for (size_t b = a+1; b < nuc.size(); ++b) {
        const Particle & A = ev[nuc[a]], & B = ev[nuc[b]];
        if ((A.id() > 0) != (B.id() > 0)) continue;
        double s = (A.p() + B.p()).m2Calc(), mA2 = A.m()*A.m(), mB2 = B.m()*B.m();
        double lam = std::max(0., (s-mA2-mB2)*(s-mA2-mB2) - 4.*mA2*mB2);
        double k = std::sqrt(lam/s);                 // = 2 |p*| = |p1 - p2| in the pair c.m.
        int np = (A.idAbs() == 2212) + (B.idAbs() == 2212);
        std::string t = np == 1 ? "pn" : (np == 2 ? "pp" : "nn");
        h["pairs_k_" + t].fill(k, 1.);
        c["pairs_" + t] += 1.;
        // DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence: b, probability, d gamma kinematics (to the end of this loop)
        if (!coal || np != 1) continue;
        // impact parameter of the straight trajectories in the pair frame, as
        // Deuteron_Coalescence::ImpactParameter2 (vertices in mm)
        Vec4 P = A.p() + B.p(), xA = A.vProd(), xB = B.vProd(), pA = A.p();
        xA.bstback(P); xB.bstback(P); pA.bstback(P);
        Vec4 dx = xA - xB; dx.e(0.);
        Vec4 n = pA; n.e(0.); if (n.pAbs() > 0.) n /= n.pAbs();
        double dl = dx.px()*n.px() + dx.py()*n.py() + dx.pz()*n.pz();
        double b2 = std::max(0., dx.pAbs2() - dl*dl)*1.e24;   // mm^2 -> fm^2
        h["pairs_b_pn"].fill(std::sqrt(b2), 1.);
        // probability, as Deuteron_Coalescence::Probabilities
        bool space = b_coal > 0.;
        double prob;
        if (coal_model == "HardDisk") prob = ((!space || b2 < b_coal*b_coal) && k < p_coal) ? 1. : 0.;
        else {
          double C_b = b_coal*b_coal/(0.75*M_PI), C_p = p_coal*p_coal/(0.75*M_PI);
          prob = 0.75*(space ? std::exp(-b2/C_b) : 1.)*std::exp(-k*k/C_p);
        }
        if (prob <= 0.) continue;
        // p n -> d gamma: isotropic two-body decay of the pair, as FinalState for pn_dgamma
        double M = std::sqrt(s);
        if (M <= m_d) continue;
        double pstar = (M*M - m_d*m_d)/(2.*M), cth = 2.*pythia.rndm.flat() - 1., phi = 2.*M_PI*pythia.rndm.flat();
        double sth = std::sqrt(std::max(0., 1. - cth*cth));
        Vec4 pd(pstar*sth*std::cos(phi), pstar*sth*std::sin(phi), pstar*cth, std::sqrt(pstar*pstar + m_d*m_d));
        pd.bst(P);
        bool anti = A.id() < 0;
        std::string sp = anti ? "dbar" : "d";
        double pabs = pd.pAbs(), cosd = pabs > 0. ? pd.pz()/pabs : 0.;
        h["expected_" + sp + "_p"].fill(pabs, prob);
        if (std::abs(pd.rap()) < 0.5 && pd.pT() > 0.)
          h["expected_" + sp + "_invyield_midy"].fill(pd.pT(), prob/(2.*M_PI*pd.pT()));
        c["expected_" + sp + "_total"] += prob;
        (anti ? coal_dbar : coal_d).push_back({pd, prob});
        if (pabs > 0.62 && pabs < 1.03 && std::abs(cosd) < 0.95) c["expected_" + sp + "_ALEPH_window"] += prob;
        if (pabs > 0.35 && pabs < 1.1) c["expected_" + sp + "_OPAL_window"] += prob;
      }
    }
    // charged pions and kaons at |y|<0.5
    for (int i = 0; i < ev.size(); ++i) {
      int id = ev[i].idAbs();
      if (!ev[i].isFinal() || (id != 211 && id != 321)) continue;
      if (std::abs(ev[i].y()) >= 0.5 || ev[i].pT() <= 0.) continue;
      std::string sp = std::string(id == 211 ? "pi" : "K") + (ev[i].id() > 0 ? "_plus" : "_minus");
      h["meson_" + sp + "_invyield_midy"].fill(ev[i].pT(), 1./(2.*M_PI*ev[i].pT()));
    }
    // event activity, from the charged particles as before the deuteron formation (as
    // CountEventActivity in Deuteron_Formation.C): deuterons out, the protons bound in them in
    int nV0 = 0, nch_mid = 0; bool inel0 = false;
    std::vector<int> activity;
    for (int i = 0; i < ev.size(); ++i)
      if (ev[i].isFinal() && ev[i].isCharged() && ev[i].idAbs() != 1000010020) activity.push_back(i);
    for (int i : nuc)
      if (!ev[i].isFinal() && ev[i].isCharged()) activity.push_back(i);
    for (int i : activity) {
      if (ev[i].pT() <= 0.) continue;                     // along the beam: no eta
      double eta = ev[i].eta();
      if ((eta >= 2.8 && eta <= 5.1) || (eta >= -3.7 && eta <= -1.7)) ++nV0;
      if (std::abs(eta) < 0.5) ++nch_mid;
      if (std::abs(eta) < 1.)  inel0 = true;
    }
    if (inel0) {
      // this event's dN/dy and spectra at |y|<0.5: deuterons formed, (anti)protons entering the formation
      std::map<std::string, double> count;
      std::map<std::string, std::vector<double>> spectrum;
      for (const std::string & sp : mult_species) { count[sp] = 0.; spectrum[sp].assign(mult_npt, 0.); }
      auto add = [&](const Particle & part, const std::string & sp) {
        if (std::abs(part.y()) >= 0.5 || part.pT() <= 0.) return;
        int ipt = int(part.pT()/mult_binwidth);
        if (ipt < mult_npt) spectrum[sp][ipt] += 1./(2.*M_PI*part.pT());
        count[sp] += 1.;
      };
      for (int i = 0; i < ev.size(); ++i)
        if (ev[i].isFinal() && ev[i].idAbs() == 1000010020) add(ev[i], ev[i].id() > 0 ? "d" : "dbar");
      // DEUTERON-STUDY (E. Demiray, 2026-11-04): coalescence deuterons in the multiplicity classes
      for (int anti = 0; anti < 2; ++anti)
        for (const auto & dw : anti ? coal_dbar : coal_d) {
          if (std::abs(dw.first.rap()) >= 0.5 || dw.first.pT() <= 0.) continue;
          std::string sp = anti ? "dbar" : "d";
          int ipt = int(dw.first.pT()/mult_binwidth);
          if (ipt < mult_npt) spectrum[sp][ipt] += dw.second/(2.*M_PI*dw.first.pT());
          count[sp] += dw.second;
        }
      for (int i : nuc)
        if (ev[i].idAbs() == 2212) add(ev[i], ev[i].id() > 0 ? "p" : "pbar");
      int row = std::min(nV0, mult_nmax-1);
      mult_events[row] += 1.;
      mult_nch_sum[row] += nch_mid;
      for (const std::string & sp : mult_species) {
        mult_count_sum[sp][row] += count[sp];
        for (int ipt = 0; ipt < mult_npt; ++ipt) {
          mult_spectrum_sum[sp][row*mult_npt+ipt]        += spectrum[sp][ipt];
          mult_spectrum_square_sum[sp][row*mult_npt+ipt] += spectrum[sp][ipt]*spectrum[sp][ipt];
        }
      }
    }
    for (const auto & x : c) { sum[x.first] += x.second; sum2[x.first] += x.second*x.second; }
  }
  pythia.stat();

  // multiplicity classes, same files as WriteMultiplicityClasses in Deuteron_Formation.C:
  // sums over the events of each row (not divided by the number of events)
  if (beams == "pp") {
    std::ofstream ev_out(outdir + "/multiplicity_events.dat");
    ev_out << "# INEL>0 events binned in the V0M charged multiplicity nV0 (V0-A 2.8<=eta<=5.1, "
           << "V0-C -3.7<=eta<=-1.7; arXiv:1902.09290)\n"
           << "# nV0  events  sum_dNch(|eta|<0.5)  sum_dN/dy(|y|<0.5): d  dbar  p  pbar\n";
    for (int row = 0; row < mult_nmax; ++row) {
      if (mult_events[row] <= 0.) continue;
      ev_out << row << " " << mult_events[row] << " " << mult_nch_sum[row];
      for (const std::string & sp : mult_species) ev_out << " " << mult_count_sum[sp][row];
      ev_out << "\n";
    }
    for (const std::string & sp : mult_species) {
      for (bool squares : {false, true}) {
        const std::vector<double> & vals = squares ? mult_spectrum_square_sum[sp] : mult_spectrum_sum[sp];
        std::ofstream o(outdir + "/multiplicity_" + sp + "_invyield_midy" + (squares ? "_sq" : "") + ".dat");
        o << "# rows: nV0 (as multiplicity_events.dat); columns: " << mult_npt << " pT bins of " << mult_binwidth
          << " GeV from 0; " << (squares ? "sum over events of the squared per-event values"
                                        : "sum over events of 1/(2 pi pT) dN/dpT dy, |y|<0.5")
          << " [GeV^-2" << (squares ? ", squared" : "") << "]\n";
        // per unit pT (the per-event values were filled per bin)
        double norm = squares ? mult_binwidth*mult_binwidth : mult_binwidth;
        for (int row = 0; row < mult_nmax; ++row) {
          if (mult_events[row] <= 0.) continue;
          o << row;
          for (int ipt = 0; ipt < mult_npt; ++ipt) o << " " << vals[row*mult_npt+ipt]/norm;
          o << "\n";
        }
      }
    }
  }

  for (const auto & x : h) x.second.write(outdir + "/" + x.first + ".dat", nacc);
  std::ofstream out(outdir + "/summary.dat");
  out << "# Pythia " << PYTHIA_VERSION << (beams == "pp" ? " non-diffractive pp at " : " e+e- -> gamma*/Z -> q qbar at ") << ecm
      << " GeV, " << (coal ? "coalescence " + coal_model + ", p_coal = " + std::to_string(p_coal) + " GeV, b_coal = " + std::to_string(b_coal) + " fm" : "DeuteronProduction 1/sigma_0 = " + std::to_string(invsig0) + " 1/b") << ", events = " << nacc << "\n"
      << "# key  per-event-mean  error\n";
  for (const auto & x : sum) {
    double m = x.second/nacc, v = std::max(0., sum2[x.first]/nacc - m*m);
    out << x.first << "  " << m << "  " << std::sqrt(v/nacc) << "\n";
  }
  return 0;
}
