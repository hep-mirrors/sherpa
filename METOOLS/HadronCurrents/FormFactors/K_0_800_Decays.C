#include "METOOLS/HadronCurrents/FormFactors/K_0_800_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Scalar_Decays.H"
#include "ATOOLS/Phys/Flavour.H"

using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////
//
// K_0*(700)/kappa, kf 9000311 (neutral)/9000321 (charged). Decays
// essentially 100% to K pi (kappa sits well below any other open
// channel), split 50/50 by isospin - identical structure to
// Kstar0_1430_0/plus_Lineshape (Kstar_Decays.C), just K0*(1430)'s
// light-nonet partner instead. Uses the same S_PP scalar partial
// width, already exercised for f0(600)->pipi (F0_Decays.C) and
// K0*(1430)->Kpi (Kstar_Decays.C).
//
// IMPORTANT: this class only supplies the DECAY CHANNELS (i.e. how
// the total width is distributed, via FixPrefactor's BR-based
// normalisation). The pole MASS and TOTAL WIDTH themselves come from
// Sherpa's own core particle database (Flavour(kf_K_0_800).HadMass()/
// .Width()), same as every other resonance in this codebase - NOT
// hardcoded here. Those values are NOT part of this HadronCurrents
// tarball and need to be set separately (same situation as the
// a1(1260) mass discussed earlier this session).
//
// Suggested starting values, given kappa's genuinely unsettled
// experimental status (see the header's caveat): M~0.845 GeV,
// Gamma~0.468 GeV are the current PDG "K0*(700)" summary-table
// values (2022 review), quoted there with large, asymmetric
// uncertainties - treat as a starting point to retune, not an
// authoritative number.
//
///////////////////////////////////////////////////////////////////////////////////

K_0_800_0_Lineshape::K_0_800_0_Lineshape() :
  Total_Width_Base(Flavour(kf_K_0_800)) {
  vector<Flavour> outflavs;
  outflavs = { Flavour(kf_K_plus), Flavour(kf_pi_plus).Bar() };
  Partial_Width_Base * K02Kpi = new S_PP(m_inflav,outflavs,0.5);
  m_channels.insert(K02Kpi);
  outflavs = { Flavour(kf_K), Flavour(kf_pi) };
  Partial_Width_Base * K02piK = new S_PP(m_inflav,outflavs,0.5);
  m_channels.insert(K02piK);
}

K_0_800_plus_Lineshape::K_0_800_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_K_0_800_plus)) {
  vector<Flavour> outflavs;
  outflavs = { Flavour(kf_K), Flavour(kf_pi_plus) };
  Partial_Width_Base * K02Kpi = new S_PP(m_inflav,outflavs,0.5);
  m_channels.insert(K02Kpi);
  outflavs = { Flavour(kf_K_plus), Flavour(kf_pi) };
  Partial_Width_Base * K02piK = new S_PP(m_inflav,outflavs,0.5);
  m_channels.insert(K02piK);
}
