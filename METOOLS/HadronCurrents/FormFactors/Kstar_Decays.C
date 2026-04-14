#include "METOOLS/HadronCurrents/FormFactors/Kstar_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Vector_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Kstar_892_0_Lineshape::Kstar_892_0_Lineshape() :
  Total_Width_Base(Flavour(kf_K_star_892)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  vector<Flavour> outflavs;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K pi, BR = 99.754%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_pi_plus).Bar() };
  Partial_Width_Base * K_star2Kpi = new V_PP(m_inflav,outflavs,0.49877);
  m_channels.insert(K_star2Kpi);
  outflavs = { Flavour(kf_K), Flavour(kf_pi) };
  Partial_Width_Base * K_star2piK = new V_PP(m_inflav,outflavs,0.49877);
  m_channels.insert(K_star2piK);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K gamma, BR = 0.246%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K), Flavour(kf_photon) };
  Partial_Width_Base * K_star2pigamma = new V_PGamma(m_inflav,outflavs,0.00246);
  m_channels.insert(K_star2pigamma);
}

Kstar_892_plus_Lineshape::Kstar_892_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_K_star_892_plus)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  vector<Flavour> outflavs;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K pi, BR = 99.902
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_pi) };
  Partial_Width_Base * K_star2Kpi = new V_PP(m_inflav,outflavs,0.49951);
  m_channels.insert(K_star2Kpi);
  outflavs = { Flavour(kf_K), Flavour(kf_pi_plus) };
  Partial_Width_Base * K_star2piK = new V_PP(m_inflav,outflavs,0.49951);
  m_channels.insert(K_star2piK);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi gamma, BR = 0.05%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_photon) };
  Partial_Width_Base * K_star2pigamma = new V_PGamma(m_inflav,outflavs,0.00098);
  m_channels.insert(K_star2pigamma); 
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// Below for the K_star(1410)'s
//
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// Overall approximate BR's (from PDG, all in %):
// -   K^*(892) pi:      > 21
// -   K pi                6.6
// -   K rho             < 7%
// -   K phi               seen
// For the moment I will assume K^* pi is (100-6.6)% and ignore K rho and K phi 
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
Kstar_1410_0_Lineshape::Kstar_1410_0_Lineshape() :
  Total_Width_Base(Flavour(kf_K_star_1410)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour K_0(kf_K), K_plus(kf_K_plus);
  vector<Flavour> outflavs;
  vector<Propagator_Base *> props;
  axis   Qaxis(100,0.,1.);
  double m_min;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K pi, BR = 6.6%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { K_plus, pi_minus };
  Partial_Width_Base * K_star2Kpi = new V_PP(m_inflav,outflavs,0.033);
  m_channels.insert(K_star2Kpi);
  outflavs = { K_0, pi_0 };
  Partial_Width_Base * K_star2piK = new V_PP(m_inflav,outflavs,0.033);
  m_channels.insert(K_star2piK);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K^*(892) pi, BR = 93.4%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K^*(892) pi, BR = 93.4%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour Kstar_0(kf_K_star_892), Kstar_plus(kf_K_star_892_plus);
  outflavs = { Kstar_plus, pi_minus };
  Partial_Width_Base * K_star2Kstarpi = new V_VP(m_inflav,outflavs,0.467);
  m_channels.insert(K_star2Kstarpi);
  outflavs = { Kstar_0, pi_0 };
  Partial_Width_Base * K_star2piKstar = new V_VP(m_inflav,outflavs,0.467);
  m_channels.insert(K_star2piKstar);
}

Kstar_1410_plus_Lineshape::Kstar_1410_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_K_star_1410_plus)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour K_0(kf_K), K_plus(kf_K_plus);
  vector<Flavour> outflavs;
  vector<Propagator_Base *> props;
  axis   Qaxis(100,0.,1.);
  double m_min;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K pi, BR = 6.6%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { K_plus, pi_minus };
  Partial_Width_Base * K_star2Kpi = new V_PP(m_inflav,outflavs,0.033);
  m_channels.insert(K_star2Kpi);
  outflavs = { K_0, pi_0 };
  Partial_Width_Base * K_star2piK = new V_PP(m_inflav,outflavs,0.033);
  m_channels.insert(K_star2piK);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K^*(892) pi, BR = 93.4%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour Kstar_0(kf_K_star_892), Kstar_plus(kf_K_star_892_plus);
  outflavs = { Kstar_plus, pi_0 };
  Partial_Width_Base * K_star2Kstarpi = new V_VP(m_inflav,outflavs,0.467);
  m_channels.insert(K_star2Kstarpi);
  outflavs = { Kstar_0, pi_plus };
  Partial_Width_Base * K_star2piKstar = new V_VP(m_inflav,outflavs,0.467);
  m_channels.insert(K_star2piKstar);
}

