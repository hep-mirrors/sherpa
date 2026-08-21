#include "METOOLS/HadronCurrents/FormFactors/K1_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Vector_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// K1(1270): on-shell rho(770)K (70%) and K*(892)pi (30%) (Frank,
// private communication). "On-shell" here means we use the simple
// 2-body V_VP partial width (fixed daughter masses, no integration
// over the rho/K* lineshape) rather than a full 3-body V_VoffP/V_PPP
// treatment - i.e. the rho and K* are treated as if stable for the
// purpose of the K1 total width, consistent with how V_VP is already
// used for K*(1410)->K*(892)pi in Kstar_Decays.C.
//
// Isospin split follows the same 50/50 pattern already used for
// K*(892)/K*(1410)'s Kpi and K*pi channels: the neutral K1(1270) (same
// isospin slot as K^0) decays to {K^+ rho^-, K^0 rho^0} and
// {K^{*+} pi^-, K^{*0} pi^0}; the charged K1(1270)^+ (same slot as
// K^+) decays to {K^0 rho^+, K^+ rho^0} and {K^{*0} pi^+, K^{*+} pi^0}.
//
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

K1_1270_0_Lineshape::K1_1270_0_Lineshape() :
  Total_Width_Base(Flavour(kf_K_1_1270)) {
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour K_0(kf_K), K_plus(kf_K_plus);
  Flavour rho_0(kf_rho_770), rho_plus(kf_rho_770_plus), rho_minus(rho_plus.Bar());
  Flavour Kstar_0(kf_K_star_892), Kstar_plus(kf_K_star_892_plus),
          Kstar_minus(Kstar_plus.Bar());
  vector<Flavour> outflavs;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: rho(770) K, BR = 70%, split 50/50 by isospin
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { K_plus, rho_minus };
  Partial_Width_Base * K1_2Krho = new V_VP(m_inflav,outflavs,0.35);
  m_channels.insert(K1_2Krho);
  outflavs = { K_0, rho_0 };
  Partial_Width_Base * K1_2piKrho = new V_VP(m_inflav,outflavs,0.35);
  m_channels.insert(K1_2piKrho);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K^*(892) pi, BR = 30%, split 50/50 by isospin
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Kstar_plus, pi_minus };
  Partial_Width_Base * K1_2Kstarpi = new V_VP(m_inflav,outflavs,0.15);
  m_channels.insert(K1_2Kstarpi);
  outflavs = { Kstar_0, pi_0 };
  Partial_Width_Base * K1_2piKstar = new V_VP(m_inflav,outflavs,0.15);
  m_channels.insert(K1_2piKstar);
}

K1_1270_plus_Lineshape::K1_1270_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_K_1_1270_plus)) {
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour K_0(kf_K), K_plus(kf_K_plus);
  Flavour rho_0(kf_rho_770), rho_plus(kf_rho_770_plus), rho_minus(rho_plus.Bar());
  Flavour Kstar_0(kf_K_star_892), Kstar_plus(kf_K_star_892_plus);
  vector<Flavour> outflavs;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: rho(770) K, BR = 70%, split 50/50 by isospin
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { K_0, rho_plus };
  Partial_Width_Base * K1_2Krho = new V_VP(m_inflav,outflavs,0.35);
  m_channels.insert(K1_2Krho);
  outflavs = { K_plus, rho_0 };
  Partial_Width_Base * K1_2piKrho = new V_VP(m_inflav,outflavs,0.35);
  m_channels.insert(K1_2piKrho);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K^*(892) pi, BR = 30%, split 50/50 by isospin
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Kstar_0, pi_plus };
  Partial_Width_Base * K1_2Kstarpi = new V_VP(m_inflav,outflavs,0.15);
  m_channels.insert(K1_2Kstarpi);
  outflavs = { Kstar_plus, pi_0 };
  Partial_Width_Base * K1_2piKstar = new V_VP(m_inflav,outflavs,0.15);
  m_channels.insert(K1_2piKstar);
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// K1(1400): exclusively K*(892)pi (Frank, private communication),
// same on-shell V_VP treatment and 50/50 isospin split as above.
//
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

K1_1400_0_Lineshape::K1_1400_0_Lineshape() :
  Total_Width_Base(Flavour(kf_K_1_1400)) {
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour Kstar_0(kf_K_star_892), Kstar_plus(kf_K_star_892_plus);
  vector<Flavour> outflavs;
  outflavs = { Kstar_plus, pi_minus };
  Partial_Width_Base * K1_2Kstarpi = new V_VP(m_inflav,outflavs,0.5);
  m_channels.insert(K1_2Kstarpi);
  outflavs = { Kstar_0, pi_0 };
  Partial_Width_Base * K1_2piKstar = new V_VP(m_inflav,outflavs,0.5);
  m_channels.insert(K1_2piKstar);
}

K1_1400_plus_Lineshape::K1_1400_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_K_1_1400_plus)) {
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour Kstar_0(kf_K_star_892), Kstar_plus(kf_K_star_892_plus);
  vector<Flavour> outflavs;
  outflavs = { Kstar_0, pi_plus };
  Partial_Width_Base * K1_2Kstarpi = new V_VP(m_inflav,outflavs,0.5);
  m_channels.insert(K1_2Kstarpi);
  outflavs = { Kstar_plus, pi_0 };
  Partial_Width_Base * K1_2piKstar = new V_VP(m_inflav,outflavs,0.5);
  m_channels.insert(K1_2piKstar);
}
