#include "METOOLS/HadronCurrents/FormFactors/A0_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Scalar_Decays.H"
#include "METOOLS/HadronCurrents/FormFactors/Line_Shapes.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// a0(980): decays to eta pi (100% below the eta' pi threshold) plus
// eta' pi with the SAME coupling once that channel opens (Frank,
// private communication - "I think"). Since eta' pi is kinematically
// CLOSED at the a0(980) pole mass (m_a0=0.980 GeV < m_eta'+m_pi~1.10
// GeV), a normal BR-based construction is ill-defined there: S_PP's
// FixPrefactor() computes m_prefactor = m_partialwidth/(*this)(pole),
// and (*this)(pole)=0 exactly for a channel closed at the pole
// (Partial_Width_Base::operator()'s threshold guard), giving 0/0.
// Instead: build the eta pi channel normally (BR=1, well-defined),
// then build the eta' pi channel and OVERWRITE its (otherwise-NaN)
// prefactor with the eta pi channel's fitted value via SetPrefactor()/
// Prefactor() - this is exactly "same coupling", independent of
// whichever placeholder BR the constructor was given.
//
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

A0_980_0_Lineshape::A0_980_0_Lineshape() :
  Total_Width_Base(Flavour(kf_a_0_980)) {
  Flavour pi_0(kf_pi), eta(kf_eta), eta_prime(kf_eta_prime_958);
  vector<Flavour> outflavs;
  outflavs = { eta, pi_0 };
  S_PP * a0_2etapi = new S_PP(m_inflav,outflavs,1.0);
  m_channels.insert(a0_2etapi);
  outflavs = { eta_prime, pi_0 };
  S_PP * a0_2etaprimepi = new S_PP(m_inflav,outflavs,0.0); // placeholder BR
  a0_2etaprimepi->SetPrefactor(a0_2etapi->Prefactor()); // "same coupling"
  m_channels.insert(a0_2etaprimepi);
}

A0_980_plus_Lineshape::A0_980_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_a_0_980_plus)) {
  Flavour pi_plus(kf_pi_plus), eta(kf_eta), eta_prime(kf_eta_prime_958);
  vector<Flavour> outflavs;
  outflavs = { eta, pi_plus };
  S_PP * a0_2etapi = new S_PP(m_inflav,outflavs,1.0);
  m_channels.insert(a0_2etapi);
  outflavs = { eta_prime, pi_plus };
  S_PP * a0_2etaprimepi = new S_PP(m_inflav,outflavs,0.0);
  a0_2etaprimepi->SetPrefactor(a0_2etapi->Prefactor());
  m_channels.insert(a0_2etaprimepi);
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// a0(1450): {pi eta = 10%, pi eta' = 3%, K Kbar = 10%,
// omega pi pi = 67%} (Frank, private communication). The K Kbar and
// omega-pi-pi channels need an assumed CHARGE COMBINATION that wasn't
// specified explicitly - flagged below; please correct if wrong.
//
// The omega pi pi channel (67% - the dominant one!) has no dedicated
// 3-body S->VPP phase-space/Dalitz treatment in this code; it is
// implemented via Constant_PW (Scalar_Decays.H), which contributes a
// CONSTANT partial width above the omega+2pi threshold and exactly
// zero below it - i.e. it keeps the total-width budget honest but has
// NO real dynamical s-dependence. A proper treatment would need a
// dedicated 3-body phase-space integration analogous to the a1(1260)/
// K1(1270,1400) machinery elsewhere in this package - flag if wanted.
//
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

A0_1450_0_Lineshape::A0_1450_0_Lineshape() :
  Total_Width_Base(Flavour(kf_a_0_1450)) {
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  Flavour eta(kf_eta), eta_prime(kf_eta_prime_958);
  Flavour K_0(kf_K), K_plus(kf_K_plus);
  Flavour omega(kf_omega_782);
  vector<Flavour> outflavs;
  outflavs = { eta, pi_0 };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.10));
  outflavs = { eta_prime, pi_0 };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.03));
  // K Kbar, 10% total - ASSUMED split 50/50 between K+K- and K0K0bar
  // (same isospin-split convention used for phi(1020) elsewhere in
  // this package); please correct if the physical split differs.
  outflavs = { K_plus, K_plus.Bar() };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.05));
  outflavs = { K_0, K_0.Bar() };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.05));
  // omega pi pi, 67% - ASSUMED to be omega pi+ pi- (the charge
  // combination allowed for a neutral parent); see the Constant_PW
  // flag above for why this is only a placeholder for the total width,
  // not a real 3-body treatment.
  outflavs = { omega, pi_plus, pi_minus };
  m_channels.insert(new Constant_PW(m_inflav,outflavs,0.67));
}

A0_1450_plus_Lineshape::A0_1450_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_a_0_1450_plus)) {
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus);
  Flavour eta(kf_eta), eta_prime(kf_eta_prime_958);
  Flavour K_0(kf_K), K_plus(kf_K_plus);
  Flavour omega(kf_omega_782);
  vector<Flavour> outflavs;
  outflavs = { eta, pi_plus };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.10));
  outflavs = { eta_prime, pi_plus };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.03));
  // K Kbar, 10% - ASSUMED entirely K+ K0bar (the only charge-
  // conserving combination for the charged parent, so no 50/50 split
  // here unlike the neutral state above).
  outflavs = { K_plus, K_0.Bar() };
  m_channels.insert(new S_PP(m_inflav,outflavs,0.10));
  // omega pi pi, 67% - ASSUMED omega pi+ pi0 (charge +1 combination).
  outflavs = { omega, pi_plus, pi_0 };
  m_channels.insert(new Constant_PW(m_inflav,outflavs,0.67));
}
