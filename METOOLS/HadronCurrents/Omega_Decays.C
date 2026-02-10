#include "METOOLS/HadronCurrents/Omega_Decays.H"
#include "METOOLS/HadronCurrents/Vector_Decays.H"
#include "METOOLS/HadronCurrents/Line_Shapes.H"
#include "ATOOLS/Phys/Flavour.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Omega_782_Lineshape::Omega_782_Lineshape() :
  Total_Width_Base(Flavour(kf_omega_782)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar()), pi_0(kf_pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Relevant containers and other inputs
  ///////////////////////////////////////////////////////////////////////////////////
  vector<Flavour>           outflavs;
  vector<Propagator_Base *> props;
  double m_min;
  axis   Q_axis(0,0,1);
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi^+ pi^- pi^0, BR = 89.2 %
  //          It is modelled by assuming omega -> rho pi, with rho -> pi pi.
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { pi_plus, pi_minus, pi_0 };
  Summed_Propagator * chargedRhos   = new Summed_Propagator();
  chargedRhos->Add(new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)),
				   resonance_type::running),
		   Complex(1.,0.));
  Summed_Propagator * neutralRhos   = new Summed_Propagator();
  neutralRhos->Add(new BreitWigner(LineShapes->Get(Flavour(kf_rho_770)),
				   resonance_type::running),
		   Complex(1.,0.));
  props    = { neutralRhos, chargedRhos, chargedRhos };
  m_min    = pi_plus.Mass(true)+pi_0.Mass(true)+pi_minus.Mass(true);
  Q_axis   = axis(100,m_min,m_min+3.);
  Partial_Width_Base * omega2pipipi = new V_PPP(m_inflav,outflavs,0.899);
  omega2pipipi->Init(props,Q_axis);
  m_channels.insert(omega2pipipi);
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi gamma: BR = 8.3%
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_pi), Flavour(kf_photon) };
  Partial_Width_Base * omega2pigamma = new V_PGamma(m_inflav,outflavs,0.08555);
  m_channels.insert(omega2pigamma);
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi: BR = 1.5%
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_pi_plus), Flavour(kf_pi_plus).Bar() };
  Partial_Width_Base * omega2pipi = new V_PP(m_inflav,outflavs,0.015);
  m_channels.insert(omega2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: eta gamma: BR = 0.045%
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_eta), Flavour(kf_photon) };
  Partial_Width_Base * omega2etagamma = new V_PGamma(m_inflav,outflavs,0.00045);
  m_channels.insert(omega2pigamma);
}

Omega_1420_Lineshape::Omega_1420_Lineshape() :
  Total_Width_Base(Flavour(kf_omega_1420)) {
  /*
  Flavour pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar()), pi_0(kf_pi);
  Flavour rho_0(kf_rho_770);
  Flavour rho_plus(kf_rho_770_plus), rho_minus(rho_plus.Bar());
  vector<Flavour> outflavs, props;
  // pi pi pi: BR = 69.9%  (from rho-pi final states)
  outflavs = { pi_plus, pi_minus, pi_0 };
  props    = { rho_0, rho_plus, rho_minus };
  double m_min = pi_plus.Mass(true)+pi_0.Mass(true)+pi_minus.Mass(true);
  axis   Q_axis(100,m_min,m_min+3.);
  Partial_Width_Base * omega2pipipi = new V_PPP(m_inflav,outflavs,0.699);
  omega2pipipi->Init(props,Q_axis);
  m_channels.insert(omega2pipipi);
  // omega(782) pi pi: 30.1%: todo.
  // Idea: model it through omega(782) + rho(770)
  */
}

Omega_1600_Lineshape::Omega_1600_Lineshape() :
  Total_Width_Base(Flavour(kf_omega_1600)) {
  /*
  Flavour pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar()), pi_0(kf_pi);
  Flavour rho_0(kf_rho_770);
  Flavour rho_plus(kf_rho_770_plus), rho_minus(rho_plus.Bar());
  vector<Flavour> outflavs, props;
  // pi pi pi: BR = 65%  (from rho-pi final states)
  outflavs = { pi_plus, pi_minus, pi_0 };
  props    = { rho_0, rho_plus, rho_minus };
  double m_min = pi_plus.Mass(true)+pi_0.Mass(true)+pi_minus.Mass(true);
  axis   Q_axis(100,m_min,m_min+3.);
  Partial_Width_Base * omega2pipipi = new V_PPP(m_inflav,outflavs,0.65);
  omega2pipipi->Init(props,Q_axis);
  m_channels.insert(omega2pipipi);
  // omega(782) pi pi: 35%: todo.
  // Idea: model it through omega(782) + rho(770)
  */
}
