#include "METOOLS/HadronCurrents/A1_Decays.H"
#include "METOOLS/HadronCurrents/Vector_Decays.H"
#include "METOOLS/HadronCurrents/Line_Shapes.H"
#include "ATOOLS/Phys/Flavour.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

///////////////////////////////////////////////////////////////////////////////////
// To do: There are two things worth doing here:
// 1. Making parameters a part of the input structure, and translating the
//    form factors to currents.
// 2. Fitting the f_0 branching ratio (and the others), and of the relative admixtures
//    of higher resonances to the propagators to the tau->3 pion and e^+e^- -> 3
//    pions mass spectra
///////////////////////////////////////////////////////////////////////////////////


A1_1260_0_Lineshape::A1_1260_0_Lineshape() :
  Total_Width_Base(Flavour(kf_a_1_1260)) {
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
  // Channel: {rho(770) + rho(1450)} pi, rho -> pi pi
  //          PDG has 60% for rho(770) and 56% for rho(1450)
  //          Will add them here arrive at 92.5%
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * rho1450_width = new Total_Width_Base(Flavour(kf_rho_1450_plus));
  Summed_Propagator * chargedRhos   = new Summed_Propagator();
  chargedRhos->Add(new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)),
				   resonance_type::running),
		   Complex(1.,0.));
  chargedRhos->Add(new BreitWigner(rho1450_width,
  				   resonance_type::fixed),
  		   Complex(-0.145,0.));
  outflavs  = { pi_plus, pi_0, pi_minus };
  props     = { chargedRhos, chargedRhos };
  m_min     = pi_plus.Mass(true)+pi_0.Mass(true)+pi_minus.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12pipipi = new V_PPP(m_inflav,outflavs,0.925);
  a12pipipi->Init(props,Q_axis);
  m_channels.insert(a12pipipi);
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: {f_0(500) + f_0(1370)} pi, f_0 -> pi pi
  //          PDG has 18% for f_0(500) and 7% for f_0(1370)
  //          However: f_0's decay into 2/3-1/3 mix of charged and neutral pions --
  //                   as they are iso-scalars -- and a_1^0 -> 3 pi_0's hasn't
  //                   been seen (BR<1%).  In addition a recent theory paper argues
  //                   that the overall BR is about 3-4% (see equation 38 in
  //                   https://arxiv.org/pdf/2107.07439).
  //          Will assume an overall BR for the f_0's of 4.5%.
  // Note: due to spin structure of propagators, interference of f_0/rho_0 will
  //       probably be quite small - and ignored here
  //
  // Note(2): the lines below are for a Breit-Wigner smeared f0, mainly to
  //          check factors etc.
  //
  //    Flavour f_0(kf_f_0_600);
  //    Summed_Propagator * f0s = new Summed_Propagator();
  //    f0s->Add(new BreitWigner(LineShapes->Get(f_0),
  //				 resonance_type::running),
  //		 Complex(1.,0.));
  //	outflavs  = { pi_0, pi_plus, pi_minus };
  //	props     = { f0s };
  //	m_min     = pi_0.Mass(true)+2.*pi_minus.Mass(true);
  //	Q_axis    = axis(500,m_min,m_min+3.);
  //	Partial_Width_Base * a12f0_500_pi = new V_PoffP(m_inflav,outflavs,0.03);
  //	a12f0_500_pi->Init(props,Q_axis);
  //	m_channels.insert(a12f0_500_pi);
  //
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour f_0(kf_f_0_600);
  Total_Width_Base  * f0_1370_width = new Total_Width_Base(Flavour(kf_f_0_1370));
  Summed_Propagator * f0s = new Summed_Propagator();
  f0s->Add(new BreitWigner(LineShapes->Get(f_0),
			   resonance_type::running),
	   Complex(1.,0.));
  f0s->Add(new BreitWigner(f0_1370_width,
			   resonance_type::fixed),
	   Complex(-0.3,0.));
  outflavs  = { pi_0, pi_plus, pi_minus };
  props     = { f0s };
  m_min     = pi_plus.Mass(true)+pi_0.Mass(true)+pi_minus.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12f0_500_pi = new V_PPP_scalar(m_inflav,outflavs,0.03);
  a12f0_500_pi->Init(props,Q_axis);
  m_channels.insert(a12f0_500_pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Same again, but for the neutral pions
  ///////////////////////////////////////////////////////////////////////////////////
  Summed_Propagator * f0Ns = new Summed_Propagator();
  f0Ns->Add(new BreitWigner(LineShapes->Get(f_0),
			    resonance_type::running),
	    Complex(1.,0.));
  f0Ns->Add(new BreitWigner(f0_1370_width,
			    resonance_type::fixed),
	   Complex(-0.3,0.));
  outflavs  = { pi_0, pi_0, pi_0 };
  props     = { f0Ns, f0Ns, f0Ns };
  m_min     = 3.*pi_0.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12f0N_500_pi = new V_PPP_scalar(m_inflav,outflavs,0.015);
  a12f0N_500_pi->Init(props,Q_axis);
  m_channels.insert(a12f0N_500_pi);
  // todo: add 3% K*(892)K etc.
}

A1_1260_plus_Lineshape::A1_1260_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_a_1_1260_plus)) {
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
  // Channel: {rho(770) + rho(1450)} pi, rho -> pi pi
  //          PDG has 60% for rho(770) and 56% for rho(1450)
  //          Will add them here arrive at 92.5%
  // This channel comes in two varieties: pi^+ pi^+ pi^- and pi^+ pi^0 pi^0 each
  // with equal probability.  Starting with the latter.
  ///////////////////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * rho1450_width = new Total_Width_Base(Flavour(kf_rho_1450_plus));
  Summed_Propagator * chargedRhos   = new Summed_Propagator();
  chargedRhos->Add(new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)),
				   resonance_type::running),
		   Complex(1.,0.));
  chargedRhos->Add(new BreitWigner(rho1450_width,
  				   resonance_type::fixed),
  		   Complex(-0.145,0.));
  outflavs  = { pi_0, pi_plus, pi_0 };
  props     = { chargedRhos, chargedRhos };
  m_min     = pi_plus.Mass(true)+2.*pi_0.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12pipipi = new V_PPP(m_inflav,outflavs,0.4625);
  a12pipipi->Init(props,Q_axis);
  m_channels.insert(a12pipipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Same again, but for the pi^+ pi^- pi^+ channel
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * rho1450N_width = new Total_Width_Base(Flavour(kf_rho_1450_plus));
  Summed_Propagator * neutralRhos    = new Summed_Propagator();
  neutralRhos->Add(new BreitWigner(LineShapes->Get(Flavour(kf_rho_770_plus)),
				   resonance_type::running),
		   Complex(1.,0.));
  neutralRhos->Add(new BreitWigner(rho1450_width,
  				   resonance_type::fixed),
  		   Complex(-0.145,0.));
  outflavs  = { pi_plus, pi_minus, pi_plus };
  props     = { neutralRhos, neutralRhos };
  m_min     = 2.*pi_plus.Mass(true)+pi_minus.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12pipipiN = new V_PPP(m_inflav,outflavs,0.4625);
  a12pipipiN->Init(props,Q_axis);
  m_channels.insert(a12pipipiN);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: {f_0(500) + f_0(1370)} pi, f_0 -> pi pi
  // See above, neutral a_1 width.
  // Starting with the f_0's decaying into neutral pions
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour f_0(kf_f_0_600);
  Total_Width_Base  * f0_1370_width = new Total_Width_Base(Flavour(kf_f_0_1370));
  Summed_Propagator * f0s = new Summed_Propagator();
  f0s->Add(new BreitWigner(LineShapes->Get(f_0),
			   resonance_type::running),
	   Complex(1.,0.));
  f0s->Add(new BreitWigner(f0_1370_width,
			   resonance_type::fixed),
	   Complex(-0.3,0.));
  outflavs  = { pi_plus, pi_0, pi_0 };
  props     = { f0s };
  m_min     = pi_plus.Mass(true)+2.*pi_0.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12f0_500_pi = new V_PPP_scalar(m_inflav,outflavs,0.03);
  a12f0_500_pi->Init(props,Q_axis);
  m_channels.insert(a12f0_500_pi);
  
  Summed_Propagator * f0Ns = new Summed_Propagator();
  f0Ns->Add(new BreitWigner(LineShapes->Get(f_0),
			    resonance_type::running),
	    Complex(1.,0.));
  f0Ns->Add(new BreitWigner(f0_1370_width,
			    resonance_type::fixed),
	   Complex(-0.3,0.));
  outflavs  = { pi_plus, pi_plus, pi_minus };
  props     = { f0Ns, f0Ns };
  m_min     = 3.*pi_0.Mass(true);
  Q_axis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * a12f0N_500_pi = new V_PPP_scalar(m_inflav,outflavs,0.015);
  a12f0N_500_pi->Init(props,Q_axis);
  m_channels.insert(a12f0N_500_pi);
  // todo: add 2% K*(892)K
}
