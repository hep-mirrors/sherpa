#include "METOOLS/HadronCurrents/Rho_Decays.H"
#include "METOOLS/HadronCurrents/Vector_Decays.H"
#include "METOOLS/HadronCurrents/Line_Shapes.H"


using namespace METOOLS;
using namespace ATOOLS;
using namespace std;

Rho_770_0_Lineshape::Rho_770_0_Lineshape() :
  Total_Width_Base(Flavour(kf_rho_770)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  vector<Flavour> outflavs;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi, BR = 99.92%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_pi_plus), Flavour(kf_pi_plus).Bar() };
  Partial_Width_Base * rho2pipi = new V_PP(m_inflav,outflavs,0.9992);
  m_channels.insert(rho2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi gamma, BR = 0.05%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_pi), Flavour(kf_photon) };
  Partial_Width_Base * rho2pigamma = new V_PGamma(m_inflav,outflavs,0.0005);
  m_channels.insert(rho2pigamma);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: eta gamma, BR = 0.03%
  // Note: maybe it should use some off-shell matrix element and integrator
  //       here, as the eta is not strictly on-shell at all.
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_eta), Flavour(kf_photon) };
  Partial_Width_Base * rho2etagamma = new V_PGamma(m_inflav,outflavs,0.0003);
  m_channels.insert(rho2pigamma);
}

Rho_770_plus_Lineshape::Rho_770_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_rho_770_plus)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  vector<Flavour> outflavs;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi, BR = 99.95%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_pi_plus), Flavour(kf_pi) };
  Partial_Width_Base * rho2pipi = new V_PP(m_inflav,outflavs,0.9995);
  m_channels.insert(rho2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi gamma, BR = 0.05%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_pi_plus), Flavour(kf_photon) };
  Partial_Width_Base * rho2pigamma = new V_PGamma(m_inflav,outflavs,0.0005);
  m_channels.insert(rho2pigamma); 
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// Below for the rho(1450)'s
//
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
// Overall approximate BR's (from PDG, all in %):
// -   omega pi:      21 
// All other 4 pi final states are normalised to BR(4 pi - omega pi):
// -   a_1 pi:        27 +/-  8   --> 37.5
// -   h_1 pi:         8 +/-  4   --> 11
// -   pi(1300) pi:   37 +/- 13   --> 51.5
// -   rho rho:       11 +/-  5   -->  0
// -   rho (pi pi)_S: 17 +/-  9   -->  0
// At the moment, we won't implement the latter two modes, and therefore
// scale up the a_1, h_1, and pi(1300) modes by about 1.4
// Also:
// - BR(pi pi)/(4 pi - omega pi):  37   +/- 10
// - BR(K K)/BR(pi pi):            30.7 +/- 8.4 +/- 8.2
// Adding those two to the mix of "omega-subtracted" rho decays results in
// about 48% (37%+11%), ignoring all others, we arrive at 148% which need to
// be rescaled to 100%-21%, i.e. all must be multipled by 0.533.
// Therefore [in square brackets ref. values from theory in 2208.10329]:
// -   a_1 pi:        27 +/-  8   --> 37.5    --> 20.0%  [ 3.75%]
// -   h_1 pi:         8 +/-  4   --> 11      -->  5.8%  [ 2.25%]
// -   pi(1300) pi:   37 +/- 13   --> 51.5    --> 27.5%  [ --]
// -   rho rho:       11 +/-  5   -->  0      -->  0.0%  [ --]
// -   rho (pi pi)_S: 17 +/-  9   -->  0      -->  0.0%  [ --]
// -   pi pi:                                 --> 19.7%  [ 7.75%]
// -   K K:                                   -->  5.8%  [ 8.00%]
// -   omega pi                               --> 21.2%  [55.00%]
// -   K* K                                              [ 7.6%]
// -   eta rho                                           [14.50%]
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
Rho_1450_0_Lineshape::Rho_1450_0_Lineshape() :
  Total_Width_Base(Flavour(kf_rho_1450)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  vector<Flavour> outflavs;
  vector<Propagator_Base *> props;
  axis   Qaxis(100,0.,1.);
  double m_min;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi, BR = 32% of omega pi = 7%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { pi_plus, pi_minus };
  Partial_Width_Base * rho2pipi = new V_PP(m_inflav,outflavs,0.197);
  m_channels.insert(rho2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K K, BR = 30% of pi pi = 5.8% for two channels
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_K_plus).Bar() };
  Partial_Width_Base * rho2KKC = new V_PP(m_inflav,outflavs,0.029);
  m_channels.insert(rho2KKC);
  outflavs = { Flavour(kf_K), Flavour(kf_K).Bar() };
  Partial_Width_Base * rho2KKN = new V_PP(m_inflav,outflavs,0.029);
  m_channels.insert(rho2KKN);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: omega(782) pi, BR = 21%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour omega(kf_omega_782);
  Summed_Propagator * omegas = new Summed_Propagator();
  omegas->Add(new BreitWigner(LineShapes->Get(omega),
			      resonance_type::running),
	      Complex(1.,0.));
  props     = { omegas };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_0 };
  m_min     = 2.*(pi_0.HadMass()+pi_plus.HadMass());
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2omegapi = new V_VoffP(m_inflav,outflavs,0.212);
  rho2omegapi->Init(props, Qaxis);
  m_channels.insert(rho2omegapi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: a_1(1260) pi, BR = 27% of 4pi (about 100&) - a_1 pi = 21%,
  //          which comes in two(?) combinations of 11.5% each:
  //          * pi^- a_1^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ a_1^- (-> pi^- pi^0 pi^0, pi^- pi^- pi^+)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all three modes.
  //          (Since there is no evidence for a_1^0 -> pi^0 rho^0 I assume
  //           that there is also no rho(1450)^0 -> a_1^0 pi^0 decay.)
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * a11260_width = new Total_Width_Base(Flavour(kf_a_1_1260_plus));
  Summed_Propagator * a1s = new Summed_Propagator();
  a1s->Add(new BreitWigner(a11260_width,resonance_type::fixed),
	   Complex(1.,0.));
  props     = { a1s };
  outflavs  = { pi_minus, pi_plus, pi_minus, pi_plus };
  m_min     = 4.*pi_plus.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2a1pi = new V_AoffP(m_inflav,outflavs,0.20);
  rho2a1pi->Init(props, Qaxis);
  m_channels.insert(rho2a1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: h_1(1170) pi, BR = 8% 
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * h11170_width = new Total_Width_Base(Flavour(kf_h_1_1170));
  Summed_Propagator * h1s = new Summed_Propagator();
  h1s->Add(new BreitWigner(h11170_width, resonance_type::fixed),
	   Complex(1.,0.));
  props     = { h1s };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_0 };
  m_min     = 2.*(pi_0.HadMass()+pi_plus.HadMass());
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2h1pi = new V_AoffP(m_inflav,outflavs,0.058);
  rho2h1pi->Init(props, Qaxis);
  m_channels.insert(rho2h1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi(1300) pi, BR = 37% 
  //          which comes in two(?) combinations of 18.5% each:
  //          * pi^- pi(1300)^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ pi(1300)^- (-> pi^- pi^0 pi^0, pi^- pi^- pi^+)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all modes.
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * pi1300_width = new Total_Width_Base(Flavour(kf_pi_1300_plus));
  Summed_Propagator * pi1300s = new Summed_Propagator();
  pi1300s->Add(new BreitWigner(pi1300_width, resonance_type::fixed),
	       Complex(1.,0.));
  props    = { pi1300s };
  outflavs = { pi_plus, pi_minus, pi_plus, pi_minus };
  m_min    = 4.*pi_plus.Mass(true);
  Qaxis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * rho2pi1300pi = new V_PoffP(m_inflav,outflavs,0.275);
  rho2pi1300pi->Init(props, Qaxis);
  m_channels.insert(rho2pi1300pi);
}

Rho_1450_plus_Lineshape::Rho_1450_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_rho_1450_plus)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  vector<Flavour> outflavs;
  vector<Propagator_Base *> props;
  axis   Qaxis(100,0.,1.);
  double m_min;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi, BR = 32% of omega pi = 7%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { pi_plus, pi_0 };
  Partial_Width_Base * rho2pipi = new V_PP(m_inflav,outflavs,0.197);
  m_channels.insert(rho2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K K, BR = 30% of pi pi = 2.3% 
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_K).Bar()};
  Partial_Width_Base * rho2KK = new V_PP(m_inflav,outflavs,0.058);
  m_channels.insert(rho2KK);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: omega(782) pi, BR = 21%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour omega(kf_omega_782);
  Summed_Propagator * omegas = new Summed_Propagator();
  omegas->Add(new BreitWigner(LineShapes->Get(omega),
			      resonance_type::running),
	      Complex(1.,0.));
  props     = { omegas };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_plus };
  m_min     = pi_0.HadMass()+3.*pi_plus.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2omegapi = new V_VoffP(m_inflav,outflavs,0.212);
  rho2omegapi->Init(props, Qaxis);
  m_channels.insert(rho2omegapi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: a_1(1260) pi, BR = 27% of 4pi (about 100&) - a_1 pi = 21%,
  //          which comes in two(?) combinations of 11.5% each:
  //          * pi^0 a_1^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ a_1^0 (-> pi^+ pi^- pi^0)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all modes.
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * a11260_width = new Total_Width_Base(Flavour(kf_a_1_1260_plus));
  Summed_Propagator * a1s = new Summed_Propagator();
  a1s->Add(new BreitWigner(a11260_width,resonance_type::fixed),
	   Complex(1.,0.));
  props     = { a1s };
  outflavs  = { pi_plus, pi_minus, pi_plus, pi_0 };
  m_min     = 3.*pi_plus.HadMass()+pi_0.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2a1pi = new V_AoffP(m_inflav,outflavs,0.20);
  rho2a1pi->Init(props, Qaxis);
  m_channels.insert(rho2a1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: h_1(1170) pi, BR = 8% 
  // Note: we ignore interferences of pions from the h_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the h_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * h11170_width = new Total_Width_Base(Flavour(kf_h_1_1170));
  Summed_Propagator * h1s = new Summed_Propagator();
  h1s->Add(new BreitWigner(h11170_width, resonance_type::fixed),
	   Complex(1.,0.));
  props     = { h1s };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_plus };
  m_min     = pi_0.HadMass()+3.*pi_plus.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2h1pi = new V_AoffP(m_inflav,outflavs,0.058);
  rho2h1pi->Init(props, Qaxis);
  m_channels.insert(rho2h1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi(1300) pi, BR = 37% 
  //          which comes in two(?) combinations of 18.5% each:
  //          * pi^0 pi(1300)^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ pi(1300)^0 (-> pi^+ pi^- pi^0)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all modes.
  // Note: we ignore interferences of pions from the pi(1300) with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the pi(1300).
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * pi1300_width = new Total_Width_Base(Flavour(kf_pi_1300));
  Summed_Propagator * pi1300s = new Summed_Propagator();
  pi1300s->Add(new BreitWigner(pi1300_width, resonance_type::fixed),
	       Complex(1.,0.));
  props    = { pi1300s };
  outflavs = { pi_plus, pi_minus, pi_0, pi_0 };
  m_min    = 2.*(pi_0.Mass(true)+pi_plus.Mass(true));
  Qaxis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * rho2pi1300pi = new V_PoffP(m_inflav,outflavs,0.275);
  rho2pi1300pi->Init(props, Qaxis);
  m_channels.insert(rho2pi1300pi);
}


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
//
// Below for the rho(1700)'s
//
// BR(pi^+ pi^-) is given as 10.8% in most recent analysis - alternatively from
// Gamma(pi^+ pi^-) * BR(e^+e^-)=0.13 keV and Gamma(e^+e^-)=7.6 keV we arrive at
// BR(pi^+ pi^-) = 1.7%.  I will use the following methodology, i.e. inferring
// the BR from Gamma * BR(e^+e^-).
// This allows me to arrive at (for the rho^0(1700)):
// -   2(pi^+pi^-)           --> 34.2%
// -   pi^+ pi^- pi^0 pi^0   --> 17.1% (with isospin factor from the all-charged)
// Incidentally, this makes me arrive at the 52% the PDG quotes for the overall
// 4-pion decay channel.  However, the absence of any sizable proportion of
// two-body decays means this is maybe not tenable.  I will therefore rescsale
// the PDG inspired values to arrive at 100%, by a value of 2.21.
//
// Below I add things to the channels where appropriate as a first "PDG-inspired"
// set of branching ratios, rescaling them with the factor 2.21 to the second
// "normalised" BRs and compare with some theory estimates in square brackets.
// -   a_1 pi       -->  8.2%   --> 18.0% [37.6%]
// -   h_1 pi       -->  8.7%   --> 19.1% [34.0%]
// -   pi(1300) pi  --> 15.3%   --> 33.7% [ --  ]
// -   rho rho      -->  4.6%   --> 10.1% [10.6%]
//     (The rho pi pi channel from the PDG adds up to 46% of the 52% 4 pi states;
//     arriving at a total of 37% here, indicating that things are missing
//     or that the central values are to be moved a bit.
//     Also, for the time being I'll ignore the rho rho channel.)
// -   omega pi     -->  1.0%   -->  2.2% [ 5.6%]
// -   K K          -->  0.5%   -->  1.1% [ 3.9%]
// -   K* K         -->  4.0%   -->  8.8% [ 2.7%]
// -   eta rho      -->  1.1%   -->  2.4% [ 3.5%]
//     (for eta rho I use the normalised Gamma of 84 eV) 
// -   pi^+ pi^-    -->  1.7%   -->  3.7% [ --  ]
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

Rho_1700_0_Lineshape::Rho_1700_0_Lineshape() :
  Total_Width_Base(Flavour(kf_rho_1700)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  vector<Flavour> outflavs;
  vector<Propagator_Base *> props;
  axis   Qaxis(100,0.,1.);
  double m_min;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi, BR = 1.7%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { pi_plus, pi_minus };
  Partial_Width_Base * rho2pipi = new V_PP(m_inflav,outflavs,0.017);
  m_channels.insert(rho2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K K, BR = 1.1% for two channels
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_K_plus).Bar() };
  Partial_Width_Base * rho2KKC = new V_PP(m_inflav,outflavs,0.0055);
  m_channels.insert(rho2KKC);
  outflavs = { Flavour(kf_K), Flavour(kf_K).Bar() };
  Partial_Width_Base * rho2KKN = new V_PP(m_inflav,outflavs,0.0055);
  m_channels.insert(rho2KKN);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: omega(782) pi, BR = 2.2%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour omega(kf_omega_782);
  Summed_Propagator * omegas = new Summed_Propagator();
  omegas->Add(new BreitWigner(LineShapes->Get(omega),
			      resonance_type::running),
	      Complex(1.,0.));
  props     = { omegas };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_0 };
  m_min     = 2.*(pi_0.HadMass()+pi_plus.HadMass());
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2omegapi = new V_VoffP(m_inflav,outflavs,0.022);
  rho2omegapi->Init(props, Qaxis);
  m_channels.insert(rho2omegapi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: a_1(1260) pi, BR = 18% 
  //          which comes in two(?) combinations of 9% each:
  //          * pi^- a_1^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ a_1^- (-> pi^- pi^0 pi^0, pi^- pi^- pi^+)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all three modes.
  //          (Since there is no evidence for a_1^0 -> pi^0 rho^0 I assume
  //           that there is also no rho(1700)^0 -> a_1^0 pi^0 decay.)
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * a11260_width = new Total_Width_Base(Flavour(kf_a_1_1260_plus));
  Summed_Propagator * a1s = new Summed_Propagator();
  a1s->Add(new BreitWigner(a11260_width,resonance_type::fixed),
	   Complex(1.,0.));
  props     = { a1s };
  outflavs  = { pi_minus, pi_plus, pi_minus, pi_plus };
  m_min     = 4.*pi_plus.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2a1pi = new V_AoffP(m_inflav,outflavs,0.18);
  rho2a1pi->Init(props, Qaxis);
  m_channels.insert(rho2a1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: h_1(1170) pi, BR = 19.1% 
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * h11170_width = new Total_Width_Base(Flavour(kf_h_1_1170));
  Summed_Propagator * h1s = new Summed_Propagator();
  h1s->Add(new BreitWigner(h11170_width, resonance_type::fixed),
	   Complex(1.,0.));
  props     = { h1s };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_0 };
  m_min     = 2.*(pi_0.HadMass()+pi_plus.HadMass());
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2h1pi = new V_AoffP(m_inflav,outflavs,0.191);
  rho2h1pi->Init(props, Qaxis);
  m_channels.insert(rho2h1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi(1300) pi, BR = 37% 
  //          which comes in two(?) combinations of 18.5% each:
  //          * pi^- pi(1300)^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ pi(1300)^- (-> pi^- pi^0 pi^0, pi^- pi^- pi^+)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all modes.
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * pi1300_width = new Total_Width_Base(Flavour(kf_pi_1300_plus));
  Summed_Propagator * pi1300s = new Summed_Propagator();
  pi1300s->Add(new BreitWigner(pi1300_width, resonance_type::fixed),
	       Complex(1.,0.));
  props    = { pi1300s };
  outflavs = { pi_plus, pi_minus, pi_plus, pi_minus };
  m_min    = 4.*pi_plus.Mass(true);
  Qaxis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * rho2pi1300pi = new V_PoffP(m_inflav,outflavs,0.275);
  rho2pi1300pi->Init(props, Qaxis);
  m_channels.insert(rho2pi1300pi);
}

Rho_1700_plus_Lineshape::Rho_1700_plus_Lineshape() :
  Total_Width_Base(Flavour(kf_rho_1700_plus)) {
  ///////////////////////////////////////////////////////////////////////////////////
  // Declaration of basic flavours and propagators to increase
  // readability
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour pi_0(kf_pi), pi_plus(kf_pi_plus), pi_minus(pi_plus.Bar());
  vector<Flavour> outflavs;
  vector<Propagator_Base *> props;
  axis   Qaxis(100,0.,1.);
  double m_min;
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi pi, BR = 3.7%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { pi_plus, pi_0 };
  Partial_Width_Base * rho2pipi = new V_PP(m_inflav,outflavs,0.037);
  m_channels.insert(rho2pipi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: K K, BR = 1.1%
  ///////////////////////////////////////////////////////////////////////////////////
  outflavs = { Flavour(kf_K_plus), Flavour(kf_K).Bar()};
  Partial_Width_Base * rho2KK = new V_PP(m_inflav,outflavs,0.011);
  m_channels.insert(rho2KK);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: omega(782) pi, BR = 2.2%
  // Note: we ignore interferences of pions from the omega with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width line shape of the omega.
  ///////////////////////////////////////////////////////////////////////////////////
  Flavour omega(kf_omega_782);
  Summed_Propagator * omegas = new Summed_Propagator();
  omegas->Add(new BreitWigner(LineShapes->Get(omega),
			      resonance_type::running),
	      Complex(1.,0.));
  props     = { omegas };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_plus };
  m_min     = pi_0.HadMass()+3.*pi_plus.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2omegapi = new V_VoffP(m_inflav,outflavs,0.022);
  rho2omegapi->Init(props, Qaxis);
  m_channels.insert(rho2omegapi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: a_1(1260) pi, BR = 18%,
  //          which comes in two(?) combinations of 9% each:
  //          * pi^0 a_1^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ a_1^0 (-> pi^+ pi^- pi^0)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all modes.
  // Note: we ignore interferences of pions from the a_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the a_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * a11260_width = new Total_Width_Base(Flavour(kf_a_1_1260_plus));
  Summed_Propagator * a1s = new Summed_Propagator();
  a1s->Add(new BreitWigner(a11260_width,resonance_type::fixed),
	   Complex(1.,0.));
  props     = { a1s };
  outflavs  = { pi_plus, pi_minus, pi_plus, pi_0 };
  m_min     = 3.*pi_plus.HadMass()+pi_0.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2a1pi = new V_AoffP(m_inflav,outflavs,0.18);
  rho2a1pi->Init(props, Qaxis);
  m_channels.insert(rho2a1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: h_1(1170) pi, BR = 19.1% 
  // Note: we ignore interferences of pions from the h_1 with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the h_1.
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * h11170_width = new Total_Width_Base(Flavour(kf_h_1_1170));
  Summed_Propagator * h1s = new Summed_Propagator();
  h1s->Add(new BreitWigner(h11170_width, resonance_type::fixed),
	   Complex(1.,0.));
  props     = { h1s };
  outflavs  = { pi_0, pi_plus, pi_minus, pi_plus };
  m_min     = pi_0.HadMass()+3.*pi_plus.HadMass();
  Qaxis     = axis(100,m_min,m_min+3.);
  Partial_Width_Base * rho2h1pi = new V_AoffP(m_inflav,outflavs,0.191);
  rho2h1pi->Init(props, Qaxis);
  m_channels.insert(rho2h1pi);
  ///////////////////////////////////////////////////////////////////////////////////
  // Channel: pi(1300) pi, BR = 33.7% 
  //          which comes in two(?) combinations of 16.85 each:
  //          * pi^0 pi(1300)^+ (-> pi^+ pi^0 pi^0, pi^+ pi^+ pi^-)
  //          * pi^+ pi(1300)^0 (-> pi^+ pi^- pi^0)
  //          Because they are nearly identical kinematically, we will not
  //          distinguish them in any form and only use the first channel,
  //          but giving it the full BR of all modes.
  // Note: we ignore interferences of pions from the pi(1300) with the "prompt" ones
  //       and use the simple Breit-Wigner fixed-width lineshape of the pi(1300).
  ///////////////////////////////////////////////////////////////////////////////////
  Total_Width_Base  * pi1300_width = new Total_Width_Base(Flavour(kf_pi_1300));
  Summed_Propagator * pi1300s = new Summed_Propagator();
  pi1300s->Add(new BreitWigner(pi1300_width, resonance_type::fixed),
	       Complex(1.,0.));
  props    = { pi1300s };
  outflavs = { pi_plus, pi_minus, pi_0, pi_0 };
  m_min    = 2.*(pi_0.Mass(true)+pi_plus.Mass(true));
  Qaxis    = axis(500,m_min,m_min+3.);
  Partial_Width_Base * rho2pi1300pi = new V_PoffP(m_inflav,outflavs,0.337);
  rho2pi1300pi->Init(props, Qaxis);
  m_channels.insert(rho2pi1300pi);
}
