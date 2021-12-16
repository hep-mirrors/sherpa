#include "PHOTONS++/MEs/Z_Decay_EW_One_Loop_Functions.H"
#include "PHOTONS++/MEs/EW_One_Loop_Functions_Base.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "ATOOLS/Math/MyComplex.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "METOOLS/Loops/Divergence_Array.H"

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define B_1(A,B,C,M)        PV_Bubble_1(A,B,C,M)
#define B_0p(A,B,C,M)       Master_Bubble_Prime(A,B,C,M)
#define B_1p(A,B,C,M)       PV_Bubble_1_Prime(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)  Master_Triangle(A,B,C,D,E,F,M)


using namespace ATOOLS;
using namespace METOOLS;
using namespace PHOTONS;

// Vertices and Counterterms for Z->ll, all in Feynman gauge
// Taken (and translated) from Bardin/Passarino

// Note the necessary replacements due to the different metric and
// conventions:
// gamma5 -> -gamma5
// Q^2 -> -s
// A_0 -> -A_0, C_0 -> -C_0
// gamma+/- -> 2 PL/R 

// Validated against OpenLoops 02.06.2017
Z_Decay_EW_One_Loop_Functions::Z_Decay_EW_One_Loop_Functions
(const double& s, const double& m2, const double& m2p,
 const double& Qf, const double& If, const double& mu2,
 const int& ew) :
  EW_One_Loop_Functions_Base(s, mu2, ew), m_m2(m2), m_m2p(m2p),
  m_Qf(Qf), m_If(If), m_ew(ew)
{
  wf = m_m2/muW2;
  zf = m_m2/muZ2;
  wfp = m_m2p/muW2;
  zfp = m_m2p/muZ2;
  beta2 = 1. - wfp;
  kappa = -0.5*beta2*(3.-beta2)*muW2/m_s; 
}

Z_Decay_EW_One_Loop_Functions::~Z_Decay_EW_One_Loop_Functions() {
}

// Photon exchange between final state particles
// Eq. (5.568)
DivArrC Z_Decay_EW_One_Loop_Functions::FAa() {
  return -2.*(m_s-2*m_m2)*C_0(m_m2,m_m2,m_s,m_m2,0.,m_m2,m_mu2) + B_0(m_s,m_m2,m_m2,m_mu2)
    - 4.*B_ff(m_s,m_m2,m_m2) - 2.*One;
}

// Additional mass dependent photon pieces (5.567), (5.588)
DivArrC Z_Decay_EW_One_Loop_Functions::FA1() {
  return -2.*m_s*m_m2/(m_s*(4.*m_m2-m_s)/4.)*B_ff(m_s,m_m2,0.);
}

DivArrC Z_Decay_EW_One_Loop_Functions::FA3() {
  return -sqrt(m_m2)/(m_s*(4.*m_m2-m_s)/4.)*((4.*m_m2-3.*m_s)/2.*B_ff(m_s,m_m2,0.)
					     + (4.*m_m2 - m_s)*One);
}

DivArrC Z_Decay_EW_One_Loop_Functions::FV2() {
  return m_s*sqrt(m_m2)/(m_s*(4.*m_m2-m_s)/2.)*B_ff(m_s,m_m2,m_m2);
}



// Z exchange between final state particles
// Eq.(5.571)
DivArrC Z_Decay_EW_One_Loop_Functions::FZa() {
  Complex z(-m_s/muZ2);
  return 2.*muZ2/z*pow((1.-z),2.)*C_0(0.,0.,m_s,0.,muZ2,0.,m_mu2)
    + B_0(m_s,0.,0.,m_mu2)
    + (2./z-4.)*(B_0(m_s,0.,0.,m_mu2)-B_0(0.,0.,muZ2,m_mu2))
    -2.*One;
}

// W exchange between final state particles
// Eq. (5.577)
DivArrC Z_Decay_EW_One_Loop_Functions::FWa() {
  Complex w(-m_s/muW2);
  return -(-2.*beta2*kappa + 3. + pow(beta2,2.) -2.*w)*muW2*C_0(0.,0.,m_s,m_m2p,muW2,m_m2p,m_mu2)
    +(3.-beta2)/2.*B_0(m_s,m_m2p,m_m2p,m_mu2)
    +2.*(kappa-2.)*(B_0(m_s,m_m2p,m_m2p,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
    -(2.+0.5*wfp)*One;
}

// mass-dependent piece if isopartner is massive
// Eq. (5.593)
DivArrC Z_Decay_EW_One_Loop_Functions::FWabar() {
  Complex w(-m_s/muW2);
  return wfp*(-(pow(beta2,2.)/w-2.)*muW2*C_0(0.,0.,m_s,m_m2p,muW2,m_m2p,m_mu2)
	      - 0.5*B_0(m_s,m_m2p,m_m2p,m_mu2)
	      - beta2/w*(B_0(m_s,m_m2p,m_m2p,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
	      + 0.5*One);
}

// Z coupling to two Ws
// Eq. (5.581)
DivArrC Z_Decay_EW_One_Loop_Functions::FWn() {
  Complex w(-m_s/muW2);
  return -(-2.*beta2*kappa+3.+pow(beta2,2.))*muW2*C_0(0.,0.,m_s,muW2,m_m2p,muW2,m_mu2)
    -2.*(kappa-2.)*(B_0(m_s,muW2,muW2,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
    -(3.+0.5*wfp)*B_0(m_s,muW2,muW2,m_mu2)
    -0.5*wfp*One;
}

// mass-dependent piece if isopartner is massive
// Eq. (5.596)
DivArrC Z_Decay_EW_One_Loop_Functions::FWnbar() {
  Complex w(-m_s/muW2);
  return 0.5/m_cw2*wfp*(-(pow(beta2,2.)/w+4.-wfp)*muW2*C_0(0.,0.,m_s,muW2,m_m2p,muW2,m_mu2)
			+0.5*B_0(m_s,muW2,muW2,m_mu2)
			+beta2/w*(B_0(m_s,muW2,muW2,m_mu2)-B_0(0.,m_m2p,muW2,m_mu2))
			+0.5*One);
}

// Counterterm contribution to right-handed coupling (proportional to PR) 
// Denner:1991kt Eq. (A.13)
DivArrC Z_Decay_EW_One_Loop_Functions::CT_R()
{
  Complex pref = Complex(0.,1.)*m_e/(m_sw*m_cw);
  Complex gR = -m_sw2*m_Qf;

  // msg_Debugging() <<
  //   "Counterterm R UV\nZ Self Energy: " << pref*gR*1./2.*dZZZ().UV() <<
  //   "\nZA transition: " << -1./2.*Complex(0.,1.)*m_e*m_Qf*dZAZ().UV() <<
  //   "\nCharge ren: " << pref*gR*real(dZe()).UV() <<
  //   "\ncW: " << -pref*gR*dcw().UV()/m_sw2 <<
  //   "\nfermion WF: " << pref*gR*dZfermR(m_m2,m_m2p,m_Qf,m_If).UV() <<
  //   "\nsum: " << pref*gR*(1./2.*dZZZ().UV()
  // 			  +real(dZe()).UV()
  // 			  -dcw().UV()/m_sw2
  // 			  +dZfermR(m_m2,m_m2p,m_Qf,m_If).UV())
  //   -1./2.*Complex(0.,1.)*m_e*m_Qf*dZAZ().UV() << "\n\n";

  // Only fermion counterterm relevant if pure QED
  // Else full form necessary
  if (!m_ew) return pref*gR*(dZfermR(m_m2,m_m2p,m_Qf,m_If));
  else return pref*gR*(
		       1./2.*dZZZ()
		       +real(dZe())
		       -dcw()/m_sw2
		       +dZfermR(m_m2,m_m2p,m_Qf,m_If))
	 -1./2.*Complex(0.,1.)*m_e*m_Qf*dZAZ();
}

// Counterterm contribution to left-handed coupling (proportional to PL)
// Denner:1991kt Eq. (A.13)
DivArrC Z_Decay_EW_One_Loop_Functions::CT_L()
{
  Complex pref = Complex(0.,1.)*m_e/(m_sw*m_cw);
  Complex gL = m_If-m_sw2*m_Qf;

  // msg_Debugging() << 
  //   "pref: " << pref << " e " << m_e << "\n" <<
  //   "Counterterm L UV\nZ Self Energy: " << pref*gL*1./2.*dZZZ().UV() <<
  //   "\nZA transition: " << -1./2.*Complex(0.,1.)*m_Qf*dZAZ().UV() <<
  //   "\nCharge ren: " << pref*gL*real(dZe()).UV() <<
  //   "\ncW: " << pref*dcw().UV()*(m_If*pow(m_cw/m_sw,2.)+(m_Qf-m_If)) <<
  //   "\nfermion WF: " << pref*gL*dZfermL(m_m2,m_m2p,m_Qf,m_If).UV() <<
  //   "\nsum: " << (pref*gL*(1./2.*dZZZ()
  // 			   +real(dZe())
  // 			   +dZfermL(m_m2,m_m2p,m_Qf,m_If))
  // 		  +pref*dcw()*(m_If*(m_cw2-m_sw2)/m_sw2+m_Qf)
  // 		  -1./2.*Complex(0.,1.)*m_Qf*dZAZ()).UV() << "\n\n";

  // Only fermion counterterm relevant if pure QED
  // Else full form necessary
  if (!m_ew) return pref*gL*(dZfermL(m_m2,m_m2p,m_Qf,m_If));
  else return pref*gL*(
		       1./2.*dZZZ()
		       +real(dZe())
		       +dZfermL(m_m2,m_m2p,m_Qf,m_If))
	 +pref*dcw()*(m_If*(m_cw2-m_sw2)/m_sw2+m_Qf)
	 -1./2.*Complex(0.,1.)*m_e*m_Qf*dZAZ();
}
