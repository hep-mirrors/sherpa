#include "PHOTONS++/MEs/W_Decay_EW_One_Loop_Functions.H"
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


// One Loop Renormalization Constants and corresponding helper functions for W-decays

// Note the necessary replacements due to the different metric and
// conventions:
// gamma5 -> -gamma5
// Q^2 -> -s
// A_0 -> -A_0, C_0 -> -C_0
// gamma+/- -> 2 PL/R 

// Taken from Bardin and Passarino (Bardin:1999ak)
// Calculations all in Feynman gauge

W_Decay_EW_One_Loop_Functions::W_Decay_EW_One_Loop_Functions
(const double& s, const double& m2, const double& m2p,
 const double& Qf, const double& Qfp,
 const double& If, const double& Ifp,
 const double& mu2, const int& ew) :
  EW_One_Loop_Functions_Base(s, mu2, ew), m_m2(m2), m_m2p(m2p), m_Qf(Qf), m_Qfp(Qfp), m_If(If), m_Ifp(Ifp), m_ew(ew)
{
}

W_Decay_EW_One_Loop_Functions::~W_Decay_EW_One_Loop_Functions() {
}

// Photon exchange between final state particles
// Necessary if both decay products are charged
// Eqs. (5.602)-(5.603) F_{L/R}^(i)
// ind corresponds to ^(i), hel = 0 is left-handed contribution, hel = 1 right-handed
DivArrC W_Decay_EW_One_Loop_Functions::FQED(const int& ind, const int& hel) {
  double m2plus(m_m2+m_m2p),m2minus(m_m2-m_m2p),Gram(-1./4.*Kallen(-m_s,m_m2,m_m2p));
  if (ind == 1) {
    if (hel == 0) {
      return -2.*(-m_s+m2plus)*C_0(m_m2,m_m2p,m_s,m_m2,0.,m_m2p,m_mu2)
	+ B_0(m_s,m_m2,m_m2p,m_mu2)
	- 2.*One
	- (2. + m_m2p*(-m_s-m2minus)/(2.*Gram))*B_ff(m_s,m_m2p,m_m2)
	- (2. + m_m2*(-m_s+m2plus)/(2.*Gram))*B_ff(m_s,m_m2,m_m2p);
    }
    else if (hel == 1) {
      return sqrt(m_m2)*sqrt(m_m2p)/(2.*Gram)*((-m_s+m2minus)*B_ff(m_s,m_m2p,m_m2)
					       +(-m_s-m2minus)*B_ff(m_s,m_m2,m_m2p));
    }
  }	   
  else if (ind == 2) {
    if (hel == 0) {
      return sqrt(m_m2)/(4.*Gram)*((-m_s+m2plus)*B_ff(m_s,m_m2,m_m2p)-2.*m_m2p*B_ff(m_s,m_m2p,m_m2));
    }
    else if (hel == 1) {
      return sqrt(m_m2p)/(4.*Gram)*((-m_s+m2plus)*B_ff(m_s,m_m2p,m_m2)-2.*m_m2*B_ff(m_s,m_m2,m_m2p));
    }
  }
  else if (ind == 3) {
    if (hel == 0) {
      return -0.5*sqrt(m_m2)/m_s*(4.*One+1./Gram*(m_m2p*(m2minus-2.*m_s)*B_ff(m_s,m_m2p,m_m2)
						  +(-m_m2*m2minus+3./2.*m_s*m_m2p+5./2.*m_s*m_m2
						    -3./2.*pow(m_s,2.))*B_ff(m_s,m_m2,m_m2p)));
    }
    else if (hel == 1) {
      return -0.5*sqrt(m_m2)/m_s*(-4.*One-1./Gram*(m_m2p*(-m2minus-2.*m_s)*B_ff(m_s,m_m2p,m_m2)
						   +(m_m2*m2minus+3./2.*m_s*m_m2p+5./2.*m_s*m_m2
						     -3./2.*pow(m_s,2.))*B_ff(m_s,m_m2,m_m2p)));
    }
  }
  return Zero;
}

// Photon coupling to W and charged final state particle
// Eq. (5.606) - M^2 in prefactors already corrected to read s
DivArrC W_Decay_EW_One_Loop_Functions::FAn() {
  return m_Qf*(m_s*C_0(m_m2,m_m2p,m_s,0.,m_m2,muW2,m_mu2)+B_0(m_m2,m_m2,0.,m_mu2))
    -m_Qfp*(m_s*C_0(m_m2,m_m2p,m_s,muW2,m_m2p,0.,m_mu2)+B_0(m_m2p,m_m2p,0.,m_mu2))
    +(m_Qf-m_Qfp)/2.*(-2.*B_0(m_s,muW2,0.,m_mu2)+3.*B_0(0.,0.,muW2,m_mu2));
}

// Z exchange between final state particles
// Eq. (5.571)
DivArrC W_Decay_EW_One_Loop_Functions::FZa() { 
  Complex z(-m_s/muZ2);
  return 2.*muZ2/z*pow((1.-z),2.)*C_0(0.,0.,m_s,0.,muZ2,0.,m_mu2)
    + B_0(m_s,0.,0.,m_mu2)
    + (2./z-4.)*(B_0(m_s,0.,0.,m_mu2)-B_0(0.,0.,muZ2,m_mu2))
    -2.*DivArrC(0.,0.,0.,1.,0.,0.);
}

// Z exchange between W and final state particle
// Eq. (5.615)
DivArrC W_Decay_EW_One_Loop_Functions::FZn() { 
  Complex w(-m_s/muW2), z(-m_s/muZ2);
  return 0.5*(-((1./w-1.)/m_cw2-1.)*muW2*C_0(0.,0.,m_s,muW2,0.,muZ2,m_mu2)
	      +0.5*(1./z+1./w-1.)*B_0(m_s,muW2,muZ2,m_mu2)
	      -(0.5/z-1.)/muZ2*A_0(muZ2,m_mu2)
	      -(0.5/w-1.)/muW2*A_0(muW2,m_mu2)); /// check signs
}

// Renormalization counterterm for left-handed coupling
// Denner:1991kt Eq. (A.13)
DivArrC W_Decay_EW_One_Loop_Functions::CT_L() { 
  return dZe() + m_cw2/m_sw2*dcw() + 0.5*dZW() + 0.5*(dZfermL(m_m2p,m_m2,m_Qfp,m_Ifp) + dZfermL(m_m2,m_m2p,m_Qf,m_If));
}

// Renormalization counterterm for right-handed coupling
DivArrC W_Decay_EW_One_Loop_Functions::CT_R() { 
  // No right handed coupling at LO, hence no counterterm contribution
  return Zero;
}

