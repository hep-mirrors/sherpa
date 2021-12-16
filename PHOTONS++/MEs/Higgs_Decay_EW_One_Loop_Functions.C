#include "PHOTONS++/MEs/Higgs_Decay_EW_One_Loop_Functions.H"
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

// EW form factors and CT for H -> ll, all in Feynman gauge
// Taken (and translated) from Bardin/Passarino

// Note the necessary replacements due to the different metric and
// conventions:
// gamma5 -> -gamma5
// Q^2 -> -s
// A_0 -> -A_0, C_0 -> -C_0
// gamma+/- -> 2 PL/R 

Higgs_Decay_EW_One_Loop_Functions::Higgs_Decay_EW_One_Loop_Functions
(const double& s, const double& m2, const double& m2p,
 const double& Qf, const double& If, const double& mu2, 
 const int& ew) :
  EW_One_Loop_Functions_Base(s, mu2, ew), m_m2(m2), m_m2p(m2p), m_Qf(Qf), m_If(If)
{
  v = (If - 2.*m_sw2*Qf);
  double a = If;
  // This is defined in Eq. (5.450)
  sigma2 = pow(v,2.) + pow(a,2.);
  // Helper variables, defined in Eq. (5.548)
  w = -m_s/MW2;
  mu2w = MW2/(4.*m_m2-m_s);
  z = -m_s/MZ2;
  wfp = m_m2p/MW2;
  wf = m_m2/MW2;
  wh = MH2/MW2;
}

Higgs_Decay_EW_One_Loop_Functions::~Higgs_Decay_EW_One_Loop_Functions() {
}

// Coefficients needed in form factor calculations
// index i refers to respective inded in f_i
// Eq. (5.547)
Complex Higgs_Decay_EW_One_Loop_Functions::f(const int& i) {
  if (i == 1) {
    return 0.25*(((4.+wfp*(2.+wh))*(1.-wfp)-wf*(10.-4.*wfp-(1.-2.*wfp)*wh))*mu2w
		 +2.+wh*wfp-2.*wf);
  }
  else if (i == 2) {
    // changed pieces to match 5.547 with 5.623
    return 0.25*((4./m_cw4*sigma2-wf*(2./m_cw2-wh))*(1./m_cw2-2.*wf)*mu2w
		 +4.*pow(v,2.)/m_cw4-0.5*wf*(2./m_cw2-wh)); 
  }
  else if (i == 3) {
    return 0.25*wfp*((2.*(2.+wfp)*(1.-wfp)+2.*wf*(1.+2.*wfp-wf))*mu2w-1.);
  }
  else if (i == 4) {
    return 0.25/m_cw2*((sigma2-0.5)*w+wf*(4.*sigma2*(mu2w/m_cw2+0.5)-3./2.));
  }
  else if (i == 5) {
    return 0.25*((4. + wfp*(2. + wh) - wf*(6. - wh))*mu2w + 1.);
  }
  else if (i == 6) {
    // changed last term compared to 5.547 to correspond to 5.623
    return 0.25*((4./m_cw4*sigma2 - wf*(2./m_cw2 - wh))*mu2w + 0.5/m_cw2); 
  }
  else if (i == 7) {
    return -0.25*wfp*(2.*(2. + wfp - wf)*mu2w + 1.);
  }
  else if (i == 8) {
    return -0.25*((2.*(2. + wfp)*(1. - wfp) + wfp*wh - wf*(6. - 2.*wfp - wh))*mu2w + 2.);
  }
  else if (i == 9) {
    return -0.25*((4./m_cw2*sigma2*(1./m_cw2 - wf) - wf*(2./m_cw2 - wh))*mu2w + 2./m_cw2*sigma2);
  }
  else if (i == 10) {
    return 0.25/m_cw2*(sigma2 - 0.5);
  }
  else {
    msg_Out() << METHOD << "Wrong call of f coefficient. i " << i << "\n";
    return 0.;
  }
}

// Coefficients needed in form factor calculations
// index i refers to respective inded in h_i
// Eq. (5.547)
Complex Higgs_Decay_EW_One_Loop_Functions::h(const int& i) {
  if (i == 1) {
    return 3./2.*wf*wh*((0.5*wh-wf)*mu2w-0.25);
  }
  else if (i == 2) {
    return -wf*(1./8.*wh-wf*(wh*mu2w-0.5));
  }
  else if (i == 3) {
    return 3./4.*wf*wh*mu2w;
  }
  else if (i == 4) {
    return -wf*(sigma2/m_cw2 + wf)*mu2w;
  }
  else if (i == 5) {
    return -wf*(3./4.*wh - wf)*mu2w;
  }
  else {
    msg_Out() << METHOD << "Wrong call of h coefficient. i " << i << "\n";
    return 0.;
  }
}

// Photon exchange between final state particles
// Eq. (5.545)
DivArrC Higgs_Decay_EW_One_Loop_Functions::FSQED() {
  return pow(m_Qf,2.)*m_sw2*((m_s-2.*m_m2)*C_0(m_m2,m_m2,m_s,m_m2,0.,m_m2,m_mu2) 
			     -2.*B_0(m_m2,m_m2,0.,m_mu2)
			     +One
			     -4.*m_m2/(4.*m_m2-m_s)*(B_0(m_s,m_m2,m_m2,m_mu2)-B_0(m_m2,m_m2,0.,m_mu2)));
}

// Weak corrections to scalar vertex - vanish if QED only
// Eq. (5.546)
DivArrC Higgs_Decay_EW_One_Loop_Functions::FS() {

  if (m_ew == 0) return Zero;
  else
    return (-MW2*(f(1)*C_0(m_m2,m_m2,m_s,muW2,m_m2p,muW2,m_mu2) 
		  + f(2)*C_0(m_m2,m_m2,m_s,muZ2,m_m2,muZ2,m_mu2)
		  + f(3)*C_0(m_m2,m_m2,m_s,m_m2p,muW2,m_m2p,m_mu2)
		  + f(4)*C_0(m_m2,m_m2,m_s,m_m2,muZ2,m_m2,m_mu2)
		  + h(1)*C_0(m_m2,m_m2,m_s,muH2,m_m2,muH2,m_mu2)
		  + h(2)*C_0(m_m2,m_m2,m_s,m_m2,muH2,m_m2,m_mu2))
	    +f(5)*B_0(m_s,muW2,muW2,m_mu2) + f(6)*B_0(m_s,muZ2,muZ2,m_mu2) + f(7)*B_0(m_s,mut2,mut2,m_mu2)
	    +f(8)*B_0(m_m2,muW2,mut2,m_mu2) + f(9)*B_0(m_m2,muZ2,m_m2,m_mu2) + f(10)*One
	    +h(3)*B_0(m_s,muH2,muH2,m_mu2) + h(4)*B_0(m_s,m_m2,m_m2,m_mu2) + h(5)*B_0(m_m2,muH2,m_m2,m_mu2));
}

// Counterterm contribution prop. to PL
// Denner:1991kt Eq. (A.16)
DivArrC Higgs_Decay_EW_One_Loop_Functions::CT_L() {
  // Only fermionic corrections relevant if QED only
  // Else full expression needed
  if (!m_ew) return -0.5*Complex(0.,m_e)/m_sw*1./sqrt(MW2)*
	       (sqrt(m_m2)*dm(m_m2,m_m2p,m_Qf,m_If)/sqrt(m_m2)
		+0.5*sqrt(m_m2)*(dZfermR(m_m2,m_m2p,m_Qf,m_If)
				 +dZfermL(m_m2,m_m2p,m_Qf,m_If)));
  else return -0.5*Complex(0.,m_e)/m_sw*1./sqrt(MW2)*(sqrt(m_m2)*(dZe()+m_cw2/m_sw2*dcw()
							     +dm(m_m2,m_m2p,m_Qf,m_If)/sqrt(m_m2)
							     -0.5*dMW2()/MW2+0.5*dZH())
						 +0.5*sqrt(m_m2)*(dZfermR(m_m2,m_m2p,m_Qf,m_If)
								  +dZfermL(m_m2,m_m2p,m_Qf,m_If)));
}

// Counterterm contribution prop. to PR
// Denner:1991kt Eq. (A.16)
DivArrC Higgs_Decay_EW_One_Loop_Functions::CT_R() {
  // Only fermionic corrections relevant if QED only
  // Else full expression needed
  if (!m_ew) return -0.5*Complex(0.,m_e)/m_sw*1./sqrt(MW2)*
	       (sqrt(m_m2)*dm(m_m2,m_m2p,m_Qf,m_If)/sqrt(m_m2)
		+0.5*sqrt(m_m2)*(dZfermR(m_m2,m_m2p,m_Qf,m_If)
				 +dZfermL(m_m2,m_m2p,m_Qf,m_If)));
  else return -0.5*Complex(0.,m_e)/m_sw*1./sqrt(MW2)*(sqrt(m_m2)*(dZe()+m_cw2/m_sw2*dcw()
						   +dm(m_m2,m_m2p,m_Qf,m_If)/sqrt(m_m2)
							     -0.5*dMW2()/MW2+0.5*dZH())
						 +0.5*sqrt(m_m2)*(dZfermR(m_m2,m_m2p,m_Qf,m_If)
								  +dZfermL(m_m2,m_m2p,m_Qf,m_If)));
}
