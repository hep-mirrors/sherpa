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
#include "PHOTONS++/Main/Photons.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Org/Data_Reader.H"

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define B_1(A,B,C,M)        PV_Bubble_1(A,B,C,M)
#define B_0p(A,B,C,M)       Master_Bubble_Prime(A,B,C,M)
#define B_1p(A,B,C,M)       PV_Bubble_1_Prime(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)  Master_Triangle(A,B,C,D,E,F,M)


using namespace ATOOLS;
using namespace METOOLS;
using namespace PHOTONS;

// EW One Loop Renormalization Constants and corresponding helper functions
// Following Denner's approach (Denner:1991kt) - all equation #'s refer to this unless otherwise described
// Calculations all in Feynman gauge
// All counterterms (except massive fermion WF) validated against OpenLoops 02.06.2017


EW_One_Loop_Functions_Base::EW_One_Loop_Functions_Base(const double& s, const double& mu2,
						       const int& ew) :
  m_alpha(Photons::s_alpha), m_s(s), m_mu2(mu2), m_ew(ew)
{
  // Switchers for debugging of counterterms
  Data_Reader reader(" ",";","#","=");
  // Which parts of the self-energies to calculate. Only use in debugging.
  m_ferm = 1; 
  m_nonferm = 1;
  // Set EW parameters
  double  MW  = Flavour(kf_Wplus).Mass();
  double  MZ  = Flavour(kf_Z).Mass();
  double  MH  = Flavour(kf_h0).Mass();
  double  mt  = Flavour(kf_t).Mass();
  MW2 = pow(MW,2.);
  MZ2 = pow(MZ,2.);
  MH2 = pow(MH,2.);
  mt2 = pow(mt,2.);
  RZ = MZ2/m_s;
  RW = MW2/m_s;
  wh = MH2/MW2;
  zh = MZ2/MW2;
  double  GH  = 0.;///Flavour(kf_h0).Width();
  double  GW  = 0.;///Flavour(kf_Wplus).Width();
  double  GZ  = 0.;///Flavour(kf_Z).Width();
  double  Gt  = 0.;///Flavour(kf_t).Width();
  muW2 = Complex(MW*MW,-GW*MW);
  muZ2 = Complex(MZ*MZ,-GZ*MZ);
  muH2 = Complex(MH*MH,-GH*MH);
  mut2 = Complex(mt*mt,-Gt*mt);
  m_cw2=muW2/muZ2;
  m_sw2=1.-m_cw2;
  m_sw = sqrt((m_sw2));
  m_cw = sqrt((m_cw2));
  m_cw4 = pow(m_cw2,2.);
  m_sw4 = pow(m_sw2,2.);
    
  double GF = 1.16639e-5;
  m_e = sqrt(4.*M_PI*m_alpha);

  // Convenient helper arrays
  One = DivArrC(0.,0.,0.,1.,0.,0.);
  Zero = DivArrC(0.,0.,0.,0.,0.,0.);
}

EW_One_Loop_Functions_Base::~EW_One_Loop_Functions_Base() {
}

// Convenient helper function
// Bardin:1999ak Eq. (5.561)
// Used in derived classes in QED corrections. Note that there are other
// functions called B_ff also in use in that book!
DivArrC EW_One_Loop_Functions_Base::B_ff(const double& p2, const Complex& m12,
					 const Complex& m22) {
  return B_0(p2,m12,m22,m_mu2) - B_0(real(m12),m12,0.,m_mu2);
}

// Kallen lambda function
double EW_One_Loop_Functions_Base::Kallen(const double& x, const double& y, const double& z) {
  return (x*x + y*y + z*z - 2.*x*y - 2.*x*z - 2.*y*z);
}

// Fermion wavefunction counterterm, contribution proportional to PL = (1 - gamma5)/2
// Eq. (3.20)
DivArrC EW_One_Loop_Functions_Base::dZfermL(const double& m2, const double& m2p,
					    const double& Qf, const double& If) {
  // Do not evaluate derivatives of self energies if m2 = 0. (tend to be nan in that case ...)
  if (m2 == 0.) {
    return -(Sigma_ferm_L(m2,m2,m2p,Qf,If));
  }
  else 
    return -(Sigma_ferm_L(m2,m2,m2p,Qf,If))
      -(m2*(dSigma_ferm_R(m2,m2,m2p,Qf,If)
	    +dSigma_ferm_L(m2,m2,m2p,Qf,If)
	    +2.*dSigma_ferm_S(m2,m2,m2p,Qf,If)));
}

// Fermion wavefunction counterterm, contribution proportional to PR = (1 + gamma5)/2
// Eq. (3.20)
DivArrC EW_One_Loop_Functions_Base::dZfermR(const double& m2, const double& m2p,
					    const double& Qf, const double& If) {
  // Do not evaluate derivatives of self energies if m2 = 0. (tend to be nan in that case ...)
  if (m2 == 0.) {
    return -(Sigma_ferm_R(m2,m2,m2p,Qf,If));
  }
  else 
    return -(Sigma_ferm_R(m2,m2,m2p,Qf,If))
      -m2*(dSigma_ferm_R(m2,m2,m2p,Qf,If)
	   +dSigma_ferm_L(m2,m2,m2p,Qf,If)
	   +2.*dSigma_ferm_S(m2,m2,m2p,Qf,If));
}

// Fermion mass counterterm
// Eq. (3.20)
DivArrC EW_One_Loop_Functions_Base::dm(const double& m2, const double& m2p,
				       const double& Qf, const double& If) {
  return sqrt(m2)/2.*(Sigma_ferm_L(m2,m2,m2p,Qf,If) 
		      + Sigma_ferm_R(m2,m2,m2p,Qf,If) 
		      + 2.*Sigma_ferm_S(m2,m2,m2p,Qf,If));
}

// Fermion self energies
// Contribution prop. to PR
// Eq. (B.6)
DivArrC EW_One_Loop_Functions_Base::Sigma_ferm_R(const double& p2,
						 const double& m2, const double& m2p,
						 const double& Qf, const double& If)
{
  Complex gplus(-m_sw/m_cw*Qf);
  return -1./4.*pow(Qf,2.)*(2.*B_1(p2,m2,0.,m_mu2)+1.*One)
    +(m_ew?(-1./4.*pow(gplus,2.)*(2.*B_1(p2,m2,muZ2,m_mu2)+1.*One)
	     -1./(16.*m_sw2)*m2/muW2*(B_1(p2,m2,muZ2,m_mu2)+B_1(p2,m2,muH2,m_mu2))
	     -1./(8.*m_sw2)*m2/muW2*B_1(p2,m2p,muW2,m_mu2)):Zero);
}

// Contribution prop. to PL
// Eq. (B.7)
DivArrC EW_One_Loop_Functions_Base::Sigma_ferm_L(const double& p2,
						 const double& m2, const double& m2p,
						 const double& Qf, const double& If)
{
  Complex gminus((If-m_sw2*Qf)/(m_sw*m_cw));
  return -1./4.*pow(Qf,2.)*(2.*B_1(p2,m2,0.,m_mu2)+1.*One)
    +(m_ew?(-1./4.*pow(gminus,2.)*(2.*B_1(p2,m2,muZ2,m_mu2)+1.*One)
	    -1./(16.*m_sw2)*m2/muW2*(B_1(p2,m2,muZ2,m_mu2)+B_1(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sw2)*((2.+m2p/muW2)*B_1(p2,m2p,muW2,m_mu2)+1.*One)):Zero);
}

// Scalar Contribution
// Eq. (B.8)
DivArrC EW_One_Loop_Functions_Base::Sigma_ferm_S(const double& p2,
						 const double& m2, const double& m2p,
						 const double& Qf, const double& If)
{
  Complex gplus(-m_sw/m_cw*Qf),gminus((If-m_sw2*Qf)/(m_sw*m_cw));
  return -1./4.*pow(Qf,2.)*(4.*B_0(p2,m2,0.,m_mu2)-2.*One)
    +(m_ew?(-1./4.*gplus*gminus*(4.*B_0(p2,m2,muZ2,m_mu2)-2.*One)
	     -1./(16.*m_sw2)*m2/muW2*(B_0(p2,m2,muZ2,m_mu2)-B_0(p2,m2,muH2,m_mu2))
	     -1./(8.*m_sw2)*m2p/muW2*B_0(p2,m2p,muW2,m_mu2)):Zero);
}

// Derivatives of fermion self energies, derived from the above
// Contribution prop. to PR
DivArrC EW_One_Loop_Functions_Base::dSigma_ferm_R(const double& p2,
						  const double& m2, const double& m2p,
						  const double& Qf, const double& If)
{
  Complex gplus(-m_sw/m_cw*Qf);
  return -1./4.*pow(Qf,2.)*2.*B_1p(p2,m2,0.,m_mu2)
    +(m_ew?(-1./4.*pow(gplus,2.)*2.*B_1p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sw2)*m2/muW2*(B_1p(p2,m2,muZ2,m_mu2)+B_1p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sw2)*m2/muW2*B_1p(p2,m2p,muW2,m_mu2)):Zero);
}

// Contribution prop. to PR
DivArrC EW_One_Loop_Functions_Base::dSigma_ferm_L(const double& p2,
					       const double& m2, const double& m2p,
					       const double& Qf, const double& If)
{
  Complex gminus((If-m_sw2*Qf)/(m_sw*m_cw));
  return -1./4.*pow(Qf,2.)*2.*B_1p(p2,m2,0.,m_mu2)
    +(m_ew?(-1./4.*pow(gminus,2.)*2.*B_1p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sw2)*m2/muW2*(B_1p(p2,m2,muZ2,m_mu2)+B_1p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sw2)*(2.+m2p/muW2)*B_1p(p2,m2p,muW2,m_mu2)):Zero);
}


// Scalar Contribution
DivArrC EW_One_Loop_Functions_Base::dSigma_ferm_S(const double& p2,
					       const double& m2, const double& m2p,
					       const double& Qf, const double& If)
{
  Complex gplus(-m_sw/m_cw*Qf),gminus((If-m_sw2*Qf)/(m_sw*m_cw));
  return -1./4.*pow(Qf,2.)*4.*B_0p(p2,m2,0.,m_mu2)
    +(m_ew?(-1./4.*gplus*gminus*4.*B_0p(p2,m2,muZ2,m_mu2)
	    -1./(16.*m_sw2)*m2/muW2*(B_0p(p2,m2,muZ2,m_mu2)-B_0p(p2,m2,muH2,m_mu2))
	    -1./(8.*m_sw2)*m2p/muW2*B_0p(p2,m2p,muW2,m_mu2)):Zero);
}


// W Boson self energy
// Eq. (B.4)
DivArrC EW_One_Loop_Functions_Base::Sigma_W(const double& p2)
{
  // Note: if p2 == 0 pieces 1/p2*(B0(p2,...)-B0(0,...)) become dB0(0,...)
  int kfl[] = {11,12,13,14,15,16}; // lepton IDs
  int kfq[] = {1,2,3,4,5,6}; // quark IDs
  // fermionic contributions
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  for (int i = 0; i < 3; ++i) {
    // leptonic contribution
    Flavour flav(Flavour(kfl[2*i]));
    double m2(pow(flav.Mass(),2.));
    fermDenner += 1./(12.*m_sw2)*(-(p2-m2/2.)*B_0(p2,0.,m2,m_mu2) + 1./3.*p2*One
				  +m2*B_0(0.,m2,m2,m_mu2) 
				  + (p2!=0.?pow(m2,2.)/(2.*p2)*(B_0(p2,0.,m2,m_mu2) - B_0(0.,0.,m2,m_mu2)):pow(m2,2.)/2.*B_0p(p2,0.,m2,m_mu2)));
  }
  for (int i = 0; i < 3; ++i) {
    // quark contribution - note that CKM matrix is the identity
    Flavour flavd(Flavour(kfq[2*i])),flavu(Flavour(kfq[2*i+1]));
    // make sure to use correct value of mt2 for comparison
    double m2d(pow(flavd.Mass(),2.)),m2u(((i==2)?mt2:pow(flavu.Mass(),2.)));
    fermDenner += 1./(4.*m_sw2)*(-(p2-(m2d+m2u)/2.)*B_0(p2,m2u,m2d,m_mu2)+1./3.*p2*One
				 +m2u*B_0(0.,m2u,m2u,m_mu2)+m2d*B_0(0.,m2d,m2d,m_mu2)
				 +(p2!=0.?pow(m2u-m2d,2.)/(2.*p2)*(B_0(p2,m2u,m2d,m_mu2)-B_0(0.,m2u,m2d,m_mu2)):pow(m2u-m2d,2.)/2.*B_0p(p2,m2u,m2d,m_mu2)));
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  // If only QED (m_ew == 0) only the terms from the photon coupling to the W contribute. 
  return -((m_ferm && m_ew)?fermDenner:Zero)
    - (m_nonferm?(1./6.*((2.*muW2+5.*p2)*B_0(p2,muW2,0.,m_mu2) - 2.*muW2*B_0(0.,muW2,muW2,m_mu2)
			 -(p2!=0.?
			   (pow(muW2,2.)/p2*(B_0(p2,muW2,0.,m_mu2)
					    -B_0(0.,muW2,0.,m_mu2)))
			   :(pow(muW2,2.)*B_0p(p2,muW2,0.,m_mu2)))
			 +1./3.*p2*One)
		  +(m_ew?(1./(48.*m_sw2)*(((40.*m_cw2-1.)*p2
					   +(16.*m_cw2+54.-10./m_cw2)*muW2)*B_0(p2,muW2,muZ2,m_mu2)
					  -(16.*m_cw2+2.)*(muW2*B_0(0.,muW2,muW2,m_mu2)
							   +muZ2*B_0(0.,muZ2,muZ2,m_mu2))
					  +(4.*m_cw2-1.)*2./3.*p2*One
					  -(8.*m_cw2+1.)*(p2!=0.?
							  (pow(muW2-muZ2,2.)/p2*(B_0(p2,muW2,muZ2,m_mu2)
									       -B_0(0.,muW2,muZ2,m_mu2))):
							  (pow(muW2-muZ2,2.)*B_0p(p2,muW2,muZ2,m_mu2)))
					  )
			  +1./(48.*m_sw2)*((2.*muH2-10.*muW2-p2)*B_0(p2,muW2,muH2,m_mu2)
					   -2.*muW2*B_0(0.,muW2,muW2,m_mu2)
					   -2.*muH2*B_0(0.,muH2,muH2,m_mu2)
					   -(p2!=0.?
					     (pow(muW2-muH2,2.)/p2*(B_0(p2,muW2,muH2,m_mu2)
								  -B_0(0.,muW2,muH2,m_mu2))):
					     pow(muW2-muH2,2.)*B_0p(p2,muW2,muH2,m_mu2))
					   -2./3.*p2*One)):Zero)):Zero);
}

// Derivative of W boson self energy, derived from above
DivArrC EW_One_Loop_Functions_Base::dSigma_W(const double& p2)
{
  int kfl[] = {11,12,13,14,15,16}; // lepton IDs
  int kfq[] = {1,2,3,4,5,6}; // quark IDs
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  // fermionic contributions
  for (int i = 0; i < 3; ++i) {
    // leptonic contributions
    Flavour flav(Flavour(kfl[2*i]));
    double m2(pow(flav.Mass(),2.));
    fermDenner += 1./(12.*m_sw2)*(-(p2-m2/2.)*B_0p(p2,0.,m2,m_mu2)
				  -B_0(p2,0.,m2,m_mu2)
				  + 1./3.*One
				  + pow(m2,2.)/(2.*p2)*B_0p(p2,0.,m2,m_mu2)
				  - pow(m2,2.)/(2.*pow(p2,2.))*(B_0(p2,0.,m2,m_mu2) - B_0(0.,0.,m2,m_mu2)));
  }
  for (int i = 0; i < 3; ++i) {
    // quark contributions - note that CKM is the identity matrix
    Flavour flavd(Flavour(kfq[2*i])),flavu(Flavour(kfq[2*i+1]));
    // make sure to use correct value of mt2 for comparison
    double m2d(pow(flavd.Mass(),2.)),m2u(((i==2)?mt2:pow(flavu.Mass(),2.)));
    fermDenner += 1./(4.*m_sw2)*(-(p2-(m2d+m2u)/2.)*B_0p(p2,m2u,m2d,m_mu2)
				 -B_0(p2,m2u,m2d,m_mu2)
				 +1./3.*One
				 +pow(m2u-m2d,2.)/(2.*p2)*B_0p(p2,m2u,m2d,m_mu2)
				 -pow(m2u-m2d,2.)/(2.*pow(p2,2.))*(B_0(p2,m2u,m2d,m_mu2)-B_0(0.,m2u,m2d,m_mu2)));
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  // If only QED (m_ew == 0) only the terms from the photon coupling to the W contribute. 
  return -((m_ferm && m_ew)?fermDenner:Zero)
    - (m_nonferm?(1./6.*((2.*muW2+5.*p2)*B_0p(p2,MW2,0.,m_mu2)
			 +5.*B_0(p2,muW2,0.,m_mu2)
			 -pow(muW2,2.)/p2*B_0p(p2,MW2,0.,m_mu2)
			 +pow(muW2/p2,2.)*(B_0(p2,muW2,0.,m_mu2)
					  -B_0(0.,muW2,0.,m_mu2))
			 +1./3.*One)
		  +(m_ew?(1./(48.*m_sw2)*(((40.*m_cw2-1.)*p2
				    +(16.*m_cw2+54.-10./m_cw2)*muW2)*B_0p(p2,muW2,muZ2,m_mu2)
				   +(40.*m_cw2-1.)*B_0(p2,muW2,muZ2,m_mu2)
				   +(4.*m_cw2-1.)*2./3.*One
				   -(8.*m_cw2+1.)*pow(muW2-muZ2,2.)/p2*B_0p(p2,muW2,muZ2,m_mu2)
				   +(8.*m_cw2+1.)*pow((muW2-muZ2)/p2,2.)*(B_0(p2,muW2,muZ2,m_mu2)
									-B_0(0.,muW2,muZ2,m_mu2)))
		  +1./(48.*m_sw2)*((2.*muH2-10.*muW2-p2)*B_0p(p2,muW2,muH2,m_mu2)
				   -B_0(p2,muW2,muH2,m_mu2)
				   -pow(muW2-muH2,2.)/p2*B_0p(p2,muW2,muH2,m_mu2)
				   +pow((muW2-muH2)/p2,2.)*(B_0(p2,muW2,muH2,m_mu2)
							  -B_0(0.,muW2,muH2,m_mu2))
				   -2./3.*One)):Zero)):Zero);
}



// Z Boson self energy
// Eq. (B.3)
DivArrC EW_One_Loop_Functions_Base::Sigma_ZZ(const double& p2)
{
  // No coupling if QED only
  if (m_ew == 0) return Zero; 

  int N_f(12);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  // Fermionic contributions
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:pow(flav.Mass(),2.)));
    Complex gplus(-m_sw/m_cw*flav.Charge()),gminus((flav.IsoWeak()-m_sw2*flav.Charge())/(m_sw*m_cw));
    fermDenner += 1./6.*((flav.IsQuark())?3.:1.)
      *((pow(gplus,2.)+pow(gminus,2.))*(-(p2+2.*m2)*B_0(p2,m2,m2,m_mu2) 
					+ 2.*m2*B_0(0.,m2,m2,m_mu2) 
					+ 1./3.*p2*One)
	+3./(4.*m_sw2*m_cw2)*m2*B_0(p2,m2,m2,m_mu2));
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return -(m_ferm?fermDenner:Zero) 
    -(m_nonferm?(1./(24.*m_sw2*m_cw2)*(((18.*m_cw4+2.*m_cw2-1./2.)*p2
					+(24.*m_cw4+16.*m_cw2-10.)*muW2)*B_0(p2,muW2,muW2,m_mu2)
				       -(24.*m_cw4-8.*m_cw2+2.)*muW2*B_0(0.,muW2,muW2,m_mu2)
				       +(4.*m_cw2-1.)*1./3.*p2*One)
		 +1./(48.*m_sw2*m_cw2)*((2.*muH2-10.*muZ2-p2)*B_0(p2,muZ2,muH2,m_mu2)
					-2.*muZ2*B_0(0.,muZ2,muZ2,m_mu2) - 2.*muH2*B_0(0.,muH2,muH2,m_mu2)
					-pow(muZ2-muH2,2.)/p2*(B_0(p2,muZ2,muH2,m_mu2)-B_0(0.,muZ2,muH2,m_mu2))
					-2./3.*p2*One)):Zero);
}

// Derivative of Z Boson self energy, derived from the above
DivArrC EW_One_Loop_Functions_Base::dSigma_ZZ(const double& p2)
{
  // No coupling if QED only
  if (m_ew == 0) return Zero; 

  int N_f(12);
  // Fermionic contributions
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; //fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:pow(flav.Mass(),2.)));
    Complex gplus(-m_sw/m_cw*flav.Charge()),gminus((flav.IsoWeak()-m_sw2*flav.Charge())/(m_sw*m_cw));
    fermDenner += 1./6.*((flav.IsQuark())?3.:1.)
      *((pow(gplus,2.)+pow(gminus,2.))*(-(p2+2.*m2)*B_0p(p2,m2,m2,m_mu2) 
					-B_0(p2,m2,m2,m_mu2)
					+ 1./3.*One)
	+3./(4.*m_sw2*m_cw2)*m2*B_0p(p2,m2,m2,m_mu2));
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return -(m_ferm?fermDenner:Zero)
    -(m_nonferm?(1./(24.*m_sw2*m_cw2)*(((18.*m_cw4+2.*m_cw2-1./2.)*p2
					+(24.*m_cw4+16.*m_cw2-10.)*muW2)*B_0p(p2,muW2,muW2,m_mu2)
				       +(18.*m_cw4+2.*m_cw2-1./2.)*B_0(p2,muW2,muW2,m_mu2)
				       +(4.*m_cw2-1.)*1./3.*One)
		 +1./(48.*m_sw2*m_cw2)*((2.*muH2-10.*muZ2-p2)*B_0p(p2,muZ2,muH2,m_mu2)
					-B_0(p2,muZ2,muH2,m_mu2)
					-pow(muZ2-muH2,2.)/p2*B_0p(p2,muZ2,muH2,m_mu2)
					+pow(muZ2-muH2,2.)/pow(p2,2.)*(B_0(p2,muZ2,muH2,m_mu2)-B_0(0.,muZ2,muH2,m_mu2))
					-2./3.*One)):Zero);
}

// Higgs Boson self energy
// Eq. (B.5)
DivArrC EW_One_Loop_Functions_Base::Sigma_H(const double& p2)
{
  // No coupling if QED only
  if (m_ew == 0) return Zero; 

  int N_f(12);
  // Fermionic contributions
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.); 
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:pow(flav.Mass(),2.)));
    fermDenner += 1./(8.*m_sw2*muW2)*(flav.IsQuark()?3.:1.)*m2*(2.*A_0(m2,m_mu2) + (4.*m2-p2)*B_0(p2,m2,m2,m_mu2));
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return -(m_ferm?fermDenner:Zero) 
    - (m_nonferm?(-1./(8.*m_sw2)*((6.*muW2-2.*p2+pow(muH2,2.)/(2.*muW2))*B_0(p2,muW2,muW2,m_mu2)
				  +(3.+muH2/(2.*muW2))*A_0(muW2,m_mu2) - 6.*muW2)
		  -1./(16.*m_sw2*m_cw2)*((6.*muZ2-2.*p2+pow(muH2,2.)/(2.*muZ2))*B_0(p2,muZ2,muZ2,m_mu2)
					 +(3.+muH2/(2.*muZ2))*A_0(muZ2,m_mu2) - 6.*muZ2)
		  -3./(32.*m_sw2)*(3.*pow(muH2,2.)/muW2*B_0(p2,muH2,muH2,m_mu2)
				   +muH2/muW2*A_0(muH2,m_mu2))):Zero);
}

// Derivative of Higgs Boson self energy, derived from the above
DivArrC EW_One_Loop_Functions_Base::dSigma_H(const double& p2)
{
  // No coupling if QED only
  if (m_ew == 0) return Zero; 

  // Fermionic contributions
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:pow(flav.Mass(),2.)));
    fermDenner += 1./(8.*m_sw2*muW2)*(flav.IsQuark()?3.:1.)*m2*((4.*m2-p2)*B_0p(p2,m2,m2,m_mu2)
							       -B_0(p2,m2,m2,m_mu2));
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return -(m_ferm?fermDenner:Zero) 
    - (m_nonferm?(-1./(8.*m_sw2)*((6.*muW2-2.*p2+pow(muH2,2.)/(2.*muW2))*B_0p(p2,muW2,muW2,m_mu2)
				  -2.*B_0(p2,muW2,muW2,m_mu2))
		  -1./(16.*m_sw2*m_cw2)*((6.*muZ2-2.*p2+pow(muH2,2.)/(2.*muZ2))*B_0p(p2,muZ2,muZ2,m_mu2)
					 -2.*B_0(p2,muZ2,muZ2,m_mu2))
		  -3./(32.*m_sw2)*3.*pow(muH2,2.)/muW2*B_0p(p2,muH2,muH2,m_mu2)):Zero);
}


// Derivative of Photon self energy, derived from Eq. (B.1)
// Photon self energy itself is never needed
DivArrC EW_One_Loop_Functions_Base::dSigma_GamGam(const double& p2)
{
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  // fermionic contributions
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:pow(flav.Mass(),2.)));
    if(m_ew) fermDenner += 1./3.*(flav.IsQuark()?3.:1.)*pow(flav.Charge(),2.)*(-B_0(p2,m2,m2,m_mu2)
									       -(p2+2.*m2)*B_0p(p2,m2,m2,m_mu2)
									       +1./3.*One);
    else fermDenner += 1./3.*(flav.IsQuark()?3.:1.)*pow(flav.Charge(),2.)*(-B_0(p2,m2,m2,m_mu2)
									   +1./3.*One);
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return -(m_ferm?fermDenner:Zero) 
    - ((m_nonferm && m_ew)?(1./4.*((3.*p2+4.*muW2)*B_0p(p2,muW2,muW2,m_mu2) 
				   +3.*B_0(p2,muW2,muW2,m_mu2))):Zero);
}

// Derivative of photon self energy at p2 = 0., including some terms to resum leading logarithms
// from fermion contributions in alpha(0) scheme - see OpenLoops implementation
DivArrC EW_One_Loop_Functions_Base::dSigma_GamGam0()
{
  DivArrC PiAALightZ = Zero;
  // Contributions to photon self energy at p2 = 0.
  // Assuming leptons, bottom and top are massive
  double mb2 = pow(Flavour(5).Mass(),2.);
  double me2 = pow(Flavour(11).Mass(),2.);
  double mm2 = pow(Flavour(13).Mass(),2.);
  double mtau2 = pow(Flavour(15).Mass(),2.);
  PiAALightZ += 1./3.*
    (2.*me2*B_0(0.,me2,me2,m_mu2) - (2.*me2+muZ2)*B_0(MZ2,me2,me2,m_mu2)
     +2.*mm2*B_0(0.,mm2,mm2,m_mu2) - (2.*mm2+muZ2)*B_0(MZ2,mm2,mm2,m_mu2)
     +2.*mtau2*B_0(0.,mtau2,mtau2,m_mu2) - (2.*mtau2+muZ2)*B_0(MZ2,mtau2,mtau2,m_mu2)
     +muZ2*(1.*One+11./27.*One-11./3.*B_0(MZ2,0.,0.,m_mu2)));
  DivArrC dSigmaTop = 4./9.*(1./3.*One-B_0(0.,mut2,mut2,m_mu2)+2.*mt2*B_0p(0.,mut2,mut2,m_mu2));
  double alphaQED_0 = 1./137.03599976;
  double alphaQED_MZ = 1./128.;
  double dAlphaQED_MZ = M_PI/Photons::s_alpha*(1.-alphaQED_0/alphaQED_MZ);
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return - (m_ferm?(PiAALightZ/muZ2 + dSigmaTop):Zero) 
    - ((m_nonferm && m_ew)?(1./4.*(4.*muW2*B_0p(0.,muW2,muW2,m_mu2) 
				   +3.*B_0(0.,muW2,muW2,m_mu2))):Zero)
    +dAlphaQED_MZ*One;
}

// Z-Photon transition
// Eq. (B.2)
DivArrC EW_One_Loop_Functions_Base::Sigma_ZA(const double& p2)
{
  // No coupling to Z if QED only
  if (m_ew == 0) return Zero;

  // fermionic contributions
  int N_f(12);
  DivArrC fermDenner(0.,0.,0.,0.,0.,0.);
  int kf[] = {1,2,3,4,5,6,11,12,13,14,15,16}; // fermion IDs
  for (int i = 0; i < N_f; ++i) {
    Flavour flav = Flavour(kf[i]);
    // make sure to use correct value of mt2 for comparison
    double m2((i==5?mt2:pow(flav.Mass(),2.))),ch(flav.Charge());
    Complex gplus(-m_sw/m_cw*flav.Charge()),gminus((flav.IsoWeak()-m_sw2*flav.Charge())/(m_sw*m_cw));
    fermDenner += 1./6.*(flav.IsQuark()?3.:1.)*(-flav.Charge())*(gplus+gminus)*(-(p2+2.*m2)*B_0(p2,m2,m2,m_mu2)
								      +2.*m2*B_0(0.,m2,m2,m_mu2)
										+1./3.*p2*One);
  }
  // Total contribution including bosonic terms. 
  // Return only either part if fermionic or bosonic terms only.
  return -(m_ferm?fermDenner:Zero)
    + ((m_nonferm && m_ew)?(1./(12.*m_sw*m_cw)*(((9.*m_cw2+1./2.)*p2+(12.*m_cw2+4.)*muW2)*B_0(p2,muW2,muW2,m_mu2)
				      -(12.*m_cw2-2.)*muW2*B_0(0.,muW2,muW2,m_mu2)
				      +1./3.*p2*One)):Zero);
}

// Wavefunction Counterterms
// All the following counterterms are found in Eq. (3.19)
DivArrC EW_One_Loop_Functions_Base::dZW()
{
  return -dSigma_W(MW2);
}

DivArrC EW_One_Loop_Functions_Base::dZZZ()
{
  return -dSigma_ZZ(MZ2);
}

DivArrC EW_One_Loop_Functions_Base::dZAA()
{
  return -dSigma_GamGam(MZ2);
}

DivArrC EW_One_Loop_Functions_Base::dZAZ()
{
  return -2.*Sigma_ZA(MZ2)/MZ2 + (muZ2/MZ2 - 1.)*dZZA();
}

DivArrC EW_One_Loop_Functions_Base::dZZA()
{
  return 2.*Sigma_ZA(0.)/muZ2;
}

DivArrC EW_One_Loop_Functions_Base::dZH()
{
  return -dSigma_H(MH2);
}

// Mass Counterterms
DivArrC EW_One_Loop_Functions_Base::dMW2()
{
  return Sigma_W(MW2) + (muW2 - MW2)*dSigma_W(MW2);
}

DivArrC EW_One_Loop_Functions_Base::dMZ2()
{
  return Sigma_ZZ(MZ2) + (muZ2-MZ2)*dSigma_ZZ(MZ2);
}

DivArrC EW_One_Loop_Functions_Base::dMH2()
{
  return Sigma_H(MH2) + (muH2-MH2)*dSigma_H(MH2);
}

// Counterterm for cos(thetaW)
// Eq. (3.35)
DivArrC EW_One_Loop_Functions_Base::dcw()
{
  return 1./2.*(dMW2()/muW2-dMZ2()/muZ2);
}

// Charge renormalization counterterm
// Eq. (3.32)
DivArrC EW_One_Loop_Functions_Base::dZe()
{
  switch (Photons::s_ew_scheme) {
  case 1: {
    // alpha(MZ) scheme
    return 1./2.*dSigma_GamGam(MZ2) - m_sw/m_cw*Sigma_ZA(0.)/muZ2;
    break;
  }
  case 2: {
    // alpha(0) scheme
    return 1./2.*dSigma_GamGam0() - m_sw/m_cw*Sigma_ZA(0.)/muZ2;
    break;
  }
  case 3: {
    // Gmu scheme
    return -m_cw2/m_sw2*dcw() + 1./2.*(Sigma_W(MW2)-Sigma_W(0.))/muW2 
      - 1./(m_cw*m_sw)*Sigma_ZA(0.)/muZ2 - One*1./8.*1./m_sw2*(6.+(7.-4.*m_sw2)/(2.*m_sw2)*log(m_cw2));
    break;
  }
  default: { 
    // by default, use alpha(0) scheme
    return 1./2.*dSigma_GamGam0() - m_sw/m_cw*Sigma_ZA(0.)/muZ2;
    break;
  }
  }
}

// print out finite parts of renormalization constants
// factor 4 needed to compare to OpenLoops
void EW_One_Loop_Functions_Base::Print_Ren_Constants_Finite() {
  msg_Out() << "Finite\n" <<
    "debug_norm " << m_alpha/(4.*M_PI) << "\n" <<
    "WF_V\n" <<
    "dZAAEW " << 4.*dZAA().Finite() << "\n" <<
    "dZAZEW " << 4.*dZAZ().Finite() << "\n" <<
    "dZZAEW " << 4.*dZZA().Finite() << "\n" <<
    "dZZZEW " << 4.*dZZZ().Finite() << "\n" << 
    "dZWEW " << 4.*dZW().Finite() << "\n" << 
    "dZHEW " << 4.*dZH().Finite() << "\n" << 
    "M_V\n" << 
    "MW2 " << 4.*dMW2().Finite() << "\n" <<
    "MZ2 " << 4.*dMZ2().Finite() << "\n" <<
    "MH2 " << 4.*dMH2().Finite() << "\n" <<
    "C\n" <<
    "dZe " << 4.*dZe().Finite() << "\n" <<
    "dcw " << 4.*dcw().Finite() << "\n\n";
}

// print out UV parts of renormalization constants
// factor 4 needed to compare to OpenLoops
void EW_One_Loop_Functions_Base::Print_Ren_Constants_UV() {
  msg_Out() << "UV\n" <<
    "WF_V\n" <<
    "WF_V\n" <<
    "dZAAEW " << 4.*dZAA().UV() << "\n" <<
    "dZAZEW " << 4.*dZAZ().UV() << "\n" <<
    "dZZAEW " << 4.*dZZA().UV() << "\n" <<
    "dZZZEW " << 4.*dZZZ().UV() << "\n" << 
    "dZWEW " << 4.*dZW().UV() << "\n" << 
    "dZHEW " << 4.*dZH().UV() << "\n" << 
    "M_V\n" << 
    "MW2 " << 4.*dMW2().UV() << "\n" <<
    "MZ2 " << 4.*dMZ2().UV() << "\n" <<
    "MH2 " << 4.*dMH2().UV() << "\n" <<
    "C\n" <<
    "dZe " << 4.*dZe().UV() << "\n" <<
    "dcw " << 4.*dcw().UV() << "\n\n";
}

// print out IR parts of renormalization constants
// factor 4 needed to compare to OpenLoops
void EW_One_Loop_Functions_Base::Print_Ren_Constants_IR() {
  msg_Out() << "IR\n" <<
    "WF_V\n" <<
    "WF_V\n" <<
    "dZAAEW " << 4.*dZAA().IR() << "\n" <<
    "dZAZEW " << 4.*dZAZ().IR() << "\n" <<
    "dZZAEW " << 4.*dZZA().IR() << "\n" <<
    "dZZZEW " << 4.*dZZZ().IR() << "\n" << 
    "dZWEW " << 4.*dZW().IR() << "\n" << 
    "dZHEW " << 4.*dZH().IR() << "\n" << 
    "M_V\n" << 
    "MW2 " << 4.*dMW2().IR() << "\n" <<
    "MZ2 " << 4.*dMZ2().IR() << "\n" <<
    "MH2 " << 4.*dMH2().IR() << "\n" <<
    "C\n" <<
    "dZe " << 4.*dZe().IR() << "\n" <<
    "dcw " << 4.*dcw().IR() << "\n\n";
}
