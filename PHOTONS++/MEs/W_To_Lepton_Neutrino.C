#include "PHOTONS++/MEs/W_To_Lepton_Neutrino.H"
#include "PHOTONS++/MEs/W_Decay_EW_One_Loop_Functions.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Math/Vector.H"
#include "ATOOLS/Math/Tensor.H"
#include "ATOOLS/Math/Tensor_Build.H"
#include "ATOOLS/Math/Tensor_Contractions.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "ATOOLS/Org/Data_Reader.H"
#include "ATOOLS/Org/Message.H"
#include "ATOOLS/Math/Random.H"
#include "PHOTONS++/Main/Photons.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)  Master_Triangle(A,B,C,D,E,F,M)

//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////
////                                                                  ////
////     CAUTION :                                                    ////
////                                                                  ////
////     pvv_zero contains m_chargedoutparticles at position 2        ////
////     --> l sits at position 0                                     ////
////                                                                  ////
////     pvv_one contains m_newdipole at position 2                   ////
////     --> W sits by default at m_newdipole.at(0)                   ////
////     --> l sits on postion 1                                      ////
////                                                                  ////
//////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////

W_To_Lepton_Neutrino::W_To_Lepton_Neutrino
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FI(pvv) {
  if (m_ew == 0) {
    msg_Out() << METHOD << ": Pure QED corrections not gauge-invariant in this case." 
	      << "Setting m_ew to 1 and calculate full EW corrections.\n";
    m_ew = 1;
  }
  m_no_weight = 0;
  m_name = "W_To_Lepton_Neutrino";
  m_flavs[0]  = pvv[0][0]->Flav();
  m_masses[0] = pvv[0][0]->FinalMass();
  m_flavs[1]  = pvv[2][0]->Flav();
  m_masses[1] = pvv[2][0]->FinalMass();
  m_flavs[2]  = pvv[3][0]->Flav();
  m_masses[2] = pvv[3][0]->FinalMass();
  for (unsigned int i=3; i<9; i++) {
    m_flavs[i]  = Flavour(kf_photon);
    m_masses[i] = 0.;
  }
  m_cL = Complex(1.,0.);
  m_cR = Complex(0.,0.);
}

W_To_Lepton_Neutrino::~W_To_Lepton_Neutrino() {
}

void W_To_Lepton_Neutrino::BoostOriginalPVVToMultipoleCMS() {
  // m_pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole.at(0) into +z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum += m_olddipole[i]->Momentum();
  }
  Vec4D p1 = m_olddipole[0]->Momentum();
  p_boost = new Poincare(sum);
  p_boost->Boost(p1);
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,-1.)); // in Dipole_FI W is rotated into -z-direction!!!
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D vec = m_olddipole[i]->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_olddipole[i]->SetMomentum(vec);
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D vec = m_oldspectator[i]->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_oldspectator[i]->SetMomentum(vec);
  }
}

void W_To_Lepton_Neutrino::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  // m_moms0 - no photon
  Poincare boost(m_pvv_zero[0][0]->Momentum());
  Vec4D vec;
  vec = m_pvv_zero[0][0]->Momentum();
  boost.Boost(vec);
  m_moms0[0] = vec;
  vec = m_pvv_zero[2][0]->Momentum();
  boost.Boost(vec);
  m_moms0[1] = vec;
  vec = m_pvv_zero[3][0]->Momentum();
  boost.Boost(vec);
  m_moms0[2] = vec;
  vec = m_pvv_zero[0][0]->Momentum();

  // Store fully rescaled final state
  m_momsFull[0] = pvv_one[2][0]->Momentum();
  m_momsFull[1] = pvv_one[2][1]->Momentum();
  m_momsFull[2] = pvv_one[3][0]->Momentum();
  for (size_t i = 0; i < pvv_one[4].size(); i++) {
    m_momsFull[i+3] = pvv_one[4][i]->Momentum();
  }


  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated
  Dipole_FI::DefineDipole();

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[2][0]->Momentum();
    m_moms1[0][1] = pvv_one[2][1]->Momentum();
    m_moms1[0][2] = pvv_one[3][0]->Momentum();
    m_moms1[0][3] = pvv_one[4][0]->Momentum();
  }
  else if (pvv_one[4].size() > 1) {
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<pvv_one[4].size(); i++) {
      m_softphotons.push_back(pvv_one[4][i]);
      m_K = CalculateMomentumSum(m_softphotons);
      DetermineQAndKappa();
      CorrectMomenta();
      if (m_u == -1.) {
	// if any u = -1. is found, momenta are not corrected properly since root could not be found
	// return no weights in this case!
	m_no_weight = 1;
	continue;
      }
      m_moms1[i][0] = m_newdipole[0]->Momentum();
      m_moms1[i][1] = m_newdipole[1]->Momentum();
      m_moms1[i][2] = m_newspectator[0]->Momentum();
      m_moms1[i][3] = m_softphotons[0]->Momentum();
      m_softphotons.clear();
    }
  }
}

double W_To_Lepton_Neutrino::SmodFull(unsigned int kk) {
  // Smod calculated using fully dressed momenta 
  // Used to cancel the Smod in the dipole weight

  // Index kk denotes which photon is used in calculation. Note that m_momsFull contains all photons
  // starting at position 3, hence index of photon kk is 3+kk.
  m_moms = m_momsFull;
  Vec4D k   = m_moms[3+kk];
  Vec4D pi  = m_moms[0];
  Vec4D pj  = m_moms[1];
  double Zi = -1.;
  double Zj = -1.;
  int    ti = -1;
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double W_To_Lepton_Neutrino::Smod(unsigned int kk) {
  // Smod calculated with momenta dressed with one additional hard photon
  // Used in the subtraction in single-real matrix elements

  // Index kk denotes which photon is taken to be hard - determines which set of momenta m_moms1[kk]
  // needs to be called. Note that the photon is always at position 3 in the sets.
  m_moms = m_moms1[kk];
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms[0];
  Vec4D pj  = m_moms[1];
  double Zi = -1.;
  double Zj = -1.;
  int    ti = -1;
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double W_To_Lepton_Neutrino::Smod(unsigned int a, unsigned int b, unsigned int c) {
  // Smod calculated with momenta dressed with two additional hard photons
  // Used in subtraction in double-real matrix elements

  // Index a denotes position of first momentum, index b position of second photon. Used to get
  // correct set of dressed momenta m_moms2[a][b]
  // Index c denotes whether the photon in the calculation is the first or second of the two 
  // (c is either 0 or 1)
  m_moms = m_moms2[a][b];
  Vec4D k   = m_moms[3+c];
  Vec4D pi  = m_moms[0];
  Vec4D pj  = m_moms[1];
  double Zi = -1.;
  double Zj = -1.;
  int    ti = -1;
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_0_0() {
  m_moms = m_moms0;
  Vec4C epsW = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  // one diagram:
  // M_0 = ie/(sqrt(2)sW) * ubar(l)gamma_rho P_L v(nu) eps_rho^W
  if (m_flavs[1].IsAnti()) {
    return m_i*m_e/(m_sqrt2*m_sW)*XYZ.X(1,m_spins[1],epsW,2,m_spins[2],m_cR,m_cL);
  }
  else {
    return m_i*m_e/(m_sqrt2*m_sW)*XYZ.X(2,m_spins[2],epsW,1,m_spins[1],m_cR,m_cL);    
  }
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_0_1(const int& ewmode) {
  //return 0.;
  m_moms = m_moms0;
  // p1, unprimed quantities correspond to particle
  // p2, primed quantities correspond to antiparticle - note that the properties of the fermion are used in the EW one-loop functions rather than the antifermion (e. g. W->nu l+ -> Qfp = -1)!
  Vec4D p1,p2;
  double m2, m2p, If, Ifp, Qf, Qfp;

  // Setting up proper values for EW parameters:
  // Qf = 1; corresponding to If = 1/2 and m2 = ml2. Ifp = -1/2 and Qfp = 0, mfp = 0.
  // Need to offset the extra minus signs in .Charge() and .IsoWeak() for anti-particles
  if (!m_flavs[1].IsAnti()) {
    p1 = m_moms[1];
    p2 = m_moms[2];
    m2 = pow(m_masses[1],2.);
    m2p = pow(m_masses[2],2.);
    Qf = -m_flavs[1].Charge();
    Qfp = m_flavs[2].Charge();
    If = -m_flavs[1].IsoWeak();
    Ifp = m_flavs[2].IsoWeak();
  }
  else if (m_flavs[1].IsAnti()) {
    p1 = m_moms[1];
    p2 = m_moms[2];
    m2 = pow(m_masses[1],2.);
    m2p = pow(m_masses[2],2.);
    Qf = m_flavs[1].Charge();
    Qfp = -m_flavs[2].Charge();
    If = m_flavs[1].IsoWeak();
    Ifp = -m_flavs[2].IsoWeak();
  }  
  Vec4C epsV = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  double s((m_moms[1]+m_moms[2]).Abs2());
  double p1p2(m_moms[1]*m_moms[2]);
  double mu2(s);

  // Set up EW corrections as corrections to cR,cL
  Complex gfm((If-m_sW2*Qf)/(m_sW*m_cW)),gfpm((Ifp-m_sW2*Qfp)/(m_sW*m_cW));
  W_Decay_EW_One_Loop_Functions* OL = new W_Decay_EW_One_Loop_Functions(s, m2, m2p, Qf, Qfp, If, Ifp, mu2, m_ew);
  Complex term(0.,0.);
  Complex pref(Complex(0.,1.)*m_e/(m_sqrt2*m_sW));
  // pieces proportional to X(1,epsV,2)
  DivArrC cRDivArrX = m_alpha/M_PI*(pref*(OL->CT_R()
					  +1./4.*Qf*Qfp*OL->FQED(1,1)));
  DivArrC cLDivArrX = m_alpha/M_PI*(pref*(
					  OL->CT_L()
					  +1./4.*Qf*Qfp*OL->FQED(1,0)
					  +(Qf==-1?-1.:1.)*1./2.*OL->FAn()
					  +((m_ew>=1)?(1./(4.*m_sW2*m_cW2)*(Ifp-m_sW2*Qfp)*(If-m_sW2*Qf)*OL->FZa()
					  +(Qf==-1?1.:-1.)*1./m_sW2*(Ifp-m_sW2*Qfp-If+m_sW2*Qf)*OL->FZn()):DivArrC(0.,0.,0.,0.,0.,0.))
					  ));
  // pieces proportional to P_R/L (multiplied by eps*momenta) - no counterterm contribution
  // should vanish for massless case
  DivArrC cRDivArrY = m_alpha/M_PI*(pref*1./4.*Qf*Qfp*Complex(0.,-1.)*((epsV*(p1-p2))*OL->FQED(2,1)
								       +(epsV*(p1+p2))*OL->FQED(3,1)));
  DivArrC cLDivArrY = m_alpha/M_PI*(pref*1./4.*Qf*Qfp*Complex(0.,-1.)*((epsV*(p1-p2))*OL->FQED(2,0)
								       +(epsV*(p1+p2))*OL->FQED(3,0)));
  delete OL;
  DivArrC B(0.,0.,0.,0.,0.,0.);
  if (m_flavs[1].IsAnti()) {
    if (Qf != 0.) {
      if (m2 != 0.) { B += m_alpha/M_PI*pref*(-0.5*(s-m2p+m2)*C_0(m2,m2p,s,0.,m2,muW2,mu2)
					      +0.5*s*C_0(s,s,0.,0.,s,s,mu2)
					      +0.5*m2*C_0(m2,m2,0.,0.,m2,m2,mu2)
					      +0.25*B_0(m2p,s,m2,mu2)
					      -0.125*B_0(0.,m2,m2,mu2)
					      -0.125*B_0(0.,muW2,muW2,mu2))*
	  XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
      }
      else if (m2 == 0.) {
	B += m_alpha/M_PI*pref*(-0.5*(s-m2p+m2)*C_0(m2,m2p,s,0.,m2,muW2,mu2)
				+0.5*s*C_0(s,s,0.,0.,s,s,mu2)
				+0.25*DivArrC(0.,1.,0.,0.,0.,0.)
				+0.25*B_0(m2p,s,m2,mu2)
				-0.125*B_0(0.,m2,m2,mu2)
				-0.125*B_0(0.,s,s,mu2))*
	  XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
      }
    }
    else {
      if (m2p != 0.) { B += m_alpha/M_PI*pref*(-0.5*(s-m2+m2p)*C_0(m2p,m2,s,0.,m2p,muW2,mu2)
					       +0.5*s*C_0(s,s,0.,0.,s,s,mu2)
					       +0.5*m2p*C_0(m2p,m2p,0.,0.,m2p,m2p,mu2)
					       +0.25*B_0(m2,s,m2p,mu2)
					       -0.125*B_0(0.,m2p,m2p,mu2)
					       -0.125*B_0(0.,muW2,muW2,mu2))*
	  XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
      }
      else if (m2p == 0.) {
	B += m_alpha/M_PI*pref*(-0.5*(s-m2+m2p)*C_0(m2p,m2,s,0.,m2p,muW2,mu2)
				+0.5*s*C_0(s,s,0.,0.,s,s,mu2)
				+0.25*DivArrC(0.,1.,0.,0.,0.,0.)
				+0.25*B_0(m2,s,m2p,mu2)
				-0.125*B_0(0.,m2p,m2p,mu2)
				-0.125*B_0(0.,muW2,muW2,mu2))*
	  XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
      }
    }
    // Assemble finite parts of results
    term += XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.Finite(),cLDivArrX.Finite())
      +XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Finite(),cLDivArrY.Finite()) + B.Finite();
    
    if (msg_LevelIsDebugging()) {
      m_res = DivArrC(0.,0.,0.,0.,0.,0.);
      m_res = 
	DivArrC(XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.UV(),cLDivArrX.UV()),
		XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.IR(),cLDivArrX.IR()),
		XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.IR2(),cLDivArrX.IR2()),
		XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.Finite(),cLDivArrX.Finite()),
		XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.Epsilon(),cLDivArrX.Epsilon()),
		XYZ.X(1,m_spins[1],epsV,2,m_spins[2],cRDivArrX.Epsilon2(),cLDivArrX.Epsilon2()))
	+DivArrC(XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.UV(),cLDivArrY.UV()),
		 XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.IR(),cLDivArrY.IR()),
		 XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.IR2(),cLDivArrY.IR2()),
		 XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Finite(),cLDivArrY.Finite()),
		 XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Epsilon(),cLDivArrY.Epsilon()),
		 XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArrY.Epsilon2(),cLDivArrY.Epsilon2()))
	+B
	;
    }
  }
  else {
    if (Qf != 0.) {
      if (m2 != 0.) { B += m_alpha/M_PI*pref*(-0.5*(s-m2p+m2)*C_0(m2,m2p,s,0.,m2,muW2,mu2)
					      +0.5*s*C_0(s,s,0.,0.,s,s,mu2)
					      +0.5*m2*C_0(m2,m2,0.,0.,m2,m2,mu2)
					      +0.25*B_0(m2p,s,m2,mu2)
					      -0.125*B_0(0.,m2,m2,mu2)
					      -0.125*B_0(0.,muW2,muW2,mu2))*
	  XYZ.X(2,m_spins[2],epsV,1,m_spins[1],m_cR,m_cL);
      }
      else if (m2 == 0.) {
	B += m_alpha/M_PI*pref*(-0.5*(s-m2p+m2)*C_0(m2,m2p,s,0.,m2,muW2,mu2)
				+0.5*s*C_0(s,s,0.,0.,s,s,mu2)
				+0.25*DivArrC(0.,1.,0.,0.,0.,0.)
				+0.25*B_0(m2p,s,m2,mu2)
				-0.125*B_0(0.,m2,m2,mu2)
				-0.125*B_0(0.,s,s,mu2))*
	  XYZ.X(2,m_spins[2],epsV,1,m_spins[1],m_cR,m_cL);
      }
    }
    else {
      if (m2p != 0.) { B += m_alpha/M_PI*pref*(-0.5*(s-m2+m2p)*C_0(m2p,m2,s,0.,m2p,muW2,mu2)
					       +0.5*s*C_0(s,s,0.,0.,s,s,mu2)
					       +0.5*m2p*C_0(m2p,m2p,0.,0.,m2p,m2p,mu2)
					       +0.25*B_0(m2,s,m2p,mu2)
					       -0.125*B_0(0.,m2p,m2p,mu2)
					       -0.125*B_0(0.,muW2,muW2,mu2))*
	  XYZ.X(2,m_spins[2],epsV,1,m_spins[1],m_cR,m_cL);
      }
      else if (m2p == 0.) {
	B += m_alpha/M_PI*pref*(-0.5*(s-m2+m2p)*C_0(m2p,m2,s,0.,m2p,muW2,mu2)
				+0.5*s*C_0(s,s,0.,0.,s,s,mu2)
				+0.25*DivArrC(0.,1.,0.,0.,0.,0.)
				+0.25*B_0(m2,s,m2p,mu2)
				-0.125*B_0(0.,m2p,m2p,mu2)
				-0.125*B_0(0.,muW2,muW2,mu2))*
	  XYZ.X(2,m_spins[2],epsV,1,m_spins[1],m_cR,m_cL);
      }
    }
    term += XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.Finite(),cLDivArrX.Finite())
      +XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.Finite(),cLDivArrY.Finite()) + B.Finite();
    
    if (msg_LevelIsDebugging()) {
      m_res = DivArrC(0.,0.,0.,0.,0.,0.);
      m_res = 
	DivArrC(XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.UV(),cLDivArrX.UV()),
		XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.IR(),cLDivArrX.IR()),
		XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.IR2(),cLDivArrX.IR2()),
		XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.Finite(),cLDivArrX.Finite()),
		XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.Epsilon(),cLDivArrX.Epsilon()),
		XYZ.X(2,m_spins[2],epsV,1,m_spins[1],cRDivArrX.Epsilon2(),cLDivArrX.Epsilon2()))
	+DivArrC(XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.UV(),cLDivArrY.UV()),
		 XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.IR(),cLDivArrY.IR()),
		 XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.IR2(),cLDivArrY.IR2()),
		 XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.Finite(),cLDivArrY.Finite()),
		 XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.Epsilon(),cLDivArrY.Epsilon()),
		 XYZ.Y(2,m_spins[2],1,m_spins[1],cRDivArrY.Epsilon2(),cLDivArrY.Epsilon2()))
	+B
	;
    }
    
  }
  return term;
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_0_2() {
  return 0.;
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_1_05(unsigned int i, const double& xi) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsW   = Polarization_Vector(m_moms[0])[m_spins[0]];
  Vec4C epsP   = conj(Polarization_Vector(m_moms[3])[m_spins[3]]);
  //Vec4C epsP   = m_moms[3]; // To check Ward identity
  Vec4D q      = m_moms[1]+m_moms[3];       // fermion propagator momenta
  double q2    = q.Abs2();
  Vec4D Q      = m_moms[0]-m_moms[3];       // boson propagator momenta
  double Q2    = Q.Abs2();
  double m     = m_masses[1];               // fermion mass/propagator pole
  double m2    = sqr(m);
  double M     = m_masses[0];               // boson mass/propagator pole
  double M2    = sqr(M);
  m_moms[4]    = m_moms[5] = q;             // enter those into m_moms
  m_flavs[4]   = m_flavs[1];                // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[1].Bar();
  XYZFunc XYZ(6,m_moms,m_flavs,false);
  // two diagrams
  // M_1 = -ie^2/(sqrt(2)sW) * 1/((pl+k)^2-m^2)
  //       * ubar(l)gamma^mu(-pl-k+m)gamma^nu P_L v(nu) eps_nu^W eps_mu^y*
  // M_2 = ie^2/(sqrt(2)sW) * 1/(pW-k)^2-M^2)
  //       * ubar(l)gamma_rho P_L v(nu)
  //       * [-2g^{rho,nu}pW^mu + g^{rho,mu}(pW-2k)^nu
  //          + 2g^{nu,mu}k^rho + 1/pW^2(pW-k)^rho pW^nu pW^mu]
  //       * eps_nu^W eps_mu^y*
  // Note that we need to use pW^2 instead of MW^2 in the propagator since W may be off-shell!
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Lorentz_Ten3C ten31,ten32,ten33,ten34,ten35;
  if (m_flavs[1].IsAnti()) {
    for (unsigned int s=0; s<=1; s++) {
      r1 += XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)
	*XYZ.X(4,s,epsW,2,m_spins[2],m_cR,m_cL);
      r2 += XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)
	*XYZ.X(5,s,epsW,2,m_spins[2],m_cR,m_cL);
    }
    Vec4D p = m_moms[0];
    Vec4D k = m_moms[3];
    // index ordering rho(1),nu(2),mu(3)
    // -2g^{rho,nu}pW^mu
    ten31 = BuildTensor(MetricTensor(),-2.*p);
    // g^{rho,mu}(pW-2.*k)^nu
    ten32 = BuildTensor(MetricTensor(),p-2.*k).Transpose(2,3);
    // 2g^{nu,mu}k^rho
    ten33 = BuildTensor(MetricTensor(),2.*k).Transpose(1,3);
    // 1/pW^2(pW-k)^rho pW^nu pW^mu
    ten34 = 1./(p*p)*BuildTensor(p-k,p,p);
    
    Lorentz_Ten3C ten = ten31+ten32+ten33+ten34;
    // v^\sigma = L^\sigma\mu\lambda epsW_\mu epsP_\lambda
    Vec4C v3 = Contraction(Contraction(ten,3,epsP),2,epsW);
    r3 = XYZ.X(1,m_spins[1],v3,2,m_spins[2],m_cR,m_cL);
    r1 *= (1.+m/sqrt(q2))/(2.*(q2-m2));
    r2 *= (1.-m/sqrt(q2))/(2.*(q2-m2));
    r3 *= -1./(Q2-p*p);
  }
  else {
    for (unsigned int s=0; s<=1; s++) {
      r1 += XYZ.X(4,s,epsP,1,m_spins[1],1.,1.)
	*XYZ.X(2,m_spins[2],epsW,4,s,m_cR,m_cL);
      r2 += XYZ.X(5,s,epsP,1,m_spins[1],1.,1.)
	*XYZ.X(2,m_spins[2],epsW,5,s,m_cR,m_cL);
    }
    Vec4D p = m_moms[0];
    Vec4D k = m_moms[3];
    // index ordering rho(1),nu(2),mu(3)
    // -2g^{rho,nu}pW^mu
    ten31 = BuildTensor(MetricTensor(),-2.*p);
    // g^{rho,mu}(pW-2.*k)^nu
    ten32 = BuildTensor(MetricTensor(),p-2.*k).Transpose(2,3);
    // 2g^{nu,mu}k^rho
    ten33 = BuildTensor(MetricTensor(),2.*k).Transpose(1,3);
    // 1/pW^2(pW-k)^rho pW^nu pW^mu
    ten34 = 1./(p*p)*BuildTensor(p-k,p,p);
    
    Lorentz_Ten3C ten = ten31+ten32+ten33+ten34;
    // v^\sigma = L^\sigma\mu\lambda epsW_\mu epsP_\lambda
    Vec4C v3 = Contraction(Contraction(ten,3,epsP),2,epsW);
    r3 = XYZ.X(2,m_spins[2],v3,1,m_spins[1],m_cR,m_cL);
    r1 *= (1.+m/sqrt(q2))/(2.*(q2-m2));
    r2 *= (1.-m/sqrt(q2))/(2.*(q2-m2));
    r3 *= -1./(Q2-p*p);
  }

  // erase intermediate entries from m_flavs
  m_flavs[4] = m_flavs[5] = Flavour(kf_none);
  // Note there was a factor of 2 wrong here - only applies to r1 and r2, since they involve spin sums
  return -(m_i*m_e*m_e)/(m_sqrt2*m_sW)*(r1+r2+r3);
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_1_15(unsigned int i, const double& xi) {
  return 0.;
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  return 0.;
}

double W_To_Lepton_Neutrino::GetBeta_0_0() {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin W
        m_spins[0] = k;
        m_spins[1] = j;
        m_spins[2] = i;
        Complex M_0_0 = InfraredSubtractedME_0_0();
        sum = sum + (M_0_0*conj(M_0_0)).real();
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  return sum;
}

double W_To_Lepton_Neutrino::GetBeta_0_1(const int& ewmode) {
  // limit mW > ml
  if (m_limit == 1) {
    return m_alpha/M_PI * (log(m_M/m_masses[1])+1.) * GetBeta_0_0();
  }
  else {
    double sum = 0.;
    std::vector<Complex> sumVec (6,0.);
    for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
      for (unsigned int j=0; j<=1; j++) {         // spin l
	for (unsigned int k=0; k<=2; k++) {       // spin Z
	  m_spins[0] = k;
	  m_spins[1] = j;
	  m_spins[2] = i;
	  Complex M_0_0 = InfraredSubtractedME_0_0();
	  Complex M_0_1 = InfraredSubtractedME_0_1();
	  sum = sum + 2.*(M_0_0*conj(M_0_1)).real();
	  if (msg_LevelIsDebugging()) {
	    const std::vector<Complex>& res_copy = m_res.GetResult();
	    for (size_t r(0); r < res_copy.size(); ++r) {
	      sumVec[r] += 2.*(M_0_0*conj(res_copy[r])).real();
	    }
	  }
	}
      }
    }
    if (msg_LevelIsDebugging()) {
      for (size_t r(0); r < sumVec.size(); ++r) {
	sumVec[r] = 1./3.*sumVec[r];
      }
      m2_loop = DivArrC(sumVec[0],sumVec[1],sumVec[2],sumVec[3],sumVec[4],sumVec[5]);
      Print_Info();
    }
    // spin avarage over initial state
    sum = (1./3.)*sum;
    return sum;
  }
}

double W_To_Lepton_Neutrino::GetBeta_0_2() {
  return 1./2.*sqr(m_alpha/M_PI*(log(m_M/m_masses[1])+1.)) * GetBeta_0_0();
  // return 0.;
}

double W_To_Lepton_Neutrino::GetBeta_1_1(unsigned int a) {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin W
	for (unsigned int l=0; l<=1; l++) {     // spin gamma
	  m_spins[0] = k;
	  m_spins[1] = j;
	  m_spins[2] = i;
	  m_spins[3] = l;
	  Complex M_1_05 = InfraredSubtractedME_1_05(a);
	  sum = sum + (M_1_05*conj(M_1_05)).real();
	}
      }
    }
  }
  // spin avarage over initial state
  sum = (1./3.)*sum;
  
  msg_Debugging()<<"a " << a 
		 << "\n|M_1(k)|^2="<<1./(16.*M_PI*M_PI*M_PI)*sum
		 <<"\nS(k)|M_0|^2="<<Smod(a)*GetBeta_0_0()
		 <<"\nsum="<<1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0()
		 <<"\nS(k)="<<Smod(a)
		 <<"\n|M_0|^2="<<GetBeta_0_0()<<std::endl<<std::endl;
  
  msg_Debugging() << "Photon energy " << m_moms1[a][3][0] << "\nPerc. diff. over M_1_1=" << (1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0())/(1./(16.*M_PI*M_PI*M_PI)*sum)
		  << "\nPerc. diff. over Smod approx.=" << (1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0())/(Smod(a)*GetBeta_0_0()) << std::endl << std::endl;

  sum = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0();
  return sum;
}

double W_To_Lepton_Neutrino::GetBeta_1_2(unsigned int a) {
  return 0.;
}

double W_To_Lepton_Neutrino::GetBeta_2_2(unsigned int a, unsigned int b) {
  return 0.;
}

void W_To_Lepton_Neutrino::Print_Info() {
  msg_Out() << " --------------------------- \n"
	    << "MomOne" << " " << m_moms0[0][0] << " " << m_moms0[0][1] << " " << m_moms0[0][2] 
	       << " " << m_moms0[0][3] << " " << m_flavs[0].Mass() 
	       << " , " << m_flavs[0] << " " << sqrt(m_moms0[0]*m_moms0[0]) <<  "\n"
	    << "MomTwo" << " " << m_moms0[1][0] << " " << m_moms0[1][1] << " " << m_moms0[1][2] 
	       << " " << m_moms0[1][3] << " " << m_flavs[1].Mass() 
	       << " , " << m_flavs[1] << " " << sqrt(m_moms0[1]*m_moms0[1]) <<  "\n"
	    << "MomThree" << " " << m_moms0[2][0] << " " << m_moms0[2][1] << " " << m_moms0[2][2] 
	       << " " << m_moms0[2][3] << " " << m_flavs[2].Mass() 
	       << " , " << m_flavs[2] << " " << sqrt(m_moms0[2]*m_moms0[2]) <<  "\n"
	    << "\n"
            << "Tree         " << GetBeta_0_0() << "\n"
	    << "Loop UV:     " << m_res.UV() << "\n"
	    << "Loop IR:     " << m_res.IR() << "\n"
	    << "Loop IR2:    " << m_res.IR2() << "\n"
	    << "Loop ep^0:   " << m_res.Finite() << "\n\n"
	    << "Loop^2 UV:    " << m2_loop.UV() << "\n"
	    << "Loop^2 IR:    " << m2_loop.IR() << "\n"
	    << "Loop^2 IR2:   " << m2_loop.IR2() << "\n\n"
	    << "Loop^2 ep^0:  " << m2_loop.Finite() << "\n"
    	    << "Loop^2 ep^-1: " << m2_loop.IR()+m2_loop.UV() << "\n" 
    	    << "Loop^2 ep^-2: " << m2_loop.IR2() << "\n"
	    << "\n";
}

DECLARE_PHOTONS_ME_GETTER(W_To_Lepton_Neutrino,
                          "W_To_Lepton_Neutrino")

PHOTONS_ME_Base *ATOOLS::Getter<PHOTONS_ME_Base,Particle_Vector_Vector,
				W_To_Lepton_Neutrino>::
operator()(const Particle_Vector_Vector &pvv) const
{
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 1) && (pvv[0][0]->Flav().Kfcode() == kf_Wplus) &&
       (pvv[1].size() == 0) &&
       (pvv[2].size() == 1) && pvv[2][0]->Flav().IsLepton() &&
       (pvv[3].size() == 1) && pvv[3][0]->Flav().IsLepton() )
    return new W_To_Lepton_Neutrino(pvv);
  return NULL;
}
