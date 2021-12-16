#include "PHOTONS++/MEs/Higgs_To_Fermion_Fermion.H"
#include "PHOTONS++/MEs/Higgs_Decay_EW_One_Loop_Functions.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"
#include "ATOOLS/Org/Message.H"
#ifdef USING__YFS_NNLO
#include "PHOTONS++/MEs/RVTools/Higgs_Decay_Two_Loop_Functions.H"
#endif

#define B_0(A,B,C,M)         Master_Bubble(A,B,C,M)
#define C_0(A,B,C,D,E,F,M)   Master_Triangle(A,B,C,D,E,F,M)

using namespace PHOTONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Higgs_To_Fermion_Fermion::Higgs_To_Fermion_Fermion
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FF(pvv) {
  m_name = "Higgs_To_Fermion_Fermion";
  m_no_weight = 0;
  m_flavs[0] = pvv[1][0]->Flav();
  m_masses[0] = pvv[1][0]->FinalMass();
  // switch ordering if necessary
  m_switch = pvv[2][0]->Flav().IsAnti();
  // m_switch == true if first multipole particle is anti
  if (!m_switch) {
    m_flavs[1] = pvv[2][0]->Flav(); m_masses[1] = pvv[2][0]->FinalMass();
    m_flavs[2] = pvv[2][1]->Flav(); m_masses[2] = pvv[2][1]->FinalMass();
  }
  else {
    m_flavs[2] = pvv[2][0]->Flav(); m_masses[2] = pvv[2][0]->FinalMass();
    m_flavs[1] = pvv[2][1]->Flav(); m_masses[1] = pvv[2][1]->FinalMass();
  }
  for (unsigned int i=3; i<9; i++) {
    m_flavs[i]  = Flavour(kf_photon);
    m_masses[i] = 0.;
  }
  
  // Figure out isospin partner: Ordering in kf_code goes as downtype, uptype hence
  // If uptype, create Flavour with kf_code-1 and m_anti
  // If downtype, create Flavour with kf_code+1 and m_anti
  // Only makes sense if m_flavs[1] is fermion (which it should always be)
  Flavour IsoPartner;
  if (m_flavs[1].IsUptype()) {
    IsoPartner = ATOOLS::Flavour(m_flavs[1].Kfcode()-1,m_flavs[1].IsAnti());
  }
  else if (m_flavs[1].IsDowntype()) {
    IsoPartner = ATOOLS::Flavour(m_flavs[1].Kfcode()+1,m_flavs[1].IsAnti());
  }
  else {
    IsoPartner = ATOOLS::Flavour(0);
  }
  m2 = pow(m_flavs[1].Mass(),2.);
  Qf = m_flavs[1].Charge();
  If = m_flavs[1].IsoWeak();
  m2p = (IsoPartner.Kfcode()==6?mt2:pow(IsoPartner.Mass(),2.));
  Qfp = IsoPartner.Charge();
  Ifp = IsoPartner.IsoWeak();
  
  // Set couplings
  m_cL = -m_i*m_e*sqrt(m2)/(2.*m_sW*sqrt(MW2));
  m_cR = -m_i*m_e*sqrt(m2)/(2.*m_sW*sqrt(MW2));
}

Higgs_To_Fermion_Fermion::~Higgs_To_Fermion_Fermion() {
}

void Higgs_To_Fermion_Fermion::BoostOriginalPVVToMultipoleCMS() {
  // m_pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole[0] into +z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum += m_olddipole[i]->Momentum();
  }
  Vec4D p1 = m_olddipole[0]->Momentum();
  p_boost = new Poincare(sum);
  p_boost->Boost(p1);
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,1.));
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

void Higgs_To_Fermion_Fermion::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  // m_moms0 - no photon
  m_moms0[0] = m_pvv_zero[1][0]->Momentum();
  if (!m_switch) {
    m_moms0[1] = m_pvv_zero[2][0]->Momentum();
    m_moms0[2] = m_pvv_zero[2][1]->Momentum();
  }
  else {
    m_moms0[2] = m_pvv_zero[2][0]->Momentum();
    m_moms0[1] = m_pvv_zero[2][1]->Momentum();
  }

  // Store fully rescaled final state
  m_momsFull[0] = pvv_one[1][0]->Momentum();
  if (m_switch == false) {
    m_momsFull[1] = pvv_one[2][0]->Momentum();
    m_momsFull[2] = pvv_one[2][1]->Momentum();
  }
  else {
    m_momsFull[2] = pvv_one[2][0]->Momentum();
    m_momsFull[1] = pvv_one[2][1]->Momentum();
  }
  for (size_t i = 0; i < pvv_one[4].size(); i++) {
    m_momsFull[i+3] = pvv_one[4][i]->Momentum();
  }

  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated
  Dipole_FF::DefineDipole();

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[1][0]->Momentum();
    if (!m_switch) {
      m_moms1[0][1] = pvv_one[2][0]->Momentum();
      m_moms1[0][2] = pvv_one[2][1]->Momentum();
    }
    else {
      m_moms1[0][2] = pvv_one[2][0]->Momentum();
      m_moms1[0][1] = pvv_one[2][1]->Momentum();
    }
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
      if (!m_switch) {
        m_moms1[i][1] = m_newdipole[0]->Momentum();
        m_moms1[i][2] = m_newdipole[1]->Momentum();
      }
      else {
        m_moms1[i][2] = m_newdipole[0]->Momentum();
        m_moms1[i][1] = m_newdipole[1]->Momentum();
      }
      m_moms1[i][3] = m_softphotons[0]->Momentum();
      m_moms1[i][0] = m_moms1[i][1]+m_moms1[i][2]+m_moms1[i][3];
      m_softphotons.clear();
    }
  }

#ifdef USING__YFS_NNLO
  // m_moms2 - project multiphoton state onto two photon phase space
  // do reconstruction procedure again pretending only two photons were generated

  // not necessary if only two photons
  if (pvv_one[4].size() == 2) {
    m_moms2[0][1][0] = pvv_one[1][0]->Momentum();
    if (m_switch == false) {
      m_moms2[0][1][1] = pvv_one[2][0]->Momentum();
      m_moms2[0][1][2] = pvv_one[2][1]->Momentum();
    }
    else {
      m_moms2[0][1][2] = pvv_one[2][0]->Momentum();
      m_moms2[0][1][1] = pvv_one[2][1]->Momentum();
    }
    m_moms2[0][1][3] = pvv_one[4][0]->Momentum();
    m_moms2[0][1][4] = pvv_one[4][1]->Momentum();
  }
  else if (pvv_one[4].size() > 2) {
    BoostOriginalPVVToMultipoleCMS();
    // have pairs of momenta as hard photons produced - different momentum reconstruction for each
    // taje all possible pairs to account for exchange P1 <-> P2 - doesn't work
    for (unsigned int j=0; j<pvv_one[4].size(); j++) { 
      m_softphotons.push_back(pvv_one[4][j]);
      for (unsigned int i=0; i<j; i++) { 
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
	if (m_switch == false) {
	  m_moms2[i][j][1] = m_newdipole[0]->Momentum();
	  m_moms2[i][j][2] = m_newdipole[1]->Momentum();
	}
	else {
	  m_moms2[i][j][2] = m_newdipole[0]->Momentum();
	  m_moms2[i][j][1] = m_newdipole[1]->Momentum();
	}
	// Note choice of softphotons-index: Photon j is stored at position 0, Photon i at position 1
	// but m_moms2 ordered by Photon i, then Photon j
	m_moms2[i][j][3] = m_softphotons[1]->Momentum();
	m_moms2[i][j][4] = m_softphotons[0]->Momentum();
	m_moms2[i][j][0] = m_moms2[i][j][1]+m_moms2[i][j][2]+m_moms2[i][j][3]+m_moms2[i][j][4];
	m_softphotons.pop_back();
      }
      m_softphotons.clear();
    }
  }
#endif
}

double Higgs_To_Fermion_Fermion::SmodFull(unsigned int kk) {
  // Smod calculated using fully dressed momenta 
  // Used to cancel the Smod in the dipole weight

  // Index kk denotes which photon is used in calculation. Note that m_momsFull contains all photons
  // starting at position 3, hence index of photon kk is 3+kk.
  m_moms = m_momsFull;
  Vec4D k   = m_moms[3+kk];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

double Higgs_To_Fermion_Fermion::Smod(unsigned int kk) {
  // Smod calculated with momenta dressed with one additional hard photon
  // Used in the subtraction in single-real matrix elements

  // Index kk denotes which photon is taken to be hard - determines which set of momenta m_moms1[kk]
  // needs to be called. Note that the photon is always at position 3 in the sets.
  m_moms = m_moms1[kk];
  double sum = 0.;
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  sum = m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
  return sum;
}

double Higgs_To_Fermion_Fermion::Smod(unsigned int a, unsigned int b, unsigned int c) {
  // Smod calculated with momenta dressed with two additional hard photons
  // Used in subtraction in double-real matrix elements

  // Index a denotes position of first momentum, index b position of second photon. Used to get
  // correct set of dressed momenta m_moms2[a][b]
  // Index c denotes whether the photon in the calculation is the first or second of the two 
  // (c is either 0 or 1)
  m_moms = m_moms2[a][b];
  Vec4D k   = m_moms[3+c];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex Higgs_To_Fermion_Fermion::InfraredSubtractedME_0_0() {
  // Born level amplitude - simple Y-function
  m_moms = m_moms0;
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  return XYZ.Y(1,m_spins[1],2,m_spins[2],m_cR,m_cL);
}

Complex Higgs_To_Fermion_Fermion::InfraredSubtractedME_0_1(const int& ewmode) {
  m_moms = m_moms0;
  // Set parameters and XYZFunc
  double m2   = sqr(0.5*(m_masses[1]+m_masses[2]));
  double s    = sqr(m_masses[0]);
  double mu2  = s;
  double p1p2 = -m_moms[1]*m_moms[2];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
    
  // Set up EW corrections as corrections to cR,cL
  // If only QED corrections, all weak form factors and weak contributions in CTs are neglected
  Higgs_Decay_EW_One_Loop_Functions* OL = new Higgs_Decay_EW_One_Loop_Functions(s, m2, m2p, Qf, If, mu2, m_ew);
  Complex pref = Complex(0.,-m_e/(2.*m_sW*sqrt(MW2))); // checked
  // Contributions proportional to Y(1,2)
  DivArrC cRDivArr = m_alpha/M_PI*(
				   -pref*sqrt(m2)/(2.*m_sW2)*(OL->FSQED() + OL->FS())
				   +OL->CT_R() //checked
				   );
  DivArrC cLDivArr = m_alpha/M_PI*(
				   -pref*sqrt(m2)/(2.*m_sW2)*(OL->FSQED() + OL->FS())
				   +OL->CT_L() //checked
				   );
  if (msg_LevelIsDebugging()) {
    m_nlo_res = DivArrC(XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.UV(),cLDivArr.UV()),
			XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.IR(),cLDivArr.IR()),
			XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.IR2(),cLDivArr.IR2()),
			XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.Finite(),cLDivArr.Finite()),
			XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.Epsilon(),cLDivArr.Epsilon()),
			XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.Epsilon2(),cLDivArr.Epsilon2()))
      +m_alpha/M_PI*sqr(Qf)*(0.5*(s-2.*m2)*C_0(m2,m2,s,m2,0.,m2,mu2)
			     +m2*C_0(m2,m2,0.,0.,m2,m2,mu2)
			     +0.25*(B_0(s,m2,m2,mu2)-B_0(0.,m2,m2,mu2)))
      *XYZ.Y(1,m_spins[1],2,m_spins[2],m_cR,m_cL);
  }
  Complex term(0.,0.);
  // Assemble finite parts. Second term is infrared subtraction term B
  term += XYZ.Y(1,m_spins[1],2,m_spins[2],cRDivArr.Finite(),cLDivArr.Finite())
    +m_alpha/M_PI*sqr(Qf)*(0.5*(s-2.*m2)*C_0(m2,m2,s,m2,0.,m2,mu2)
			   +m2*C_0(m2,m2,0.,0.,m2,m2,mu2)
			   +0.25*(B_0(s,m2,m2,mu2)-B_0(0.,m2,m2,mu2))).Finite()
    *XYZ.Y(1,m_spins[1],2,m_spins[2],m_cR,m_cL);
  delete OL;
  return term;
}

Complex Higgs_To_Fermion_Fermion::InfraredSubtractedME_0_2() {
  // double virtual calculated in GetBeta_0_2()
  return 0.;
}

Complex Higgs_To_Fermion_Fermion::InfraredSubtractedME_1_05(unsigned int i, const double& xi) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4D real_moms[8];
  // Rescale momenta by xi - necessary for scaling test in RV
  for (size_t j = 0; j < 4; j++) real_moms[j] = xi*m_moms[j];
  Vec4C epsP   = conj(Polarization_Vector(real_moms[3])[m_spins[3]]);
  // Vec4C epsP   = real_moms[3]; // To check Ward identity
  Vec4D pa     = real_moms[1]+real_moms[3];       // fermion propagator momenta
  Vec4D pb     = real_moms[2]+real_moms[3];
  double m     = 0.5*(m_masses[1]+m_masses[2]);         // fermion mass/propagator pole
  real_moms[4]    = real_moms[5] = pa;            // enter those into real_moms
  real_moms[6]    = real_moms[7] = pb;
  m_flavs[4]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[7] = m_flavs[2];
  XYZFunc XYZ(8,real_moms,m_flavs,false);
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Complex r4 = Complex(0.,0.);
  // emission off fermions
  for (unsigned int s=0; s<=1; s++) { // spin of pseudo-particle in propagator representation
    r1 += XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)*XYZ.Y(4,s,2,m_spins[2],m_cR,m_cL);
    r2 += XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)*XYZ.Y(5,s,2,m_spins[2],m_cR,m_cL);
    r3 += XYZ.Y(1,m_spins[1],6,s,m_cR,m_cL)*XYZ.X(6,s,epsP,2,m_spins[2],1.,1.);
    r4 += XYZ.Y(1,m_spins[1],7,s,m_cR,m_cL)*XYZ.X(7,s,epsP,2,m_spins[2],1.,1.);
  }
  r1 *= m_e/(2.*(pa*pa-m*m))*(1.+m/sqrt(pa*pa));
  r2 *= m_e/(2.*(pa*pa-m*m))*(1.-m/sqrt(pa*pa));
  r3 *= -m_e/(2.*(pb*pb-m*m))*(1.-m/sqrt(pb*pb));
  r4 *= -m_e/(2.*(pb*pb-m*m))*(1.+m/sqrt(pb*pb));
  // erase intermediate entries from m_flavs
  m_flavs[4] = m_flavs[5] = m_flavs[6] = m_flavs[7] = Flavour(kf_none);

  return (r1+r2+r3+r4);
}

Complex Higgs_To_Fermion_Fermion::InfraredSubtractedME_1_15(unsigned int a, const double& xi) {
#ifdef USING__YFS_NNLO
  m_moms = m_moms1[a];
  Complex res = Complex(0.,0.);
  // Choose respective two loop functions instance 
  // Also define scaled momenta to be used in this scope - won't change m_moms
  Vec4D real_moms[6];
  Higgs_Decay_Two_Loop_Functions * TL;
  if (xi == 1.) TL = p_TL;
  else if (xi != 1.) TL = p_TL_scaled;
  for (size_t j = 0; j < 4; j++) real_moms[j] = xi*m_moms[j];
  Vec4C epsP   = conj(Polarization_Vector(real_moms[3])[m_spins[3]]);
  //Vec4C epsP   = real_moms[3];
  Vec4D pa     = real_moms[1]+real_moms[3];       // fermion propagator momenta
  Vec4D pb     = real_moms[2]+real_moms[3];
  double m     = 0.5*(m_masses[1]+m_masses[2]); // fermion mass/propagator pole
  double m2    = sqr(m);
  double s     = (real_moms[1]+real_moms[2]+real_moms[3]).Abs2();
  double s12   = (real_moms[1]+real_moms[2]).Abs2();
  double mu2(s);
  real_moms[4]    = real_moms[5] = real_moms[3];     // used for Standard MEs to replace \slashed{k}
  m_flavs[4] = Flavour(kf_neutrino);
  m_flavs[5] = Flavour(kf_neutrino,1); // require massless particle for calculation of standard MEs
  XYZFunc XYZ(6,real_moms,m_flavs,false);

  // Initialize results arrays
  // BubbleRes: bubble insertions on internal propagators
  // VertexRes: vertex corrections
  // BoxRes: box diagrams
  // FermionCT: CT to internal propagator
  // VertexCT: vertex CTs
  // B: contribution due to B*M_1_05
  DivArrC BubbleRes = DivArrC(0.,0.,0.,0.,0.,0.);
  DivArrC VertexRes = DivArrC(0.,0.,0.,0.,0.,0.);
  DivArrC BoxRes = DivArrC(0.,0.,0.,0.,0.,0.);
  DivArrC VertexCT(0.,0.,0.,0.,0.,0.);
  DivArrC FermionCT(0.,0.,0.,0.,0.,0.);
  DivArrC B(0.,0.,0.,0.,0.,0.);

  // define standard matrix elements - each with left and right handed projection
  // require insertion of internal particle(s) for MEs with two or more slashed vectors
  std::vector<Complex> ME4 = {Complex(0.,0.),Complex(0.,0.)};
  for (int sa = 0; sa <= 1; ++sa) {
    ME4[0] += 0.5*(XYZ.X(1,m_spins[1],epsP,4,sa,1.,1.)*XYZ.Y(4,sa,2,m_spins[2],0.,m_cL)
      		   +XYZ.X(1,m_spins[1],epsP,5,sa,1.,1.)*XYZ.Y(5,sa,2,m_spins[2],0.,m_cL));
    ME4[1] += 0.5*(XYZ.X(1,m_spins[1],epsP,4,sa,1.,1.)*XYZ.Y(4,sa,2,m_spins[2],m_cR,0.)
      		   +XYZ.X(1,m_spins[1],epsP,5,sa,1.,1.)*XYZ.Y(5,sa,2,m_spins[2],m_cR,0.));
  }


  // Standard MEs for either emission
  std::vector< std::vector<Complex> > StandardMEs = 
    {{XYZ.Y(1,m_spins[1],2,m_spins[2],0.,m_cL),XYZ.Y(1,m_spins[1],2,m_spins[2],m_cR,0.)}, // ubar1 P_i v2
     {XYZ.X(1,m_spins[1],real_moms[3],2,m_spins[2],0.,m_cL),XYZ.X(1,m_spins[1],real_moms[3],2,m_spins[2],m_cR,0.)}, // ubar1 \slashed{k} P_i v2
     {XYZ.X(1,m_spins[1],epsP,2,m_spins[2],0.,m_cL),XYZ.X(1,m_spins[1],epsP,2,m_spins[2],m_cR,0.)}, // ubar1 \slashed{epsP} P_i v2
     ME4}; // ubar1 \slashed{epsP} \slashed{k} P_i v2 

  // Calculate RV coefficients using respective polarization vector
  p_TL->Calculate_RV_Coeffs(epsP);

  for (int i = 0; i < 4; ++i) {
    for (int j = 0; j <= 1; ++j) {
      if (StandardMEs[i][j] != Complex(0.,0.)) {
	// Assemble results by multiplying respective coefficient with the correct standard ME
	// and prefactor
	BubbleRes += pow(m_e,3.)/(16.*sqr(M_PI))*StandardMEs[i][j]
	  *(p_TL->Get_RV_Bubble_Insertion_1(i,j)+p_TL->Get_RV_Bubble_Insertion_2(i,j));
	VertexRes += pow(m_e,3.)/(16.*sqr(M_PI))*StandardMEs[i][j]
	  *(p_TL->Get_RV_H_Vertex_1(i,j)+p_TL->Get_RV_H_Vertex_2(i,j)
	    +p_TL->Get_RV_P_Vertex_1(i,j)+p_TL->Get_RV_P_Vertex_2(i,j));
	BoxRes += pow(m_e,3.)/(16.*sqr(M_PI))*StandardMEs[i][j]
	  *(p_TL->Get_RV_Box_1(i,j)+p_TL->Get_RV_Box_2(i,j));
	FermionCT += m_e*m_alpha/M_PI*StandardMEs[i][j]*(p_TL->Get_RV_Fermion_CT(i,j));
	VertexCT += m_e*m_alpha/M_PI*StandardMEs[i][j]*(p_TL->Get_RV_Vertex_CT(i,j));
	B += m_e*m_alpha/M_PI*StandardMEs[i][j]*(p_TL->Get_RV_B(i,j));
      }
    }
  }
  // Assemble finite parts of result
  res += (B + FermionCT + VertexCT + BubbleRes + VertexRes + BoxRes).Finite();

  if (msg_LevelIsDebugging()) {
    if (xi == 1.) m_rv_res = B+BubbleRes+VertexRes+BoxRes+FermionCT+VertexCT;
    else m_rv_res_scaled = B+BubbleRes+VertexRes+BoxRes+FermionCT+VertexCT;
  }

  m_flavs[4] = m_flavs[5] = m_flavs[6] = m_flavs[7] = m_flavs[8] = m_flavs[9] = Flavour(kf_none);
  return res;
#else 
  return 0.;
#endif
}

Complex Higgs_To_Fermion_Fermion::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  // Calculate each diagram separately (in contrast to IR-subtracted 1_05!!!)
  // => pa, pb have different meanings (different prop. moms. in same diagram)!
  // Need to also calculate each diagram with photon momenta swapped since we only call i < j
#ifdef USING__YFS_NNLO
  Complex res = Complex(0.,0.);
  m_moms = m_moms2[i][j]; // set to set of momenta to be used: assuming i and j are the hard photons produced
  Vec4C epsP1 = conj(Polarization_Vector(m_moms[3])[m_spins[3]]); // polarization vecs.
  Vec4C epsP2 = conj(Polarization_Vector(m_moms[4])[m_spins[4]]);
  // Vec4C epsP1 = m_moms[3]; // To check Ward identity
  // Vec4C epsP2 = m_moms[4];

  // Diagram 1: both emissions of fermion (= m_moms[1]) leg
  Vec4D pa     = m_moms[1]+m_moms[3];           // fermion propagator momenta
  Vec4D pb     = m_moms[1]+m_moms[3]+m_moms[4];     
  Vec4D pc     = m_moms[1]+m_moms[4];
  double m     = 0.5*(m_masses[1]+m_masses[2]); // fermion mass/propagator pole
  m_moms[5]    = m_moms[7] = pa;            // enter those into m_moms
  m_moms[6]    = m_moms[8] = pb;
  m_flavs[5]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[7]   = m_flavs[8] = m_flavs[2];
  XYZFunc XYZ11(9,m_moms,m_flavs,false);
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Complex r4 = Complex(0.,0.);
  for (unsigned int sa=0; sa<=1; sa++) {
    for (unsigned int sb=0; sb<=1; sb++) {
      r1 += XYZ11.X(1,m_spins[1],epsP1,5,sa,1.,1.)
	*XYZ11.X(5,sa,epsP2,6,sb,1.,1.)
	*XYZ11.Y(6,sb,2,m_spins[2],m_cR,m_cL);
      r2 += XYZ11.X(1,m_spins[1],epsP1,5,sa,1.,1.)
	*XYZ11.X(5,sa,epsP2,8,sb,1.,1.)
	*XYZ11.Y(8,sb,2,m_spins[2],m_cR,m_cL);
      r3 += XYZ11.X(1,m_spins[1],epsP1,7,sa,1.,1.)
	*XYZ11.X(7,sa,epsP2,6,sb,1.,1.)
	*XYZ11.Y(6,sb,2,m_spins[2],m_cR,m_cL);
      r4 += XYZ11.X(1,m_spins[1],epsP1,7,sa,1.,1.)
	*XYZ11.X(7,sa,epsP2,8,sb,1.,1.)
	*XYZ11.Y(8,sb,2,m_spins[2],m_cR,m_cL);
    }
  }
  r1 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pa*pa))*(1.+m/sqrt(pb*pb));
  r2 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pa*pa))*(1.-m/sqrt(pb*pb));
  r3 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pa*pa))*(1.+m/sqrt(pb*pb));
  r4 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pa*pa))*(1.-m/sqrt(pb*pb));
  res += r1+r2+r3+r4;
  
  // reset intermediate results. momenta/flavours will be overwritten
  r1 = Complex(0.,0.);
  r2 = Complex(0.,0.);
  r3 = Complex(0.,0.);
  r4 = Complex(0.,0.);

  // Diagram 1: Swap P1 <-> P2 - amounts to swapping pa <-> pc, epsP1 <-> epsP2
  m_moms[5]    = m_moms[7] = pc;            // enter those into m_moms
  m_moms[6]    = m_moms[8] = pb;
  m_flavs[5]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[7]   = m_flavs[8] = m_flavs[2];
  XYZFunc XYZ12(9,m_moms,m_flavs,false);
  for (unsigned int sc=0; sc<=1; sc++) {
    for (unsigned int sb=0; sb<=1; sb++) {
      r1 += XYZ12.X(1,m_spins[1],epsP2,5,sc,1.,1.)
	*XYZ12.X(5,sc,epsP1,6,sb,1.,1.)
	*XYZ12.Y(6,sb,2,m_spins[2],m_cR,m_cL);
      r2 += XYZ12.X(1,m_spins[1],epsP2,5,sc,1.,1.)
	*XYZ12.X(5,sc,epsP1,8,sb,1.,1.)
	*XYZ12.Y(8,sb,2,m_spins[2],m_cR,m_cL);
      r3 += XYZ12.X(1,m_spins[1],epsP2,7,sc,1.,1.)
	*XYZ12.X(7,sc,epsP1,6,sb,1.,1.)
	*XYZ12.Y(6,sb,2,m_spins[2],m_cR,m_cL);
      r4 += XYZ12.X(1,m_spins[1],epsP2,7,sc,1.,1.)
	*XYZ12.X(7,sc,epsP1,8,sb,1.,1.)
	*XYZ12.Y(8,sb,2,m_spins[2],m_cR,m_cL);
    }
  }
  r1 *= sqr(m_e)/(4.*(pc*pc-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pc*pc))*(1.+m/sqrt(pb*pb));
  r2 *= sqr(m_e)/(4.*(pc*pc-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pc*pc))*(1.-m/sqrt(pb*pb));
  r3 *= sqr(m_e)/(4.*(pc*pc-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pc*pc))*(1.+m/sqrt(pb*pb));
  r4 *= sqr(m_e)/(4.*(pc*pc-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pc*pc))*(1.-m/sqrt(pb*pb));
  res += r1+r2+r3+r4;

  // reset intermediate results. momenta/flavours will be overwritten
  r1 = Complex(0.,0.);
  r2 = Complex(0.,0.);
  r3 = Complex(0.,0.);
  r4 = Complex(0.,0.);

    
  // Diagram 2: one emission of fermion leg, other off antifermion leg - comes with a minus sign
  pa	 = m_moms[1]+m_moms[3];	      // fermion propagator momenta
  pb	 = m_moms[2]+m_moms[4];
  pc	 = m_moms[1]+m_moms[4];	      // fermion propagator momenta
  Vec4D pd = m_moms[2]+m_moms[3];
  m	 = 0.5*(m_masses[1]+m_masses[2]); // fermion mass/propagator pole
  m_moms[5]    = m_moms[7] = pa;            // enter those into m_moms
  m_moms[6]    = m_moms[8] = pb;
  m_flavs[5]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[7]   = m_flavs[8] = m_flavs[2];
  XYZFunc XYZ21(9,m_moms,m_flavs,false);
  for (unsigned int sa=0; sa<=1; sa++) {
    for (unsigned int sb=0; sb<=1; sb++) {
      r1 += XYZ21.X(1,m_spins[1],epsP1,5,sa,1.,1.)
	*XYZ21.Y(5,sa,6,sb,m_cR,m_cL)
	*XYZ21.X(6,sb,epsP2,2,m_spins[2],1.,1.);
      r2 += XYZ21.X(1,m_spins[1],epsP1,5,sa,1.,1.)
	*XYZ21.Y(5,sa,8,sb,m_cR,m_cL)
	*XYZ21.X(8,sb,epsP2,2,m_spins[2],1.,1.);
      r3 += XYZ21.X(1,m_spins[1],epsP1,7,sa,1.,1.)
	*XYZ21.Y(7,sa,6,sb,m_cR,m_cL)
	*XYZ21.X(6,sb,epsP2,2,m_spins[2],1.,1.);
      r4 += XYZ21.X(1,m_spins[1],epsP1,7,sa,1.,1.)
	*XYZ21.Y(7,sa,8,sb,m_cR,m_cL)
	*XYZ21.X(8,sb,epsP2,2,m_spins[2],1.,1.);
    }
  }
  r1 *= -sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pa*pa))*(1.-m/sqrt(pb*pb));
  r2 *= -sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pa*pa))*(1.+m/sqrt(pb*pb));
  r3 *= -sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pa*pa))*(1.-m/sqrt(pb*pb));
  r4 *= -sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pa*pa))*(1.+m/sqrt(pb*pb));
  res += r1+r2+r3+r4;

  // reset intermediate results. momenta/flavours will be overwritten
  r1 = Complex(0.,0.);
  r2 = Complex(0.,0.);
  r3 = Complex(0.,0.);
  r4 = Complex(0.,0.);

  // Diagram 2: Swap P1 <-> P2 - amounts to epsP1 <-> epsP2, pa -> pc, pb -> pd (k1 <-> k2)
  m_moms[5]    = m_moms[7] = pc;            // enter those into m_moms
  m_moms[6]    = m_moms[8] = pd;
  m_flavs[5]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[7]   = m_flavs[8] = m_flavs[2];
  XYZFunc XYZ22(9,m_moms,m_flavs,false);
  for (unsigned int sc=0; sc<=1; sc++) {
    for (unsigned int sd=0; sd<=1; sd++) {
      r1 += XYZ22.X(1,m_spins[1],epsP2,5,sc,1.,1.)
	*XYZ22.Y(5,sc,6,sd,m_cR,m_cL)
	*XYZ22.X(6,sd,epsP1,2,m_spins[2],1.,1.);
      r2 += XYZ22.X(1,m_spins[1],epsP2,5,sc,1.,1.)
	*XYZ22.Y(5,sc,8,sd,m_cR,m_cL)
	*XYZ22.X(8,sd,epsP1,2,m_spins[2],1.,1.);
      r3 += XYZ22.X(1,m_spins[1],epsP2,7,sc,1.,1.)
	*XYZ22.Y(7,sc,6,sd,m_cR,m_cL)
	*XYZ22.X(6,sd,epsP1,2,m_spins[2],1.,1.);
      r4 += XYZ22.X(1,m_spins[1],epsP2,7,sc,1.,1.)
	*XYZ22.Y(7,sc,8,sd,m_cR,m_cL)
	*XYZ22.X(8,sd,epsP1,2,m_spins[2],1.,1.);
    }
  }
  r1 *= -sqr(m_e)/(4.*(pc*pc-m*m)*(pd*pd-m*m))*(1.+m/sqrt(pc*pc))*(1.-m/sqrt(pd*pd));
  r2 *= -sqr(m_e)/(4.*(pc*pc-m*m)*(pd*pd-m*m))*(1.+m/sqrt(pc*pc))*(1.+m/sqrt(pd*pd));
  r3 *= -sqr(m_e)/(4.*(pc*pc-m*m)*(pd*pd-m*m))*(1.-m/sqrt(pc*pc))*(1.-m/sqrt(pd*pd));
  r4 *= -sqr(m_e)/(4.*(pc*pc-m*m)*(pd*pd-m*m))*(1.-m/sqrt(pc*pc))*(1.+m/sqrt(pd*pd));
  res += r1+r2+r3+r4;

  // reset intermediate results. momenta/flavours will be overwritten
  r1 = Complex(0.,0.);
  r2 = Complex(0.,0.);
  r3 = Complex(0.,0.);
  r4 = Complex(0.,0.);
    
  // Diagram 3: both emissions off antifermion leg
  pa	 = m_moms[2]+m_moms[3]+m_moms[4];	// fermion propagator momenta
  pb	 = m_moms[2]+m_moms[4];
  pc	 = m_moms[2]+m_moms[3];
  m	 = 0.5*(m_masses[1]+m_masses[2]); // fermion mass/propagator pole
  m_moms[5]    = m_moms[7] = pa;            // enter those into m_moms
  m_moms[6]    = m_moms[8] = pb;
  m_flavs[5]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[7]   = m_flavs[8] = m_flavs[2];
  XYZFunc XYZ31(9,m_moms,m_flavs,false);
  for (unsigned int sa=0; sa<=1; sa++) {
    for (unsigned int sb=0; sb<=1; sb++) {
      r1 += XYZ31.Y(1,m_spins[1],5,sa,m_cR,m_cL)
	*XYZ31.X(5,sa,epsP1,6,sb,1.,1.)
	*XYZ31.X(6,sb,epsP2,2,m_spins[2],1.,1.);
      r2 += XYZ31.Y(1,m_spins[1],5,sa,m_cR,m_cL)
	*XYZ31.X(5,sa,epsP1,8,sb,1.,1.)
	*XYZ31.X(8,sb,epsP2,2,m_spins[2],1.,1.);
      r3 += XYZ31.Y(1,m_spins[1],7,sa,m_cR,m_cL)
	*XYZ31.X(7,sa,epsP1,6,sb,1.,1.)
	*XYZ31.X(6,sb,epsP2,2,m_spins[2],1.,1.);
      r4 += XYZ31.Y(1,m_spins[1],7,sa,m_cR,m_cL)
	*XYZ31.X(7,sa,epsP1,8,sb,1.,1.)
	*XYZ31.X(8,sb,epsP2,2,m_spins[2],1.,1.);
    }
  }
  r1 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pa*pa))*(1.-m/sqrt(pb*pb));
  r2 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.-m/sqrt(pa*pa))*(1.+m/sqrt(pb*pb));
  r3 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pa*pa))*(1.-m/sqrt(pb*pb));
  r4 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pb*pb-m*m))*(1.+m/sqrt(pa*pa))*(1.+m/sqrt(pb*pb));
  res += r1+r2+r3+r4;

  // reset intermediate results. momenta/flavours will be overwritten
  r1 = Complex(0.,0.);
  r2 = Complex(0.,0.);
  r3 = Complex(0.,0.);
  r4 = Complex(0.,0.);

  // Diagram 3: Swap P1 <-> P2 - amounts to epsP1 <-> epsP2, pb <-> pc
  m_moms[5]    = m_moms[7] = pa;            // enter those into m_moms
  m_moms[6]    = m_moms[8] = pc;
  m_flavs[5]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[7]   = m_flavs[8] = m_flavs[2];
  XYZFunc XYZ32(9,m_moms,m_flavs,false);
  for (unsigned int sa=0; sa<=1; sa++) {
    for (unsigned int sc=0; sc<=1; sc++) {
      r1 += XYZ32.Y(1,m_spins[1],5,sa,m_cR,m_cL)
	*XYZ32.X(5,sa,epsP2,6,sc,1.,1.)
	*XYZ32.X(6,sc,epsP1,2,m_spins[2],1.,1.);
      r2 += XYZ32.Y(1,m_spins[1],5,sa,m_cR,m_cL)
	*XYZ32.X(5,sa,epsP2,8,sc,1.,1.)
	*XYZ32.X(8,sc,epsP1,2,m_spins[2],1.,1.);
      r3 += XYZ32.Y(1,m_spins[1],7,sa,m_cR,m_cL)
	*XYZ32.X(7,sa,epsP2,6,sc,1.,1.)
	*XYZ32.X(6,sc,epsP1,2,m_spins[2],1.,1.);
      r4 += XYZ32.Y(1,m_spins[1],7,sa,m_cR,m_cL)
	*XYZ32.X(7,sa,epsP2,8,sc,1.,1.)
	*XYZ32.X(8,sc,epsP1,2,m_spins[2],1.,1.);
    }
  }
  r1 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pc*pc-m*m))*(1.-m/sqrt(pa*pa))*(1.-m/sqrt(pc*pc));
  r2 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pc*pc-m*m))*(1.-m/sqrt(pa*pa))*(1.+m/sqrt(pc*pc));
  r3 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pc*pc-m*m))*(1.+m/sqrt(pa*pa))*(1.-m/sqrt(pc*pc));
  r4 *= sqr(m_e)/(4.*(pa*pa-m*m)*(pc*pc-m*m))*(1.+m/sqrt(pa*pa))*(1.+m/sqrt(pc*pc));
  res += r1+r2+r3+r4;

  // erase intermediate entries from m_flavs
  m_flavs[5] = m_flavs[6] = m_flavs[7] = m_flavs[8] = Flavour(kf_none);
  return res;
#else 
  return 0.;
#endif
}

double Higgs_To_Fermion_Fermion::GetBeta_0_0() {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=0; k++) {       // spin scalar
        m_spins[0] = k;
        m_spins[1] = j;
        m_spins[2] = i;
        Complex M_0_0 = InfraredSubtractedME_0_0();
        sum = sum + (M_0_0*conj(M_0_0)).real();
      }
    }
  }
  // spin avarage over initial state gives factor 1
  return sum;
}

double Higgs_To_Fermion_Fermion::GetBeta_0_1(const int& ewmode) {
  if (m_limit == 1) {
    return m_alpha/M_PI*(2.*log(m_M/(0.5*(m_masses[1]+m_masses[2])))-1.)
           *GetBeta_0_0();
  }
  else {
    int mode = m_ew&ewmode;
    double sum = 0.;
    std::vector<Complex> sumVec (6,0.);
    for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
      for (unsigned int j=0; j<=1; j++) {         // spin l
	for (unsigned int k=0; k<=0; k++) {       // spin scalar
	  m_spins[0] = k;
	  m_spins[1] = j;
	  m_spins[2] = i;
	  Complex M_0_0 = InfraredSubtractedME_0_0();
	  Complex M_0_1 = InfraredSubtractedME_0_1(mode);
	  sum = sum + 2.*(M_0_0*conj(M_0_1)).real();
	  
	  if (msg_LevelIsDebugging()) {
	    const std::vector<Complex>& res_copy = m_nlo_res.GetResult();
	    for (size_t r(0); r < res_copy.size(); ++r) {
	      sumVec[r] += 2.*(M_0_0*conj(res_copy[r])).real();
	    }
	  }
	}
      }
    }
    // spin avarage over initial state gives factor 1 - factor 3 for quark colours
    if (msg_LevelIsDebugging()) {
      m2_loop = DivArrC(sumVec[0],sumVec[1],sumVec[2],sumVec[3],sumVec[4],sumVec[5]);
      Print_Info();
    }
    return sum;
  }
}

double Higgs_To_Fermion_Fermion::GetBeta_0_2() {
  // rough estimate of leading log following YFS2
#ifdef USING__YFS_NNLO
  if (m_nnlo_qed == 1 || m_nnlo_qed == 4) {
    m_moms = m_moms0;
    return 1./2.*sqr(m_alpha/M_PI*2.*log(m_M/(0.5*(m_masses[1]+m_masses[2]))))
      *GetBeta_0_0();  
  }
  else return 0.;
#else 
  return 0.;
#endif
}

double Higgs_To_Fermion_Fermion::GetBeta_1_1(unsigned int a) {
  double sum = 0.;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=0; k++) {       // spin scalar
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
  // spin avarage over initial state gives factor 1
  if (msg_LevelIsDebugging()) {
    m_r2_res = sum;
    m_r2s_res = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0();
    msg_Debugging()<<"a " << a 
		   << "\n|M_1(k)|^2="<<1./(16.*M_PI*M_PI*M_PI)*sum
		   <<"\nS(k)|M_0|^2="<<Smod(a)*GetBeta_0_0()
		   <<"\nsum="<<1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0()
		   <<"\nS(k)="<<Smod(a)
		   <<"\n|M_0|^2="<<GetBeta_0_0()<<std::endl<<std::endl;
    
    msg_Debugging() << "Photon energy " << m_moms1[a][3][0] << "\nPerc. diff. over M_1_1=" << (1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0())/(1./(16.*M_PI*M_PI*M_PI)*sum)
		    << "\nPerc. diff. over Smod approx.=" << (1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0())/(Smod(a)*GetBeta_0_0()) << std::endl << std::endl;
  }

  sum = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0();
  return sum;
}

double Higgs_To_Fermion_Fermion::GetBeta_1_2(unsigned int a) {
  /// Berends, Neerven, Burgers / do reduction by self
#ifdef USING__YFS_NNLO
  if (m_nnlo_qed == 1) {
    m_moms = m_moms1[a];
    double m     = 0.5*(m_masses[1]+m_masses[2]); // fermion mass/propagator pole
    double m2    = sqr(m);
    double s     = (m_moms[1]+m_moms[2]+m_moms[3]).Abs2();
    double s1k   = (m_moms[1]+m_moms[3]).Abs2();
    double s2k   = (m_moms[2]+m_moms[3]).Abs2();
    // Calculate using soft-collinear approximation if m2/s1k or m2/s2k large enough
    double dmod = 0.;
    double softsub = 0.;
    double coll = 0.;
    if (m2/s1k > m_dev || m2/s2k > m_dev) {
      for (int i = 1; i < 3; i++) {
	Vec4D pi = m_moms[i];
	int j = (i==1)?2:1;
	Vec4D pj = m_moms[j];
	double pik = m_moms[3]*pi;
	double pjk = m_moms[3]*pj;
	double pipj = pi*pj;
	double Dsoft = 1./(pik)*(2.*(pipj)/(pik+pjk)-(pi*pi)/(pik));
	double y = (pik)/(pik+pjk+pipj);
	double z = (pipj)/(pjk+pipj);
	double p = s;
	double P = p - pi*pi - pj*pj - m_moms[3]*m_moms[3];
	double R = sqrt(sqr(2.*(pj*pj)+P-P*y)-4.*p*(pj*pj))/sqrt(sqr(P)-sqr(pi*pi)-sqr(pj*pj) - 2.*P*(pj*pj) - 2.*P*(pi*pi) - 2.*(pj*pj)*(pi*pi));
	coll += 1./((pik)*R)*(2./(1.-z*(1.-y))-1.-z-(pi*pi)/(pik));
	softsub += Dsoft;
	dmod += 1./((pik)*R)*(2./(1.-z*(1.-y))-1.-z-(pi*pi)/(pik))-Dsoft;
      }
      return -m_alpha/(4.*sqr(M_PI))*dmod*GetBeta_0_1(0);
    }
    double mu2(s);
    double sum = 0.;
    double sumscaled = 0.;
    Complex realsum = Complex(0.,0.);
    DivArrC ampsum = DivArrC(0.,0.,0.,0.,0.,0.);
    std::vector<Complex> sumVec (6,0.);
    // Set up (non-)rescaled two loop functions
    p_TL = new Higgs_Decay_Two_Loop_Functions(m, s, m_moms[1], m_moms[2], m_moms[3], m_cL, m_cR, mu2);
    p_TL_scaled = new Higgs_Decay_Two_Loop_Functions(m_xi*m, sqr(m_xi)*s, m_xi*m_moms[1], m_xi*m_moms[2], m_xi*m_moms[3], m_cL, m_cR, sqr(m_xi)*mu2);
    for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
      for (unsigned int j=0; j<=1; j++) {         // spin l
    	for (unsigned int k=0; k<=0; k++) {       // spin Z
    	  for (unsigned int l=0; l<=1; l++) {     // spin gamma
    	    m_spins[0] = k;
    	    m_spins[1] = j;
    	    m_spins[2] = i;
    	    m_spins[3] = l;
    	    Complex M_1_15 = InfraredSubtractedME_1_15(a);
    	    Complex M_1_05 = InfraredSubtractedME_1_05(a);
    	    Complex M_1_15scaled = InfraredSubtractedME_1_15(a,m_xi);
	    Complex M_1_05scaled = InfraredSubtractedME_1_05(a,m_xi);
	    realsum = realsum + M_1_05;
	    ampsum = ampsum + m_rv_res;
    	    sum = sum + (M_1_15*conj(M_1_05)+M_1_05*conj(M_1_15)).real();
    	    sumscaled = sumscaled + (M_1_15scaled*conj(M_1_05scaled)+M_1_05scaled*conj(M_1_15scaled)).real();
	    if (msg_LevelIsDebugging()) {
	      ampsum = ampsum + m_rv_res;
	      const std::vector<Complex>& res_copy = m_rv_res.GetResult();
	      for (size_t r(0); r < res_copy.size(); ++r) {
		sumVec[r] += 2.*(res_copy[r]*conj(M_1_05)+M_1_05*conj(res_copy[r])).real();
	      }
	    }
    	  }
    	}
      }
    }
    delete p_TL;
    delete p_TL_scaled;
    sum = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_1(0);
    sumscaled = 1./(16.*M_PI*M_PI*M_PI)*sumscaled - Smod(a)*GetBeta_0_1(0);
    if (abs(sum/sumscaled-1.) > m_dev) return 0.;
    if (msg_LevelIsDebugging()) {
      m_rv_res = ampsum;
      m_r_res = realsum;
      m2_rvfull = 1./(16.*M_PI*M_PI*M_PI)*DivArrC(sumVec[0],sumVec[1],sumVec[2],sumVec[3],sumVec[4],sumVec[5]) - Smod(a)*m2_loop;
      m2_rv = 1./(16.*M_PI*M_PI*M_PI)*DivArrC(sumVec[0],sumVec[1],sumVec[2],sumVec[3],sumVec[4],sumVec[5]);
      m_smod = Smod(a);
      Print_Info();
    }
    return sum;
  }
  else
    return 0.;
#else 
  return 0.;
#endif
}

double Higgs_To_Fermion_Fermion::GetBeta_2_2(unsigned int a, unsigned int b) {
#ifdef USING__YFS_NNLO
  if (m_nnlo_qed == 1) {
    double sum = 0.;
    for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
      for (unsigned int j=0; j<=1; j++) {         // spin l
	for (unsigned int k=0; k<=0; k++) {       // spin Higgs
	  for (unsigned int l=0; l<=1; l++) {     // spin gamma1
	    for (unsigned int m=0; m<=1; m++) {     // spin gamma2
	      m_spins[0] = k;
	      m_spins[1] = j;
	      m_spins[2] = i;
	      m_spins[3] = l;
	      m_spins[4] = m;
	      Complex M_2_1 = InfraredSubtractedME_2_1(a,b);
	      sum = sum + (M_2_1*conj(M_2_1)).real();
	    }
	  }
	}
      }
    }
    msg_Debugging()<<"|M_2(k)|^2="<<sqr(1./(16.*M_PI*M_PI*M_PI))*sum
		   <<", S(k)|M_1|^2 + S(k1)*S(k2)*|M_0|^2="<<Smod(a,b,0)*GetBeta_1_1(b) + Smod(a,b,1)*GetBeta_1_1(a) + Smod(a,b,0)*Smod(a,b,1)*GetBeta_0_0()
		   <<"before , S(k)|M_1|^2 + S(k1)*S(k2)*|M_0|^2="<<Smod(a)*GetBeta_1_1(b) + Smod(a)*Smod(b)*GetBeta_0_0()<<std::endl;

    sum = sqr(1./(16.*M_PI*M_PI*M_PI))*sum - Smod(a,b,0)*GetBeta_1_1(b) - Smod(a,b,1)*GetBeta_1_1(a) - Smod(a,b,0)*Smod(a,b,1)*GetBeta_0_0();
    return sum;
  }
  else 
    return 0.;
#else 
  return 0.;
#endif
}

void Higgs_To_Fermion_Fermion::Print_Info() {
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
	    << "NLO Loop UV:     " << m_nlo_res.UV() << "\n"
	    << "NLO Loop IR:	 " << m_nlo_res.IR() << "\n"
	    << "NLO Loop IR2:	 " << m_nlo_res.IR2() << "\n"
	    << "NLO Loop ep^0:	 " << m_nlo_res.Finite() << "\n\n"
	    << "NLO Loop^2 UV:	 " << m2_loop.UV() << "\n"
	    << "NLO Loop^2 IR:	 " << m2_loop.IR() << "\n"
	    << "NLO Loop^2 IR2:	 " << m2_loop.IR2() << "\n\n"
	    << "NLO Loop^2 ep^0:	 " << m2_loop.Finite() << "\n"
	    << "NLO Loop^2 ep^-1: " << m2_loop.IR()+/*3.**/m2_loop.UV() << "\n" 
	    << "NLO Loop^2 ep^-2: " << m2_loop.IR2() << "\n"
	    << "\n";

  if (m_nnlo_qed == 1) {
    msg_Out() << "RV UV:	   " << m_rv_res.UV() << "\n"
	      << "RV IR:	   " << m_rv_res.IR() << "\n"
	      << "RV IR2:    " << m_rv_res.IR2() << "\n"
	      << "RV ep^0:   " << m_rv_res.Finite() << "\n\n"
	      << "RV^2 UV:	   " << m2_rv.UV() << "\n"
	      << "RV^2 IR:	   " << m2_rv.IR() << "\n"
	      << "RV^2 IR2:   " << m2_rv.IR2() << "\n" 
	      << "Smod:      " << m_smod << "\n\n"
	      << "RV^2 ep^0:  " << m2_rv.Finite() << "\n"
	      << "RV^2 ep^-1: " << m2_rv.IR()+m2_rv.UV() << "\n" 
	      << "RV^2 ep^-2: " << m2_rv.IR2() << "\n"
	      << "\n"
	      << "RV^2 - Smod UV:	   " << m2_rvfull.UV() << "\n"
	      << "RV^2 - Smod IR:	   " << m2_rvfull.IR() << "\n"
	      << "RV^2 - Smod IR2:	    " << m2_rvfull.IR2() << "\n" 
	      << "Smod:      " << m_smod << "\n\n"
	      << "RV^2 - Smod ep^0:  " << m2_rvfull.Finite() << "\n"
	      << "RV^2 - Smod ep^-1: " << m2_rvfull.IR()+m2_rvfull.UV() << "\n" 
	      << "RV^2 - Smod ep^-2: " << m2_rvfull.IR2() << "\n"
	      << "\n";
  }
}

DECLARE_PHOTONS_ME_GETTER(Higgs_To_Fermion_Fermion,
                          "Higgs_To_Fermion_Fermion")

PHOTONS_ME_Base *ATOOLS::Getter<PHOTONS_ME_Base,Particle_Vector_Vector,
				Higgs_To_Fermion_Fermion>::
operator()(const Particle_Vector_Vector &pvv) const
{
  // same mass restriction can be lifted if M_0_1 is computed for general case
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 0) &&
       (pvv[1].size() == 1) && pvv[1][0]->Flav().Kfcode() == kf_h0 &&
       (pvv[2].size() == 2) && pvv[2][0]->Flav().IsFermion() &&
       (pvv[2][0]->Flav() == pvv[2][1]->Flav().Bar()) &&
       (pvv[3].size() == 0) )
    return new Higgs_To_Fermion_Fermion(pvv);
  return NULL;
}
