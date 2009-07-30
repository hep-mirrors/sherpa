#include "PHOTONS++/MEs/W_To_Lepton_Neutrino.H"
#include "ATOOLS/Math/Poincare.H"
#include "HELICITIES/Main/XYZFuncs.H"
#include "HELICITIES/Main/Polarization_Tools.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;

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
  m_name = "W_To_Lepton_Neutrino";
  m_flavs[0]  = pvv[0][0]->Flav();
  m_masses[0] = pvv[1][0]->FinalMass();
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
  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[2][0]->Momentum();
    m_moms1[0][1] = pvv_one[2][1]->Momentum();
    m_moms1[0][2] = pvv_one[3][0]->Momentum();
    m_moms1[0][3] = pvv_one[4][0]->Momentum();
  }
  else {
    Dipole_FI::DefineDipole();
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<pvv_one[4].size(); i++) {
      m_softphotons.push_back(pvv_one[4][i]);
      m_K = CalculateMomentumSum(m_softphotons);
      CorrectMomenta();
      m_moms1[i][0] = m_newdipole[0]->Momentum();
      m_moms1[i][1] = m_newdipole[1]->Momentum();
      m_moms1[i][2] = m_newspectator[0]->Momentum();
      m_moms1[i][3] = m_softphotons[0]->Momentum();
      m_softphotons.clear();
    }
  }
}

double W_To_Lepton_Neutrino::Smod(unsigned int kk) {
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

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_0_0() {
  m_moms = m_moms0;
  Vec4C epsW = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,1,false);
  return m_i*m_e/(m_sqrt2*m_sW)*XYZ.X(1,m_spins[1],epsW,2,m_spins[2],m_cL,m_cR);
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_0_1() {
  return 0.;
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_0_2() {
  return 0.;
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_1_05(unsigned int i) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsW   = Polarization_Vector(m_moms[0])[m_spins[0]];
  Vec4C epsP   = conj(Polarization_Vector(m_moms[3])[m_spins[3]]);
  Vec4D q      = m_moms[1]+m_moms[3];       // fermion propagator momenta
  Vec4D Q      = m_moms[0]+m_moms[3];       // boson propagator momenta
  double m     = m_masses[1];               // fermion mass/propagator pole
  double M     = m_masses[0];               // boson mass/propagator pole
  m_moms[4]    = m_moms[5] = q;             // enter those into m_moms
  m_flavs[4]   = m_flavs[1];                // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[1].Bar();
  XYZFunc XYZ(6,m_moms,m_flavs,1,false);
  m_flavs[4] = m_flavs[5] = Flavour(kf_none);
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Complex r4 = Complex(0.,0.);
  Complex r5 = Complex(0.,0.);
  Complex r6 = Complex(0.,0.);
  Complex r7 = Complex(0.,0.);
  Complex r8 = Complex(0.,0.);
  for (unsigned int s=0; s<=1; s++) {
    r1 += XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)
          *XYZ.X(4,s,epsW,2,m_spins[2],m_cL,m_cR);
    r2 += XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)
          *XYZ.X(5,s,epsW,2,m_spins[2],m_cL,m_cR);
  }
  Vec4D p = m_moms[0];
  Vec4D k = m_moms[3];
  r3 = XYZ.X(1,m_spins[1],epsW,2,m_spins[2],m_cL,m_cR);
  r4 = XYZ.X(1,m_spins[1],p+k,2,m_spins[2],m_cL,m_cR);
  r5 = XYZ.X(1,m_spins[1],epsP,2,m_spins[2],m_cL,m_cR);
  r6 = r7 = r8 = XYZ.X(1,m_spins[1],p-k,2,m_spins[2],m_cL,m_cR);
  r1 *= (1.+m/sqrt(q*q))/(q*q-m*m);
  r2 *= (1.-m/sqrt(q*q))/(q*q-m*m);
  r3 *= (epsP*(-2.*p+k))/(Q*Q-M*M);
  r4 *= (epsW*epsP)/(Q*Q-M*M);
  r5 *= (epsW*(p-2.*k))/(Q*Q-M*M);
  r6 *= -(epsW*(p-k))*(epsP*(-2.*p+k))/(Q*Q-M*M);
  r7 *= -(epsW*epsP)*((p-k)*(p+k))/(Q*Q-M*M);
  r8 *= -(epsW*(p-2.*k))*(epsP*(p-k))/(Q*Q-M*M);

  return (m_i*m_e*m_e)/(2.*m_sqrt2*m_sW)*(r1+r2+r3+r4+r5+r6+r7+r8);
}

Complex W_To_Lepton_Neutrino::InfraredSubtractedME_1_15(unsigned int i) {
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

double W_To_Lepton_Neutrino::GetBeta_0_1() {
  // limit mW > ml
  return m_alpha/M_PI * (log(m_M/m_masses[1])+1.) * GetBeta_0_0();
}

double W_To_Lepton_Neutrino::GetBeta_0_2() {
  return 0.;
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
  sum = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0();
//   msg_Out()<<sum<<endl;
  return sum;
}

double W_To_Lepton_Neutrino::GetBeta_1_2(unsigned int i) {
  return 0.;
}

double W_To_Lepton_Neutrino::GetBeta_2_2(unsigned int i, unsigned int j) {
  return 0.;
}

DECLARE_PHOTONS_ME_GETTER(W_To_Lepton_Neutrino_Getter,
                          "W_To_Lepton_Neutrino")

PHOTONS_ME_Base * W_To_Lepton_Neutrino_Getter::operator()
(const Particle_Vector_Vector &pvv) const
{
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 1) && (pvv[0][0]->Flav().Kfcode() == kf_Wplus) &&
       (pvv[1].size() == 0) &&
       (pvv[2].size() == 1) && pvv[2][0]->Flav().IsLepton() &&
       (pvv[3].size() == 1) && pvv[3][0]->Flav().IsLepton() )
    return new W_To_Lepton_Neutrino(pvv);
  return NULL;
}
