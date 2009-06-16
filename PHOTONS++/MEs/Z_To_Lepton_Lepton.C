#include "PHOTONS++/MEs/Z_To_Lepton_Lepton.H"

using namespace PHOTONS;
using namespace ATOOLS;
using namespace HELICITIES;
using namespace std;

Z_To_Lepton_Lepton::Z_To_Lepton_Lepton(Particle_Vector_Vector pvv_zero, Particle_Vector_Vector pvv_one) : Dipole_FF(pvv_zero) {
  m_flavs[0] = (pvv_zero.at(1)).at(0)->Flav();
  // switch ordering if necessary
  m_switch = (pvv_zero.at(2)).at(0)->Flav().IsAnti();
  // m_switch == true if first multipole particle is anti
  if (m_switch == false) {
    m_flavs[1] = (pvv_zero.at(2)).at(0)->Flav();
    m_flavs[2] = (pvv_zero.at(2)).at(1)->Flav();
  }
  else {
    m_flavs[2] = (pvv_zero.at(2)).at(0)->Flav();
    m_flavs[1] = (pvv_zero.at(2)).at(1)->Flav();
  }
  for (unsigned int i=3; i<9; i++) {
    m_flavs[i] = Flavour(kf_photon);
  }

  m_cL = m_i*m_e/(2*m_sW*m_cW)*(-1.+2.*m_sW*m_sW);
  m_cR = m_i*m_e/(2*m_sW*m_cW)*(2.*m_sW*m_sW);

  m_pvv_zero = pvv_zero;
  m_pvv_one  = pvv_one;

  FillMomentumArrays();
}

Z_To_Lepton_Lepton::~Z_To_Lepton_Lepton() {
}

void Z_To_Lepton_Lepton::BoostOriginalPVVToMultipoleCMS() {
  // m_pvv_one already in multipole CMS
  // m_pvv_zero in arbitrary frame -> boost m_olddipole into its CMS
  // and rotate m_olddipole.at(0) into +z direction
  Vec4D sum(0.,0.,0.,0.);
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    sum = sum + m_olddipole.at(i)->Momentum();
  }
  Vec4D p1 = m_olddipole.at(0)->Momentum();
  p_boost = new Poincare(sum);
  p_boost->Boost(p1);
  p_rot   = new Poincare(p1,Vec4D(0.,0.,0.,1.));
  for (unsigned int i=0; i<m_olddipole.size(); i++) {
    Vec4D vec = m_olddipole.at(i)->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_olddipole.at(i)->SetMomentum(vec);
  }
  for (unsigned int i=0; i<m_oldspectator.size(); i++) {
    Vec4D vec = m_oldspectator.at(i)->Momentum();
    p_boost->Boost(vec);
    p_rot->Rotate(vec);
    m_oldspectator.at(i)->SetMomentum(vec);
  }
}

void Z_To_Lepton_Lepton::FillMomentumArrays() {
  // m_moms0 - no photon
  m_moms0[0] = m_pvv_zero.at(1).at(0)->Momentum();
  if (m_switch == false) {
    m_moms0[1] = m_pvv_zero.at(2).at(0)->Momentum();
    m_moms0[2] = m_pvv_zero.at(2).at(1)->Momentum();
  }
  else {
    m_moms0[2] = m_pvv_zero.at(2).at(0)->Momentum();
    m_moms0[1] = m_pvv_zero.at(2).at(1)->Momentum();
  }
  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated

  // not necessary if only one photon
  if (m_pvv_one.at(4).size() == 1) {
    m_moms1[0][0] = m_pvv_one.at(1).at(0)->Momentum();
    if (m_switch == false) {
      m_moms1[0][1] = m_pvv_one.at(2).at(0)->Momentum();
      m_moms1[0][2] = m_pvv_one.at(2).at(1)->Momentum();
    }
    else {
      m_moms1[0][2] = m_pvv_one.at(2).at(0)->Momentum();
      m_moms1[0][1] = m_pvv_one.at(2).at(1)->Momentum();
    }
    m_moms1[0][3] = m_pvv_one.at(4).at(0)->Momentum();
  }
  else {
    Dipole_FF::DefineDipole();
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<m_pvv_one.at(4).size(); i++) {
      m_softphotons.push_back(m_pvv_one.at(4).at(i));
      m_K = CalculateMomentumSum(m_softphotons);
      CorrectMomenta();
      if (m_switch == false) {
        m_moms1[i][1] = m_newdipole.at(0)->Momentum();
        m_moms1[i][2] = m_newdipole.at(1)->Momentum();
      }
      else {
        m_moms1[i][2] = m_newdipole.at(0)->Momentum();
        m_moms1[i][1] = m_newdipole.at(1)->Momentum();
      }
      m_moms1[i][3] = m_softphotons.at(0)->Momentum();
      m_moms1[i][0] = m_moms1[i][1]+m_moms1[i][2]+m_moms1[i][3];
      m_softphotons.clear();
    }
  }
}

double Z_To_Lepton_Lepton::Smod(unsigned int kk) {
  m_moms = m_moms1[kk];
  double sum = 0;
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = -1;
  double Zj = +1;
  int    ti = +1;
  int    tj = +1;
  sum = m_alpha/(4*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k))*(pi/(pi*k)-pj/(pj*k));
  return sum;
}

Complex Z_To_Lepton_Lepton::InfraredSubtractedME_0_0() {
  m_moms = m_moms0;
  Vec4C epsZ = Polarization_Vector(m_moms[0]).at(m_spins[0]);
  XYZFunc XYZ(3,m_moms,m_flavs,1,false);
  return  XYZ.X(1,m_spins[1],epsZ,2,m_spins[2],m_cL,m_cR);
}

Complex Z_To_Lepton_Lepton::InfraredSubtractedME_0_1() {
  return 0;
}

Complex Z_To_Lepton_Lepton::InfraredSubtractedME_0_2() {
  return 0;
}

Complex Z_To_Lepton_Lepton::InfraredSubtractedME_1_05(unsigned int i) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsZ   = Polarization_Vector(m_moms[0]).at(m_spins[0]);
  Vec4C epsP   = conj(Polarization_Vector(m_moms[3]).at(m_spins[3]));
  Vec4D pa     = m_moms[1]+m_moms[3];       // fermion propagator momenta
  Vec4D pb     = m_moms[2]+m_moms[3];
  double m     = m_flavs[1].HadMass();       // fermion mass/propagator pole
  m_moms[4]    = m_moms[5] = pa;            // enter those into m_moms
  m_moms[6]    = m_moms[7] = pb;
  m_flavs[4]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[7] = m_flavs[2];
  XYZFunc XYZ(8,m_moms,m_flavs,1,false);
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Complex r4 = Complex(0.,0.);
  // emission off fermion (a)
  for (unsigned int s=0; s<=1; s++) { // spin of pseudo-particle in propagator representation
    r1 = r1 + XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)*XYZ.X(4,s,epsZ,2,m_spins[2],m_cL,m_cR);
    r2 = r2 + XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)*XYZ.X(5,s,epsZ,2,m_spins[2],m_cL,m_cR);
  }
  // emission off anti-fermion (b)
  for (unsigned int s=0; s<=1; s++) { // spin of pseudo-particle in propagator representation
    r3 = r3 + XYZ.X(1,m_spins[1],epsZ,6,s,m_cL,m_cR)*XYZ.X(6,s,epsP,2,m_spins[2],1.,1.);
    r4 = r4 + XYZ.X(1,m_spins[1],epsZ,7,s,m_cL,m_cR)*XYZ.X(7,s,epsP,2,m_spins[2],1.,1.);
  }
  // add prefactors
  r1 = m_e/(2*(pa*pa-m*m))*(1+m/sqrt(pa*pa))*r1;
  r2 = m_e/(2*(pa*pa-m*m))*(1-m/sqrt(pa*pa))*r2;
  r3 = -m_e/(2*(pb*pb-m*m))*(1-m/sqrt(pb*pb))*r3;
  r4 = -m_e/(2*(pb*pb-m*m))*(1+m/sqrt(pb*pb))*r4;
  // erase intermediate entries from m_flavs
  m_flavs[4] = m_flavs[5] = m_flavs[6] = m_flavs[7] = Flavour(kf_none);
  return (r1+r2+r3+r4);
}

Complex Z_To_Lepton_Lepton::InfraredSubtractedME_1_15(unsigned int i) {
  return 0;
}

Complex Z_To_Lepton_Lepton::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  return 0;
}

double Z_To_Lepton_Lepton::GetBeta_0_0() {
  double sum = 0;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
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

double Z_To_Lepton_Lepton::GetBeta_0_1() {
  return m_alpha/M_PI * (2*log(m_M/m_flavs[1].HadMass())+3./2.) * GetBeta_0_0();
}

double Z_To_Lepton_Lepton::GetBeta_0_2() {
  return 1./2.*pow(m_alpha/M_PI*2*log(m_M/m_flavs[1].HadMass()),2) * GetBeta_0_0();
}

double Z_To_Lepton_Lepton::GetBeta_1_1(unsigned int a) {
  double sum = 0;
  for (unsigned int i=0; i<=1; i++) {           // spin l.Bar
    for (unsigned int j=0; j<=1; j++) {         // spin l
      for (unsigned int k=0; k<=2; k++) {       // spin Z
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
  sum = 1./(2*pow(2*M_PI,3))*sum - Smod(a)*GetBeta_0_0();
  return sum;
}

double Z_To_Lepton_Lepton::GetBeta_1_2(unsigned int i) {
  return 0;
}

double Z_To_Lepton_Lepton::GetBeta_2_2(unsigned int i, unsigned int j) {
  return 0;
}
