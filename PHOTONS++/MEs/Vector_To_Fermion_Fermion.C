#include "PHOTONS++/MEs/Vector_To_Fermion_Fermion.H"
#include "ATOOLS/Math/Poincare.H"
#include "ATOOLS/Phys/Flavour.H"
#include "ATOOLS/Phys/Particle.H"
#include "METOOLS/Main/XYZFuncs.H"
#include "METOOLS/Main/Polarization_Tools.H"
#include "METOOLS/Loops/Master_Integrals.H"
#include "METOOLS/Loops/PV_Integrals.H"

#define A_0(A,M)            Master_Tadpole(A,M)
#define B_0(A,B,C,M)        Master_Bubble(A,B,C,M)
#define B_1(A,B,C,M)        PV_Bubble_1(A,B,C,M)
#define C_11(A,B,C,D,E,F,M) PV_Triangle_11(A,B,C,D,E,F,M)
#define C_12(A,B,C,D,E,F,M) PV_Triangle_12(A,B,C,D,E,F,M)
#define C_21(A,B,C,D,E,F,M) PV_Triangle_21(A,B,C,D,E,F,M)
#define C_22(A,B,C,D,E,F,M) PV_Triangle_22(A,B,C,D,E,F,M)
#define C_23(A,B,C,D,E,F,M) PV_Triangle_23(A,B,C,D,E,F,M)
#define C_24(A,B,C,D,E,F,M) PV_Triangle_24(A,B,C,D,E,F,M)

using namespace PHOTONS;
using namespace ATOOLS;
using namespace METOOLS;
using namespace std;

Vector_To_Fermion_Fermion::Vector_To_Fermion_Fermion
(const Particle_Vector_Vector& pvv) : PHOTONS_ME_Base(pvv), Dipole_FF(pvv) {
  m_name = "Vector_To_Fermion_Fermion";
  m_flavs[0]  = pvv[1][0]->Flav();
  m_masses[0] = pvv[1][0]->FinalMass();
  // switch ordering if necessary
  m_switch = pvv[2][0]->Flav().IsAnti();
  // m_switch == true if first multipole particle is anti
  // such that m_flavs[1] is particle, m_flavs[2] is anti
  if (m_switch == false) {
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

  // Hadrons' form factors for basic process
  // set to one, will have to be got from Hadrons
  // for now:
  double F_L = 0.;
  double F_R = 0.;
  if (m_flavs[0] == Flavour(kf_Z)) {
    // full EW couplings
    F_L = -1./(2.*m_sW*m_cW)*(2.*m_flavs[1].IsoWeak()
                                  -2.*m_flavs[1].Charge()*m_sW*m_sW);
    F_R = -1./(2.*m_sW*m_cW)*(-2.*m_flavs[1].Charge()*m_sW*m_sW);
  }
  else {
    // assume electromagnetic/strong decay
    F_L = m_flavs[1].Charge();
    F_R = m_flavs[1].Charge();
  }
  m_cL = m_i*m_e*F_L;
  m_cR = m_i*m_e*F_R;
}

Vector_To_Fermion_Fermion::~Vector_To_Fermion_Fermion() {
}

void Vector_To_Fermion_Fermion::BoostOriginalPVVToMultipoleCMS() {
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

void Vector_To_Fermion_Fermion::FillMomentumArrays
(const Particle_Vector_Vector& pvv_one) {
  // m_moms0 - no photon
  m_moms0[0] = m_pvv_zero[1][0]->Momentum();
  if (m_switch == false) {
    m_moms0[1] = m_pvv_zero[2][0]->Momentum();
    m_moms0[2] = m_pvv_zero[2][1]->Momentum();
  }
  else {
    m_moms0[2] = m_pvv_zero[2][0]->Momentum();
    m_moms0[1] = m_pvv_zero[2][1]->Momentum();
  }
  // m_moms1 - project multiphoton state onto one photon phase space
  // do reconstruction procedure again pretending only one photon was generated

  // not necessary if only one photon
  if (pvv_one[4].size() == 1) {
    m_moms1[0][0] = pvv_one[1][0]->Momentum();
    if (m_switch == false) {
      m_moms1[0][1] = pvv_one[2][0]->Momentum();
      m_moms1[0][2] = pvv_one[2][1]->Momentum();
    }
    else {
      m_moms1[0][2] = pvv_one[2][0]->Momentum();
      m_moms1[0][1] = pvv_one[2][1]->Momentum();
    }
    m_moms1[0][3] = pvv_one[4][0]->Momentum();
  }
  else {
    Dipole_FF::DefineDipole();
    BoostOriginalPVVToMultipoleCMS();
    for (unsigned int i=0; i<pvv_one[4].size(); i++) {
      m_softphotons.push_back(pvv_one[4][i]);
      m_K = CalculateMomentumSum(m_softphotons);
      CorrectMomenta();
      if (m_switch == false) {
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
}

double Vector_To_Fermion_Fermion::Smod(unsigned int kk) {
  m_moms = m_moms1[kk];
  Vec4D k   = m_moms[3];
  Vec4D pi  = m_moms[1];
  Vec4D pj  = m_moms[2];
  double Zi = m_flavs[1].Charge();
  double Zj = m_flavs[2].Charge();
  int    ti = +1;
  int    tj = +1;
  return m_alpha/(4.*M_PI*M_PI)*Zi*Zj*ti*tj*(pi/(pi*k)-pj/(pj*k)).Abs2();
}

Complex Vector_To_Fermion_Fermion::InfraredSubtractedME_0_0() {
  m_moms = m_moms0;
  Vec4C epsV = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  return  XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cL,m_cR);
}

Complex Vector_To_Fermion_Fermion::InfraredSubtractedME_0_1() {
  return 0.;
  m_moms = m_moms0;
  Vec4C epsV = Polarization_Vector(m_moms[0])[m_spins[0]];
  XYZFunc XYZ(3,m_moms,m_flavs,false);
  double s((m_moms[1]+m_moms[2]).Abs2());
  double p1p2(m_moms[1]*m_moms[2]);
  double m(0.5*(m_masses[1]+m_masses[2]));
  double m2(sqr(m));
  double mu2(s);
  Complex term1(0.,0.), term2(0.,0.), term3(0.,0.), term4(0.,0.);
  // ~ u_1 \Gamma_V u_2
  term1 = XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cL,m_cR);
  term1 *=((p1p2+0.5*m2)*(C_11(m2,m2,s,0.,m2,m2,mu2)+C_12(m2,m2,s,0.,m2,m2,mu2))
           +(D-2.)/4.*(C_21(m2,m2,s,0.,m2,m2,mu2)+C_22(m2,m2,s,0.,m2,m2,mu2))
           +(p1p2+(D-4.)/2.)*C_23(m2,m2,s,0.,m2,m2,mu2)
           +0.25*sqr(D-2.)*C_24(m2,m2,s,0.,m2,m2,mu2)
           +0.25*B_0(s,m2,m2,mu2)
           -0.5*B_0(m2,0.,m2,mu2)
           +(D-1.)/4.*B_0(0.,m2,m2,mu2)).Finite();
  // ~ u_1 \tilde\Gamma_V u_2
  term2 = XYZ.X(1,m_spins[1],epsV,2,m_spins[2],m_cR,m_cL);
  term2 *= 0.5*m2*(C_11(m2,m2,s,0.,m2,m2,mu2)
                    +C_12(m2,m2,s,0.,m2,m2,mu2)).Finite();
  // ~ u_1 (LL+RR) u_2
  term3 = XYZ.Y(1,m_spins[1],2,m_spins[2],m_cL,m_cR);
  term3 *=(-m*(m_moms[2]*epsV)*(C_11(m2,m2,s,0.,m2,m2,mu2)
                                  +(D-2.)/2.*C_23(m2,m2,s,0.,m2,m2,mu2))
           -m*(m_moms[1]*epsV)*(D-2.)/2.*C_21(m2,m2,s,0.,m2,m2,mu2)).Finite();
  // ~ u_1 (LR+RL) u_2
  term4 = XYZ.Y(1,m_spins[1],2,m_spins[2],m_cR,m_cL);
  term4 *=(-m*(m_moms[1]*epsV)*(C_12(m2,m2,s,0.,m2,m2,mu2)
                                  +(D-2.)/2.*C_23(m2,m2,s,0.,m2,m2,mu2))
           -m*(m_moms[2]*epsV)*(D-2.)/2.*C_22(m2,m2,s,0.,m2,m2,mu2)).Finite();
  return m_alpha/M_PI*(term1+term2+term3+term4);
}

Complex Vector_To_Fermion_Fermion::InfraredSubtractedME_0_2() {
  return 0.;
}

Complex Vector_To_Fermion_Fermion::InfraredSubtractedME_1_05(unsigned int i) {
  m_moms       = m_moms1[i];                // set to set of momenta to be used
  Vec4C epsV   = Polarization_Vector(m_moms[0])[m_spins[0]];
  Vec4C epsP   = conj(Polarization_Vector(m_moms[3])[m_spins[3]]);
  Vec4D pa     = m_moms[1]+m_moms[3];       // fermion propagator momenta
  Vec4D pb     = m_moms[2]+m_moms[3];
  double m     = (0.5*(m_masses[1]+m_masses[2]));// fermion mass/propagator pole
  m_moms[4]    = m_moms[5] = pa;            // enter those into m_moms
  m_moms[6]    = m_moms[7] = pb;
  m_flavs[4]   = m_flavs[6] = m_flavs[1];   // set to corresponding particle/antiparticle
  m_flavs[5]   = m_flavs[7] = m_flavs[2];
  XYZFunc XYZ(8,m_moms,m_flavs,false);
  Complex r1 = Complex(0.,0.);
  Complex r2 = Complex(0.,0.);
  Complex r3 = Complex(0.,0.);
  Complex r4 = Complex(0.,0.);
  for (unsigned int s=0; s<=1; s++) {
    r1 += XYZ.X(1,m_spins[1],epsP,4,s,1.,1.)
          *XYZ.X(4,s,epsV,2,m_spins[2],m_cL,m_cR);
    r2 += XYZ.X(1,m_spins[1],epsP,5,s,1.,1.)
          *XYZ.X(5,s,epsV,2,m_spins[2],m_cL,m_cR);
    r3 += XYZ.X(1,m_spins[1],epsV,6,s,m_cL,m_cR)
          *XYZ.X(6,s,epsP,2,m_spins[2],1.,1.);
    r4 += XYZ.X(1,m_spins[1],epsV,7,s,m_cL,m_cR)
          *XYZ.X(7,s,epsP,2,m_spins[2],1.,1.);
  }
  // add prefactors
  r1 *= m_e/(2.*(pa*pa-m*m))*(1+m/sqrt(pa*pa));
  r2 *= m_e/(2.*(pa*pa-m*m))*(1-m/sqrt(pa*pa));
  r3 *= -m_e/(2.*(pb*pb-m*m))*(1-m/sqrt(pb*pb));
  r4 *= -m_e/(2.*(pb*pb-m*m))*(1+m/sqrt(pb*pb));
  // erase intermediate entries from m_flavs
  m_flavs[4] = m_flavs[5] = m_flavs[6] = m_flavs[7] = Flavour(kf_none);
  return (r1+r2+r3+r4);
}

Complex Vector_To_Fermion_Fermion::InfraredSubtractedME_1_15(unsigned int i) {
  return 0.;
}

Complex Vector_To_Fermion_Fermion::InfraredSubtractedME_2_1(unsigned int i, unsigned int j) {
  return 0.;
}

double Vector_To_Fermion_Fermion::GetBeta_0_0() {
  double sum = 0.;
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

double Vector_To_Fermion_Fermion::GetBeta_0_1() {
  // limit mV >> mf
  return m_alpha/M_PI*(2.*log(m_M/(0.5*(m_masses[1]+m_masses[2])))+3./2.)
           *GetBeta_0_0();
}

double Vector_To_Fermion_Fermion::GetBeta_0_2() {
  return 0.;
}

double Vector_To_Fermion_Fermion::GetBeta_1_1(unsigned int a) {
  double sum = 0.;
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
  sum = 1./(16.*M_PI*M_PI*M_PI)*sum - Smod(a)*GetBeta_0_0();
  return sum;
}

double Vector_To_Fermion_Fermion::GetBeta_1_2(unsigned int i) {
  return 0.;
}

double Vector_To_Fermion_Fermion::GetBeta_2_2(unsigned int i, unsigned int j) {
  return 0.;
}

DECLARE_PHOTONS_ME_GETTER(Vector_To_Fermion_Fermion_Getter,
                          "Vector_To_Fermion_Fermion")

PHOTONS_ME_Base * Vector_To_Fermion_Fermion_Getter::operator()
(const Particle_Vector_Vector &pvv) const
{
  // same mass restriction can be lifted if M_0_1 is computed for general case
  // Z needs to be excluded due to its different L/R properties
  if ( (pvv.size() == 4) &&
       (pvv[0].size() == 0) &&
       (pvv[1].size() == 1) && pvv[1][0]->Flav().IsVector() &&
       (pvv[1][0]->Flav().Kfcode() != kf_Z) &&
       (pvv[2].size() == 2) && pvv[2][0]->Flav().IsFermion() &&
       (pvv[2][0]->Flav() == pvv[2][1]->Flav().Bar()) &&
       (pvv[3].size() == 0) )
    return new Vector_To_Fermion_Fermion(pvv);
  return NULL;
}
